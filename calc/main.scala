
object Main {

val holidayDataPath = "holiday.txt";
val meteorDataPath = "meteor.txt";

val wordHistoryPath = "../etc/word-history.txt";

val dataPath = "../data.txt";

val tokyoLng = 139.7 / Const.PI57;
val tokyoLat = 35.7 / Const.PI57;

System.err.println("Started calculating...");

val hcs = new Hcs(tokyoLng, tokyoLat);

val words = new Words();
words.loadHistory(wordHistoryPath);

//==============================================================================
// イベント計算
//==============================================================================

val sunsetTimesData: IndexedSeq[(Double, Double, Array[Double])] = { // time, tdb, bpnMatrix
  val altHor = -0.90 / Const.PI57;
  (0 until Period.period).map { d =>
    val time = Lib2.findCrossingBoundaryTime(altHor, true, false, Period.startTime + d + 16.0 / 24.0, 24 * 6, 4 * 6) { time =>
      Acs.Horizontal.calcPlanetAlt(time, JplData.Sun, hcs);
    }
    val tdb = TimeLib.mjdutcToTdb(time);
    val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
    (time, tdb, bpnMatrix);
  }
}

def sunsetTimes(day: Int): Double = sunsetTimesData(day)._1;

def calcPlanetOnSunsetTime(day: Int, targetPlanet: JplData.TargetPlanet): (Double, Double) = { // azi, alt
  val (time, tdb, bpnMatrix) = sunsetTimesData(day);
  val xyz = JplData.calcPlanetFromEarth(tdb, targetPlanet);
  val xyz2 = VectorLib.multiplyMV(bpnMatrix, xyz);
  hcs.trueEquatorialXyzToAziAltFromUtc(xyz2, time);
}

def calcRiseSetCulmination(targetPlanet: JplData.TargetPlanet, altHor: Double):
  IndexedSeq[(Double, Double, Double, Double)] = { // 出, 南中, 没, 南中高度
  val riseSet = MathLib.findMaxMinCrossingListContinuous(Period.startTime, Period.endTime, 0.25, 24 * 6) { time =>
    Acs.Horizontal.calcPlanetAlt(time, targetPlanet, hcs) - altHor;
  }.filter(_._2 % 2 == 0);
  var result: List[(Double, Double, Double, Double)] = Nil;
  (0 until (riseSet.size - 1)).foreach { i =>
    val (time0, flag0) = riseSet(i);
    if (flag0 == 0) {
      val (time2, flag2) = riseSet(i + 1);
      val time1 = Lib2.findCrossingBoundaryTime(0.0, false, false,
        time0, 24 * 6, ((time2 - time0) * 24 * 6).toInt) { time =>
        Acs.Horizontal.calcPlanetAzi(time, targetPlanet, hcs) - Const.PI;
      }
      val alt = Acs.Horizontal.calcPlanetAlt(time1, targetPlanet, hcs);
      result = (time0, time1, time2, alt) :: result;
    }
  }
  result.reverse.toIndexedSeq;
}

val moonRiseSetTimesData: IndexedSeq[(Double, Double, Double, Double, Int)] = {
  val altHor = +0.4 / Const.PI57; // だいたい
  val data: IndexedSeq[(Double, Double, Double, Double)] = calcRiseSetCulmination(JplData.Moon, altHor);
  val data2 = Array.fill[Int](data.size)(-1);
  (0 until 8).foreach { phase =>
    val q = phase.toDouble / 8 * Const.PI2;
    MathLib.findMaxMinListDiscrete(0, data.size, 7) { p =>
      Math.abs(MathLib.circleDiff1(calcMoonPhase(data(p)._2) / 28 * Const.PI2, q));
    }.foreach { case (p, flag) =>
      if (flag < 0) {
        data2(p) = phase;
      }
    }
  }
  (0 until data.size).map { i =>
    val t = data(i);
    (t._1, t._2, t._3, t._4, data2(i));
  }
}
val moonRiseSetTimesDataTouch: Array[Boolean] = new Array[Boolean](moonRiseSetTimesData.size);

def touchMoonRiseSetStr(time0: Double): Option[String] = {
  val p = MathLib.binarySearchBy(moonRiseSetTimesData)(t => t._1 - time0);
  def isNextDay(time: Double): Boolean = {
    (time0 + 9.0 / 24).toInt < (time + 9.0 / 24).toInt;
  }
  def timeStr0(time: Double): String = TimeLib.timeToTimeNaturalString(time);
  if (p < 0 || p >= moonRiseSetTimesData.size) {
    Some("ERROR");
  } else if (moonRiseSetTimesData(p)._1 - time0 > 8.0 / 24) {
    None;
  } else {
    moonRiseSetTimesDataTouch(p) = true;
    val data = moonRiseSetTimesData(p);
    val southStr = calcMoonPhaseString(data._2);
    val format = if (isNextDay(data._1)) {
      "翌日の月の出は%sごろ、南中は%sごろ(%s)、月の入りは%sごろ";
    } else if (isNextDay(data._2)) {
      "月の出は%sごろ、南中は翌日%sごろ(%s)、月の入りは翌日%sごろ";
    } else if (isNextDay(data._3)) {
      "月の出は%sごろ、南中は%sごろ(%s)、月の入りは翌日%sごろ";
    } else {
      "月の出は%sごろ、南中は%sごろ(%s)、月の入りは%sごろ";
    }
    Some(format.format(timeStr0(data._1), timeStr0(data._2), southStr, timeStr0(data._3)));
  }
}

def moonRiseSetForMeteorShower(time0: Double): String = {
  val sopt = touchMoonRiseSetStr(time0);
  if (sopt.isDefined) {
    return sopt.get;
  }
  def isNextDay(time: Double): Boolean = {
    (time0 + 9.0 / 24).toInt < (time + 9.0 / 24).toInt;
  }
  def timeStr0(time: Double): String = TimeLib.timeToTimeNaturalString(time);
  val p1 = MathLib.binarySearchBy(moonRiseSetTimesData)(t => t._1 - time0);
  val p2 = MathLib.binarySearchBy(moonRiseSetTimesData)(t => t._3 - time0);
  val moonRiseTime = moonRiseSetTimesData(p1)._1;
  val moonSetTime = moonRiseSetTimesData(p2)._3;
  if (moonRiseTime < moonSetTime) {
    if (isNextDay(moonRiseTime)) {
      "翌日の月の出は%sごろ、月の入りは%sごろ".format(
        timeStr0(moonRiseTime), timeStr0(moonSetTime));
    } else if (isNextDay(moonSetTime)) {
      "月の出は%sごろ、月の入りは翌日の%sごろ".format(
        timeStr0(moonRiseTime), timeStr0(moonSetTime));
    } else {
      "月の出は%sごろ、月の入りは%sごろ".format(
        timeStr0(moonRiseTime), timeStr0(moonSetTime));
    }
  } else {
    if (isNextDay(moonSetTime)) {
      "翌日の月の入りは%sごろ、月の出は%sごろ".format(
        timeStr0(moonSetTime), timeStr0(moonRiseTime));
    } else if (isNextDay(moonRiseTime)) {
      "月の入りは%sごろ、月の出は翌日の%sごろ".format(
        timeStr0(moonSetTime), timeStr0(moonRiseTime));
    } else {
      "月の入りは%sごろ、月の出は%sごろ".format(
        timeStr0(moonSetTime), timeStr0(moonRiseTime));
    }
  }
}

def tweetMoonRiseSet(): Unit = {
  (0 until moonRiseSetTimesData.size).foreach { p =>
    if (!moonRiseSetTimesDataTouch(p)) {
      val time = moonRiseSetTimesData(p)._1;
      val time2 = {
        val d = time - time.toInt;
        if (d < 14.5 / 24) { // 9時～23時30分
          time;
        } else if (d < 21.0 / 24) { // 23時30分～6時
          time - d + 14.5 / 24; // 23時30分
        } else { // 6時～9時
          time;
        }
      }
      val riseset = touchMoonRiseSetStr(time2).get;
      tweetsManager.putTweet(time2, riseset);
    }
  }
}

def moonPhaseNaturalString(moonPhase: Double): String = {
  if (moonPhase < 1.75) {
    "新月後の細い月";
  } else if (moonPhase <= 2.75) {
    "三日月";
  } else if (moonPhase < 6.0) {
    "";
  } else if (moonPhase < 7.0) {
    "もうすぐ上弦の月";
  } else if (moonPhase == 7.0) {
    "上弦の月";
  } else if (moonPhase <= 8.0) {
    "上弦の月を過ぎた月";
  } else if (moonPhase < 13.0) {
    "";
  } else if (moonPhase < 14.0) {
    "もうすぐ満月";
  } else if (moonPhase == 14.0) {
    "満月";
  } else if (moonPhase <= 15.0) {
    "満月を過ぎた月";
  } else {
    "";
  }
}

def calcMoonLng(time: Double): Double = {
  val utc = time;
  val tdb = TimeLib.mjdutcToTdb(utc);
  val sun = JplData.calcPlanetFromEarth(tdb, JplData.Sun);
  val moon = JplData.calcPlanetFromEarth(tdb, JplData.Moon);
  val bpnMatrix = Bpn.icrsToTrueEclipticMatrix(tdb);
  val sun2 = VectorLib.multiplyMV(bpnMatrix, sun);
  val moon2 = VectorLib.multiplyMV(bpnMatrix, moon);
  val sunLng = VectorLib.xyzToLng(sun2);
  val moonLng = VectorLib.xyzToLng(moon2);
  VectorLib.calcLngDiff(moonLng, sunLng);
}
def calcMoonPhase(time: Double): Double = {
  calcMoonLng(time) / Const.PI2 * 28.0;
}
def calcMoonPhaseRatio(time: Double): Double = {
  0.5 * (1.0 - Math.cos(calcMoonLng(time)));
}
def calcMoonPhaseString(time: Double): String = {
  val p = MathLib.binarySearchBy(moonRiseSetTimesData)(t => t._2 - time);
  val phaseTerm = if (p >= 0 && p < moonRiseSetTimesData.size && moonRiseSetTimesData(p)._2 == time) moonRiseSetTimesData(p)._5 else -1;

  val moonPhase = calcMoonPhase(time);
  val moonPhaseRatio = 100.0 * calcMoonPhaseRatio(time);
  val moonDayCountStr = {
    val moonPhaseTerms2 = moonPhaseTerms.filter(_._2 == 0);
    val idx = MathLib.binarySearchBy(moonPhaseTerms2)(t => t._1 - time);
    if (idx <= 0) {
      "ERROR";
    } else {
      val newMoonTime = moonPhaseTerms2(idx - 1)._1;
      "月齢%.1f".format(time - newMoonTime);
    }
  }
  val phaseStr = if (moonPhase > 27.6 || moonPhase > 27.0 && phaseTerm == 0) {
    "、新月\uD83C\uDF11の直前";
  } else if (moonPhase < 0.4 || moonPhase < 1.0 && phaseTerm == 0) {
    "、新月\uD83C\uDF11の直後";
  } else if (moonPhase < 1.75) {
    "、新月後の細い月";
  } else if (moonPhase < 2.75) {
    "、三日月";
  } else if (moonPhase > 3.1 && moonPhase < 3.9 || moonPhase > 2.5 && moonPhase < 4.5 && phaseTerm == 1) {
    "、新月と上弦の中間\uD83C\uDF12";
  } else if (moonPhase > 6.6 && moonPhase < 7.4 || moonPhase > 6.0 && moonPhase < 8.0 && phaseTerm == 2) {
    "、上弦の月\uD83C\uDF13";
  } else if (moonPhase > 6.0 && moonPhase <= 6.6) {
    "、もうすぐ上弦の月";
  } else if (moonPhase > 7.4 && moonPhase < 8.0) {
    "、上弦を過ぎた月";
  } else if (moonPhase > 10.1 && moonPhase < 10.9 || moonPhase > 9.5 && moonPhase < 11.5 && phaseTerm == 3) {
    "、上弦と満月の中間\uD83C\uDF14";
  } else if (moonPhase > 13.6 && moonPhase < 14.0 || moonPhase > 13.0 && moonPhase < 14.0 && phaseTerm == 4) {
    "、満月\uD83C\uDF15の直前";
  } else if (moonPhase >= 14.0 && moonPhase < 14.4 || moonPhase >= 14.0 && moonPhase < 15.0 && phaseTerm == 4) {
    "、満月\uD83C\uDF15の直後";
  } else if (moonPhase > 13.0 && moonPhase <= 13.6) {
    "、もうすぐ満月";
  } else if (moonPhase > 14.4 && moonPhase < 15.0) {
    "、満月を過ぎた月";
  } else if (moonPhase > 17.1 && moonPhase < 17.9 || moonPhase > 16.5 && moonPhase < 18.5 && phaseTerm == 5) {
    "、満月と下弦の中間\uD83C\uDF16";
  } else if (moonPhase > 20.6 && moonPhase < 21.4 || moonPhase > 20.0 && moonPhase < 22.0 && phaseTerm == 6) {
    "、下弦の月\uD83C\uDF17";
  } else if (moonPhase > 24.1 && moonPhase < 24.9 || moonPhase > 23.5 && moonPhase < 25.5 && phaseTerm == 7) {
    "、下弦と新月の中間\uD83C\uDF18";
  } else {
    "";
  }
  "月相%.1f/28、輝面比%.1f%%、%s%s".format(moonPhase, moonPhaseRatio, moonDayCountStr, phaseStr);
}
val moonPhaseTerms: IndexedSeq[(Double, Int)] = { // time, term
  MathLib.findCyclicPhaseListContinuous(8, Period.startTime, Period.endTime, 3.0, 24)(calcMoonLng);
}

val innerPlanets = IndexedSeq(
  ("金星", JplData.Venus, 5),
  ("水星", JplData.Mercury, 3)
);
val innerPlanetsSunsetAziAltList: IndexedSeq[IndexedSeq[(Double, Double)]] = innerPlanets.map { planet =>
  (0 until Period.period).map { day =>
    calcPlanetOnSunsetTime(day, planet._2);
  }
}

val outerPlanets = IndexedSeq(
  ("火星", JplData.Mars, 2),
  ("木星", JplData.Jupiter, 4),
  ("土星", JplData.Saturn, 6),
);
val outerPlanetsNightAziAltList: IndexedSeq[IndexedSeq[(Double, Double)]] = outerPlanets.map { planet =>
  val period2 = if (planet._1 == "火星") {
    Period.period + 700;
  } else {
    Period.period + 300;
  }
  (0 until period2).map { day =>
    val time = Period.startTime + day + 21.0 / 24.0;
    Acs.Horizontal.calcPlanetAziAlt(time, planet._2, hcs);
  }
}

//==============================================================================

def isLunchTime(time: Double): Boolean = {
  val day = (time - Period.startTime).toInt;
  val s = sunsetTimes(day);
  //time - time.toInt >= 3.0 / 24 && time - time.toInt < 4.0 / 24; // 12時～13時
  time - time.toInt >= 1.0 / 24 && time - time.toInt < 7.0 / 24; // 10時～16時
}

def isNightTime0(time: Double): Boolean = {
  val day = (time - Period.startTime).toInt;
  val s = sunsetTimes(day);
  //time - time.toInt >= 5.0 / 24 && time - time.toInt < 15.0 / 24;
  time > s && time - time.toInt < 14.0 / 24;
}

def isNightTime2(time: Double): Boolean = {
  val day = (time - Period.startTime).toInt;
  val s = sunsetTimes(day);
  time >= s && time - time.toInt < 16.0 / 24;
}

def isNightTime3(time: Double): Boolean = {
  val day = (time - Period.startTime).toInt;
  val s = sunsetTimes(day);
  time - time.toInt >= 11.0 / 24 && time - time.toInt < 14.1 / 24;
}

def getConjunctionTweetTime(time: Double, xyz: Array[Double]): Option[Double] = {
  val altThres0 = 10 / Const.PI57;
  val sun_distanceThres = 20.0 / Const.PI57;
  val utc = time;
  val ut1 = time; // 近似的
  val tdb = TimeLib.mjdutcToTdb(utc);
  val sun_xyz = JplData.calcPlanetFromEarth(tdb, JplData.Sun);
  val sun_distance = VectorLib.angularDistance(sun_xyz, xyz);
  if (sun_distance < sun_distanceThres) {
    None;
  } else {
    val time21 = time.toInt + 0.5;
    val tdb21 = TimeLib.mjdutcToTdb(time21);
    val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb21);
    val xyz2 = VectorLib.multiplyMV(bpnMatrix, xyz);
    val (azi, alt) = hcs.trueEquatorialXyzToAziAltFromUtc(xyz2, time21);
    if (alt >= altThres0) {
      // 21時に見える
      Some(time21);
    } else {
      val (azi, alt) = hcs.trueEquatorialXyzToAziAltFromUtc(xyz2, time21 + 2.0 / 24);
      if (alt >= altThres0) {
        // 23時に見える
        Some(time21 + 2.0 / 24);
      } else {
        val day = (time - Period.startTime).toInt;
        val sunsetTime = sunsetTimes(day);
        val (azi, alt) = hcs.trueEquatorialXyzToAziAltFromUtc(xyz2, sunsetTime);
        if (alt >= altThres0) {
          // 日没時に見える
          Some(sunsetTime);
        } else {
          None;
        }
      }
    }
  }
}

//==============================================================================
// ツイート管理
//==============================================================================

val tweetsManager = new Tweets();

//==============================================================================
// 祝日・記念日
//==============================================================================

{
  val yearStart = TimeLib.timeToDateTimeString(Period.startTime).substring(0, 4).toInt;
  val yearEnd = TimeLib.timeToDateTimeString(Period.endTime).substring(0, 4).toInt;
  val source = scala.io.Source.fromFile(holidayDataPath);
  source.getLines.foreach { line =>
    if (line != "" && !line.startsWith("#")) {
      val cols = line.split(" ", 2);
      if (cols(0).startsWith("YYYY-")) {
        val msg = cols(1);
        (yearStart to yearEnd).foreach { year =>
          val timeStr = year.toString + cols(0).substring(4);
          val time = TimeLib.stringToModifiedJulianDay(timeStr + ":00+09:00");
          tweetsManager.putTweet(time, msg);
        }
      } else {
        val time = TimeLib.stringToModifiedJulianDay(cols(0) + ":00+09:00");
        val msg = cols(1);
        tweetsManager.putTweet(time, msg);
      }
    }
  }
  source.close();
}

//==============================================================================
// 特定の時間に発生する天文の現象
//==============================================================================

//--------------------------------------
// 24節気
//--------------------------------------

SolarTerms.tweets.foreach { tw =>
  tweetsManager.putTweet(tw);
}

//--------------------------------------
// 流星群
//--------------------------------------

{
  val meteorShowerData: (IndexedSeq[(Double, String)], IndexedSeq[(Double, String)]) = {
    import Ordering.Double.IeeeOrdering;
    var data1: List[(Double, String)] = Nil;
    var data2: List[(Double, String)] = Nil;
    val source = scala.io.Source.fromFile(meteorDataPath);
    source.getLines.foreach { line =>
      if (line != "" && !line.startsWith("#")) {
        val cols = line.split(" ", 2);
        if (cols(0).indexOf("T") == 10) {
          val time = TimeLib.stringToModifiedJulianDay(cols(0) + ":00+09:00");
          data2 = (time, cols(1)) :: data2;
        } else {
          val lng = cols(0).toDouble / Const.PI57;
          data1 = (lng, cols(1)) :: data1;
        }
      }
    }
    source.close();
    (data1.reverse.sortBy(_._1).toIndexedSeq, data2.reverse.toIndexedSeq);
  }

  var day: Int = 0;
  var index: Int = {
    val time = Period.startTime + day + 0.5;
    val lng = Acs.Ecliptic2000.calcPlanetLng(time, JplData.Sun);
    val index = meteorShowerData._1.indexWhere(_._1 > lng);
    if (index < 0) {
      0;
    } else {
      index;
    }
  }
  while (day < Period.period) {
    val time = Period.startTime + day + 0.5;
    val lng = Acs.Ecliptic2000.calcPlanetLng(time, JplData.Sun);
    val d = MathLib.circleAdd(lng, -meteorShowerData._1(index)._1);
    if (d >= 0) {
      val time1 = time - 1.0 + 20.0 / 60 / 24;
      val name = meteorShowerData._1(index)._2;
      val msg = "%sは今夜が極大。%s".format(name, moonRiseSetForMeteorShower(time1));
      tweetsManager.putTweet(time1, msg);
      index += 1;
      if (index == meteorShowerData._1.size) {
        index = 0;
      }
    } else {
      day += 1;
    }
  }

  meteorShowerData._2.foreach { case (time, msg) =>
    tweetsManager.putTweet(time, msg);
  }
}

//--------------------------------------
// 月相
//--------------------------------------

{
  case class MoonPhaseTermTweetContent(rawTime: Double, term: Int, distanceFlag: Int,
    cons: Option[(String, List[String])]) extends TweetContent {

    private[this] val riseset: String = {
      touchMoonRiseSetStr(rawTime) match {
        case None => "";
        case Some(s) => "。" + s;
      }
    }

    def time: Double = TimeLib.floor(rawTime, 24) + 1.0 / (24 * 4);
    def message: String = {
      val msg = if (term == 4 && distanceFlag < 0) {
        "満月\uD83C\uDF15。月が地球に近く、もっとも大きい満月です。スーパームーンとも呼ばれます。月相14/28";
      } else if (term == 4 && distanceFlag > 0) {
        "満月\uD83C\uDF15。月が地球から遠く、もっとも小さい満月です。月相14/28";
      } else {
        termStrs(term);
      }
      cons match {
        case Some((cons, hashtags)) => "%s。%sにいます%s".format(msg, cons, riseset);
        case None => "%s%s".format(msg, riseset);
      }
    }
    def hashtags: List[String] = cons.map(_._2).getOrElse(Nil);
    def starNames: List[String] = List("月");

    private[this] val termStrs = IndexedSeq(
      "新月\uD83C\uDF11。月相0/28",
      "月相3.5/28。新月と上弦の中間です\uD83C\uDF12",
      "上弦の月\uD83C\uDF13。月相7/28",
      "月相10.5/28。上弦と満月の中間です\uD83C\uDF14",
      "満月\uD83C\uDF15。月相14/28",
      "月相17.5/28。満月と下弦の中間です\uD83C\uDF16",
      "下弦の月\uD83C\uDF17。月相21/28",
      "月相24.5/28。下弦と新月の中間です\uD83C\uDF18",
    );
  }

  val altThres0 = 10 / Const.PI57;

  def calcMoonConstellation(time: Double): Option[(String, List[String])] = {
    if (isNightTime2(time)) {
      val tdb = TimeLib.mjdutcToTdb(time);
      val xyz = JplData.calcPlanetFromEarth(tdb, JplData.Moon);
      val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
      val xyz2 = VectorLib.multiplyMV(bpnMatrix, xyz);
      val (azi, alt) = hcs.trueEquatorialXyzToAziAltFromUtc(xyz2, time);
      if (alt >= altThres0) {
        val (cons, hashtags) = Constellations.icrsToConstellation(xyz);
        Some((cons, hashtags));
      } else {
        None;
      }
    } else {
      None;
    }
  }

  {
    val fullMoons = moonPhaseTerms.filter(_._2 == 4).map(_._1);
    val fullMoonsDistanceMaxMinUpDownFlags = MathLib.getMaxMinUpDownFlagListDiscrete(0, fullMoons.size, 2) { idx =>
      val utc = fullMoons(idx);
      val tdb = TimeLib.mjdutcToTdb(utc);
      val moon = JplData.calcPlanetFromEarth(tdb, JplData.Moon);
      VectorLib.distance(moon);
    }
    fullMoonsDistanceMaxMinUpDownFlags.zipWithIndex.map { case (flag, idx) =>
      val time = fullMoons(idx);
      if (flag == 1) {
        MoonPhaseTermTweetContent(time, 4, +1, calcMoonConstellation(time));
      } else if (flag == 3) {
        MoonPhaseTermTweetContent(time, 4, -1, calcMoonConstellation(time));
      } else {
        MoonPhaseTermTweetContent(time, 4, 0, calcMoonConstellation(time));
      }
    }.foreach(tweetsManager.putTweet);

    val midAltData = moonRiseSetTimesData.filter(_._5 == 4).map(t => (t._2, t._4));
    (0 until midAltData.size).foreach { i =>
      val (time, alt) = midAltData(i);
      if (i > 0 && i < midAltData.size - 1) {
        if (midAltData(i - 1)._2 <= alt && midAltData(i + 1)._2 < alt) {
          tweetsManager.putTweet(time, "満月が南中(高度%.0f°)。冬の満月は空高く上り、今日はこの冬でもっとも天頂に近い満月です".format(alt * Const.PI57));
        } else if (midAltData(i - 1)._2 >= alt && midAltData(i + 1)._2 > alt) {
          tweetsManager.putTweet(time, "満月が南中(高度%.0f°)。夏の満月は空低く、今日はこの夏でもっとも低い満月です".format(alt * Const.PI57));
        } else {
          tweetsManager.putTweet(time, "満月が南中(高度%.0f°)".format(alt * Const.PI57));
        }
      } else {
        tweetsManager.putTweet(time, "満月が南中(高度%.0f°)".format(alt * Const.PI57));
      }
    }
  }

  moonPhaseTerms.filter(_._2 != 4).filter(_._2 % 2 == 0).map { case (time, term) =>
    MoonPhaseTermTweetContent(time, term, 0, calcMoonConstellation(time));
  }.foreach(tweetsManager.putTweet);
}

//--------------------------------------
// 近日点・遠日点
//--------------------------------------

MathLib.findMaxMinListContinuous(Period.startTime, Period.endTime, 30, 24) { time =>
  val utc = time;
  val tdb = TimeLib.mjdutcToTdb(utc);
  val sun = JplData.calcPlanetFromEarth(tdb, JplData.Sun);
  VectorLib.distance(sun);
}.foreach { case (time, flag) =>
  val s = if (flag < 0) "近日点" else "遠日点";
  tweetsManager.putTweet(TimeLib.round(time, 24) - 1.0 / (24 * 4), "地球が%s通過".format(s));
}

//--------------------------------------
// 惑星の天象
//--------------------------------------

case class PlanetAstronomyTweetContent(time: Double, message: String, planetName: String) extends TweetContent {
  def hashtags: List[String] = List(planetName);
  def starNames: List[String] = List(planetName);
}

{
  def calcInnerPlanetLngEc(time: Double, targetPlanet: JplData.TargetPlanet): Double = {
    val utc = time;
    val tdb = TimeLib.mjdutcToTdb(utc);
    val sun = JplData.calcPlanetFromEarth(tdb, JplData.Sun);
    val planet = JplData.calcPlanetFromEarth(tdb, targetPlanet);
    val bpnMatrix = Bpn.icrsToTrueEclipticMatrix(tdb);
    val sun2 = VectorLib.multiplyMV(bpnMatrix, sun);
    val planet2 = VectorLib.multiplyMV(bpnMatrix, planet);
    val sunLng = VectorLib.xyzToLng(sun2);
    val planetLng = VectorLib.xyzToLng(planet2);
    val d1 = VectorLib.calcLngDiff(planetLng, sunLng);
    if (d1 >= Const.PI) d1 - Const.PI2 else d1;
  }

  {
    val (planetName, targetPlanet) = ("水星", JplData.Mercury);
    MathLib.findMaxMinCrossingListContinuous(Period.startTime, Period.endTime, 10.0, 24) { time =>
      calcInnerPlanetLngEc(time, targetPlanet);
    }.map { case (time, term) =>
      if (term == 0) {
        tweetsManager.putTweet(PlanetAstronomyTweetContent(TimeLib.floor(time, 24) + 1.0 / (24 * 4),
          "水星が外合(黄経基準)", planetName));
      } else if (term == 1) {
        tweetsManager.putTweet(PlanetAstronomyTweetContent(TimeLib.round(time, 24) - 1.0 / (24 * 4),
          "水星が東方最大離角(黄経基準)\uD83C\uDF13", planetName));
      } else if (term == 2) {
        tweetsManager.putTweet(PlanetAstronomyTweetContent(TimeLib.floor(time, 24) + 1.0 / (24 * 4),
          "水星が内合(黄経基準)", planetName));
      } else {
        tweetsManager.putTweet(PlanetAstronomyTweetContent(TimeLib.round(time, 24) - 1.0 / (24 * 4),
          "水星が西方最大離角(黄経基準)\uD83C\uDF17", planetName));
      }
    }
  }

  {
    val (planetName, targetPlanet) = ("金星", JplData.Venus);
    MathLib.findMaxMinCrossingListContinuous(Period.startTime, Period.endTime, 10.0, 24) { time =>
      calcInnerPlanetLngEc(time, targetPlanet);
    }.map { case (time, term) =>
      if (term == 0) {
        tweetsManager.putTweet(PlanetAstronomyTweetContent(TimeLib.floor(time, 24) + 1.0 / (24 * 4),
          "金星が外合(黄経基準)。数か月後に夕方の西の空に現れます", planetName));
      } else if (term == 1) {
        tweetsManager.putTweet(PlanetAstronomyTweetContent(TimeLib.round(time, 24) - 1.0 / (24 * 4),
          "金星が東方最大離角(黄経基準)\uD83C\uDF13。宵の明星として夕方に西の空にいます", planetName));
      } else if (term == 2) {
        tweetsManager.putTweet(PlanetAstronomyTweetContent(TimeLib.floor(time, 24) + 1.0 / (24 * 4),
          "金星が内合(黄経基準)。数週間後に明け方の東の空に現れます", planetName));
      } else {
        tweetsManager.putTweet(PlanetAstronomyTweetContent(TimeLib.round(time, 24) - 1.0 / (24 * 4),
          "金星が西方最大離角(黄経基準)\uD83C\uDF17。明けの明星として明け方に東の空にいます", planetName));
      }
    }
  }
}

{
  def calcOuterPlanetLngEq(time: Double, targetPlanet: JplData.TargetPlanet): Double = {
    val utc = time;
    val tdb = TimeLib.mjdutcToTdb(utc);
    val sun = JplData.calcPlanetFromEarth(tdb, JplData.Sun);
    val planet = JplData.calcPlanetFromEarth(tdb, targetPlanet);
    val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
    val sun2 = VectorLib.multiplyMV(bpnMatrix, sun);
    val planet2 = VectorLib.multiplyMV(bpnMatrix, planet);
    val sunLng = VectorLib.xyzToLng(sun2);
    val planetLng = VectorLib.xyzToLng(planet2);
    val d1 = VectorLib.calcLngDiff(planetLng, sunLng);
    d1;
  }

  val planetPhases1: List[(Double, String, String, Boolean, Option[Array[Double]], Int)] = outerPlanets.zipWithIndex.toList.
  flatMap { case ((planetName, targetPlanet, _), pi) =>
    MathLib.findCyclicPhaseListContinuous(4, Period.startTime, Period.endTime, 30, 24) { time =>
      Const.PI2 - calcOuterPlanetLngEq(time, targetPlanet);
    }.map { case (time, term) =>
      val utc = time;
      val tdb = TimeLib.mjdutcToTdb(utc);
      val xyz = JplData.calcPlanetFromEarth(tdb, targetPlanet);
      if (term == 0) {
        (time, planetName, "合(赤経基準)", false, None, pi);
      } else if (term == 1) {
        (time, planetName, "西矩(赤経基準)", false, Some(xyz), pi);
      } else if (term == 2) {
        (time, planetName, "衝(赤経基準)", false, Some(xyz), pi);
      } else {
        (time, planetName, "東矩(赤経基準)", false, Some(xyz), pi);
      }
    }
  }

  val altThres = 15 / Const.PI57;

  planetPhases1.foreach { case (time, planetName, content, timeFlag, xyzOpt, pi) =>
    val time2 = if (timeFlag) {
      TimeLib.round(time, 24) - 1.0 / (24 * 4);
    } else {
      TimeLib.floor(time, 24) + 1.0 / (24 * 4);
    }
    val nextSeasonMsg = {
      val day = (time - Period.startTime).toInt;
      val (azi, alt) = outerPlanetsNightAziAltList(pi)(day);
      if (alt < altThres) {
        val (p1, p2) = calcNextSeason(outerPlanetsNightAziAltList(pi), day);
        "。" + NightNextPlanetTweetContent(time2, planetName, p1, p2).message;
      } else {
        "";
      }
    }
    val tweet = xyzOpt match {
      case None =>
        PlanetAstronomyTweetContent(time2, "%sが%s%s".format(planetName, content, nextSeasonMsg), planetName);
      case Some(xyz) =>
        val (cons, hashtags) = Constellations.icrsToConstellation(xyz);
        PlanetAstronomyTweetContent(time2, "%sが%s。%sにいます%s".format(planetName, content, cons, nextSeasonMsg) +
          hashtags.map(" #" + _).mkString, planetName);
    }
    tweetsManager.putTweet(tweet);
  }
}

case class SunsetNextPlanetTweetContent(time: Double, planetName: String,
  nextDay1: Int, nextDay2: Int) extends TweetContent {
  def message: String = {
    if (nextDay1 < 0) {
      "%sが次に日没時に見えやすくなるのは、ERRORです".format(planetName);
    } else if (nextDay2 < 0) {
      "%sが次に日没時に見えやすくなるのは、%sからERRORです".format(planetName,
        TimeLib.monthString(time, nextDay1, nextDay1));
    } else {
      "%sが次に日没時に見えやすくなるのは、%sです".format(planetName,
        TimeLib.monthString(time, nextDay1, nextDay2));
    }
  }
  def hashtags: List[String] = List(planetName);
  def starNames: List[String] = List(planetName);
}
case class NightNextPlanetTweetContent(time: Double, planetName: String,
  nextDay1: Int, nextDay2: Int) extends TweetContent {
  def message: String = {
    if (nextDay1 < 0) {
      "%sが次に夜見えやすくなるのは、ERRORです(21時基準)".format(planetName);
    } else if (nextDay2 < 0) {
      "%sが次に夜見えやすくなるのは、%sからERRORです(21時基準)".format(planetName,
        TimeLib.monthString(time, nextDay1, nextDay1));
    } else {
      "%sが次に夜見えやすくなるのは、%sです(21時基準)".format(planetName,
        TimeLib.monthString(time, nextDay1, nextDay2));
    }
  }
  def hashtags: List[String] = List(planetName);
  def starNames: List[String] = List(planetName);
}

//==============================================================================
// 日没時の西の空
//==============================================================================

sealed trait OnSunsetTweetContent extends TweetContent {
  def day: Int;
  def message: String;
  def message2: String;
  def message3: String;
  def hashtags: List[String];
  def starNames: List[String];

  def time: Double = sunsetTimes(day);
}

case class MultiSunsetTweetContent(day: Int, tc: List[OnSunsetTweetContent]) extends TweetContent {
  def time: Double = sunsetTimes(day);
  def message: String = tc.head.message2 + "。" + tc.tail.map(_.message3).mkString("。");
  def hashtags: List[String] = tc.flatMap(_.hashtags);
  def starNames: List[String] = tc.flatMap(_.starNames);
}

case class SunsetMoonTweetContent(day: Int, azi: Double, alt: Double) extends OnSunsetTweetContent {
  private[this] val azi360: Int = (azi * Const.PI57 + 0.5).toInt;
  private[this] val alt360: Int = (alt * Const.PI57 + 0.5).toInt;
  private[this] val moonPhase: Double = calcMoonPhase(sunsetTimes(day));
  private[this] val moonStr: String = {
    val s = moonPhaseNaturalString(moonPhase);
    if (s == "") {
      s;
    } else {
      "。" + s + "です";
    }
  }
  def message: String = "月は日没時に西の空高度約%d°%s".format(alt360, moonStr);
  def message2: String = "月は日没時に西の空高度約%d°にいます%s".format(alt360, moonStr);
  def message3: String = "月は約%d°にいます%s".format(alt360, moonStr);
  def hashtags: List[String] = Nil;
  def starNames: List[String] = List("月");
}
case class SunsetPlanetTweetContent(day: Int, planetName: String,
  azi: Double, alt: Double, isIncreasing: Boolean, isDecreasing: Boolean, isMax: Boolean) extends OnSunsetTweetContent {
  def alt360: Int = (alt * Const.PI57 + 0.5).toInt;
  def message: String = if (isMax) {
    "%sは日没時最大高度で西の空高度約%d°".format(planetName, alt360);
  } else if (isIncreasing) {
    "%sは日没時の高度を徐々に上げ、西の空高度約%d°にいます".format(planetName, alt360);
  } else if (isDecreasing) {
    "%sは日没時の高度を徐々に下げ、西の空高度約%d°にいます".format(planetName, alt360);
  } else {
    "%sは西の空高度約%d°にいます".format(planetName, alt360);
  }
  def message2: String = if (isMax) {
    "%sは日没時最大高度で西の空高度約%d°です".format(planetName, alt360);
  } else if (isIncreasing) {
    "%sは日没時の高度を徐々に上げ、西の空高度約%d°にいます".format(planetName, alt360);
  } else if (isDecreasing) {
    "%sは日没時の高度を徐々に下げ、西の空高度約%d°にいます".format(planetName, alt360);
  } else {
    "%sは西の空高度約%d°にいます".format(planetName, alt360);
  }
  def message3: String = if (isMax) {
    "%sは日没時最大高度で西の空高度約%d°です".format(planetName, alt360);
  } else if (isIncreasing) {
    "%sは日没時の高度を徐々に上げ、約%d°にいます".format(planetName, alt360);
  } else if (isDecreasing) {
    "%sは日没時の高度を徐々に下げ、約%d°にいます".format(planetName, alt360);
  } else {
    "%sは約%d°にいます".format(planetName, alt360);
  }
  def hashtags: List[String] = List(planetName);
  def starNames: List[String] = List(planetName);
}
case class SunsetStarTweetContent(day: Int, starName: String,
  azi: Double, alt: Double) extends OnSunsetTweetContent {
  def alt360: Int = (alt * Const.PI57 + 0.5).toInt;
  def message: String = "%sは西の空高度約%d°にいます".format(starName, alt360);
  def message2: String = "%sは西の空高度約%d°にいます".format(starName, alt360);
  def message3: String = "%sは約%d°にいます".format(starName, alt360);
  def hashtags: List[String] = List(starName);
  def starNames: List[String] = List(starName);
}

case class SunsetTweetContent(day: Int, flag: Int) extends OnSunsetTweetContent {
  def message(level: Int): String = {
    if (flag == 1) {
      "日没はこのころが最も遅く、%sごろです".format(TimeLib.timeToTimeNaturalString(time));
    } else if (flag == 3) {
      "日没はこのころが最も早く、%sごろです".format(TimeLib.timeToTimeNaturalString(time));
    } else {
      if (day >= 7) {
        val timePrev = sunsetTimes(day - 7);
        val d = Math.round((time - (timePrev + 7.0)) * (24 * 60));
        if (d > 0) {
          "日没は%sごろで、この1週間で約%d分遅くなっています".format(TimeLib.timeToTimeNaturalString(time), d);
        } else if (d < 0) {
          "日没は%sごろで、この1週間で約%d分早くなっています".format(TimeLib.timeToTimeNaturalString(time), -d);
        } else {
          if (level == 2) {
            "日没は%sごろ".format(TimeLib.timeToTimeNaturalString(time));
          } else {
            "日没は%sごろです".format(TimeLib.timeToTimeNaturalString(time));
          }
        }
      } else {
        if (level == 2) {
          "日没は%sごろ".format(TimeLib.timeToTimeNaturalString(time));
        } else {
          "日没は%sごろです".format(TimeLib.timeToTimeNaturalString(time));
        }
      }
    }
  }
  def message: String = message(1);
  def message2: String = message(2);
  def message3: String = message2;
  def hashtags: List[String] = Nil;
  def starNames: List[String] = Nil;
}

case class SunsetFirstStarTweetContent(day: Int, msg: String, starNames: List[String]) extends OnSunsetTweetContent {
  def message: String = "#" + msg;
  def message2: String = message;
  def message3: String = message;
  def hashtags: List[String] = Nil;
}

// 日没時最大高度
{
  val planets = IndexedSeq(("金星", JplData.Venus), ("水星", JplData.Mercury));
  planets.foreach { p =>
    MathLib.findMaxMinListDiscrete(0, Period.period, 15) { day =>
      val (azi, alt) = calcPlanetOnSunsetTime(day, p._2);
      alt;
    }.foreach { case (day, flag) =>
      if (flag > 0) {
        val (azi, alt) = calcPlanetOnSunsetTime(day, p._2);
        tweetsManager.putTweet(SunsetPlanetTweetContent(day, p._1, azi, alt, false, false, true));
      }
    }
  }
}

// 新月直後の月
{
  val aziThres0 = 200 / Const.PI57;
  val aziThres1 = 315 / Const.PI57;
  val altThres0 = 10 / Const.PI57;

  (0 until moonPhaseTerms.size).foreach { i =>
    val (moonPhaseTime, term) = moonPhaseTerms(i);
    if (term == 0 && moonPhaseTime + 4 < Period.endTime) {
      (0 until 4).foreach { d =>
        val day = (moonPhaseTime + d - Period.startTime).toInt;
        val (azi, alt) = calcPlanetOnSunsetTime(day, JplData.Moon);
        if (azi >= aziThres0 && azi <= aziThres1 && alt >= altThres0) {
          val day = (moonPhaseTime + d - Period.startTime).toInt;
          val (azi, alt) = calcPlanetOnSunsetTime(day, JplData.Moon);
          tweetsManager.putTweet(SunsetMoonTweetContent(day, azi, alt));
        }
      }
    }
  }
}

// 水曜・金曜
{
  val aziThres0 = 200 / Const.PI57;
  val aziThres1 = 315 / Const.PI57;
  val altThres0 = 10 / Const.PI57;
  val diffThreashold = 0.1 / Const.PI57;

  val planets = IndexedSeq(("金星", JplData.Venus, 5), ("水星", JplData.Mercury, 3));
  innerPlanets.zipWithIndex.foreach { case (planet, pi) =>
    (1 until Period.period).foreach { day =>
      val wday = TimeLib.wday(Period.startTime + day);
      val sunsetTweets = tweetsManager.getTweets(Period.startTime + day).sunsetTweets;
      if (wday == planet._3 && !sunsetTweets.exists(_.starNames.contains(planet._1))) {
        val (azi, alt) = innerPlanetsSunsetAziAltList(pi)(day);
        if (alt >= altThres0) {
          val prevAlt = innerPlanetsSunsetAziAltList(pi)(day - 1)._2;
          val isIncreasing = (alt >= prevAlt + diffThreashold);
          val isDecreasing = (alt <= prevAlt - diffThreashold);
          tweetsManager.putTweet(SunsetPlanetTweetContent(day, planet._1, azi, alt, isIncreasing, isDecreasing, false));
        }
      }
    }
  }
}

// 日没ツイートがある場合に他の天体のツイートも追加
{
  val aziThres0 = 200 / Const.PI57;
  val aziThres1 = 315 / Const.PI57;
  val altThres0 = 10 / Const.PI57;
  val diffThreashold = 0.1 / Const.PI57;

  val planets = IndexedSeq(("金星", JplData.Venus), ("水星", JplData.Mercury));
  val planets2 = IndexedSeq(("火星", JplData.Mars), ("木星", JplData.Jupiter), ("土星", JplData.Saturn));
  (1 until Period.period).foreach { day =>
    val sunsetTweets = tweetsManager.getTweets(Period.startTime + day).sunsetTweets;
    if (sunsetTweets.nonEmpty) {
      planets.foreach { p =>
        if (!sunsetTweets.exists(_.starNames.contains(p._1))) {
          val (azi, alt) = calcPlanetOnSunsetTime(day, p._2);
          if (azi >= aziThres0 && azi <= aziThres1 && alt >= altThres0) {
            val (_, prevAlt) = calcPlanetOnSunsetTime(day - 1, p._2);
            val isIncreasing = (alt >= prevAlt + diffThreashold);
            val isDecreasing = (alt <= prevAlt - diffThreashold);
            tweetsManager.putTweet(SunsetPlanetTweetContent(day, p._1, azi, alt, isIncreasing, isDecreasing, false));
          }
        }
      };
      {
        val p = ("月", JplData.Moon);
        if (!sunsetTweets.exists(_.starNames.contains(p._1))) {
          val (azi, alt) = calcPlanetOnSunsetTime(day, p._2);
          if (azi >= aziThres0 && azi <= aziThres1 && alt >= altThres0) {
            tweetsManager.putTweet(SunsetMoonTweetContent(day, azi, alt));
          }
        }
      }
      planets2.foreach { p =>
          val (azi, alt) = calcPlanetOnSunsetTime(day, p._2);
          if (azi >= aziThres0 && azi <= aziThres1 && alt >= altThres0) {
            tweetsManager.putTweet(SunsetStarTweetContent(day, p._1, azi, alt));
          }
      };
    }
  }
}

def calcNextSeason(aziAltList: IndexedSeq[(Double, Double)], day: Int): (Int, Int) = {
  val altThres = 15 / Const.PI57;
  val p1 = aziAltList.indexWhere(_._2 >= altThres, day);
  if (p1 >= 0) {
    val p2 = aziAltList.indexWhere(_._2 < altThres, p1) - 1;
    if (p2 >= 0) {
      (p1, p2);
    } else {
      (p1, -1);
    }
  } else {
    (-1, -1);
  }
}

// 水曜・金曜だけど水星・金星が見えない場合
{
  val altThres = 15 / Const.PI57;

  innerPlanets.zipWithIndex.foreach { case (planet, pi) =>
    (1 until Period.period).foreach { day =>
      if ((Period.startTime + day).toInt % 2 == 0) {
        val wday = TimeLib.wday(Period.startTime + day);
        val sunsetTweets = tweetsManager.getTweets(Period.startTime + day).sunsetTweets;
        if (wday == planet._3 && !sunsetTweets.exists(_.starNames.contains(planet._1))) {
          val (azi, alt) = innerPlanetsSunsetAziAltList(pi)(day);
          if (alt < altThres) {
            val time = Period.startTime + day + 12.5 / 24;
            val (p1, p2) = calcNextSeason(innerPlanetsSunsetAziAltList(pi), day);
            tweetsManager.putTweet(SunsetNextPlanetTweetContent(time, planet._1, p1, p2));
          }
        }
      }
    }
  }
}

// 日没時間
{
  MathLib.getMaxMinUpDownFlagListDiscrete(0, Period.period, 90) { day =>
    sunsetTimes(day) - Period.startTime - day;
  }.zipWithIndex.foreach { case (flag, day) =>
    val wday = TimeLib.wday(Period.startTime + day);
    if (flag == 1 || flag == 3 || wday == 0) {
      tweetsManager.putTweet(SunsetTweetContent(day, flag));
    }
  }
}

//==============================================================================
// 月・惑星と恒星の会合
//==============================================================================

case class CloseStarsTweetContent(rawTime: Double, stepCountPerDay: Int, slowStarName: String, fastStarName: String,
  distance: Double, hashtags: List[String]) extends TweetContent {
  def time: Double = TimeLib.round(rawTime, stepCountPerDay);
  def message: String = "%sが%s%s%s".format(fastStarName, slowStarName, distanceStr, moonPhaseStr);
  def distanceStr: String = {
    val distance360 = distance * Const.PI57;
    if (distance360 < 1.0) {
      "に接近 (1°未満)";
    } else if (distance360 < 2.0) {
      "に接近 (2°未満)";
    } else if (distance360 < 3.0) {
      "に接近 (3°未満)";
    } else {
      "の近くにいます";
    }
  }
  def moonPhaseStr: String = {
    if (fastStarName == "月") {
      "。%s".format(calcMoonPhaseString(time));
    } else {
      "";
    }
  }
  def starNames: List[String] = List(slowStarName, fastStarName);
}

{
  val altThres0 = 10 / Const.PI57;
  val distanceThres = 3.0 / Const.PI57;
  val distanceThresMoon = 6.0 / Const.PI57;

  def calcClosestMoon(slowStarName: String, fastStarName: String, hashtags: List[String],
    slowStarXyzFunc: Double => Array[Double],
    fastStarXyzFunc: Double => Array[Double]): Unit = {
    MathLib.findMaxMinListContinuous(Period.startTime, Period.endTime, 10.0, 24 * 6) { time =>
      val utc = time;
      val tdb = TimeLib.mjdutcToTdb(utc);
      val xyz_s = slowStarXyzFunc(tdb);
      val xyz_f = fastStarXyzFunc(tdb);
      VectorLib.angularDistance1(xyz_s, xyz_f);
    }.foreach { case (time, flag) =>
      if (flag > 0 && isNightTime2(time)) {
        val tdb = TimeLib.mjdutcToTdb(time);
        val xyz_f = fastStarXyzFunc(tdb);
        val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
        val xyz_f2 = VectorLib.multiplyMV(bpnMatrix, xyz_f);
        val (azi, alt) = hcs.trueEquatorialXyzToAziAltFromUtc(xyz_f2, time);
        if (alt >= altThres0) {
          val xyz_s = slowStarXyzFunc(tdb);
          val distance = VectorLib.angularDistance(xyz_s, xyz_f);
          if (distance < distanceThresMoon) {
            tweetsManager.putTweet(CloseStarsTweetContent(time, 24 * 6, slowStarName, fastStarName, distance, hashtags));
          }
        }
      }
    }
  }
  def calcClosest2(slowStarName: String, fastStarName: String, hashtags: List[String],
    slowStarXyzFunc: Double => Array[Double],
    fastStarXyzFunc: Double => Array[Double]): Unit = {
    MathLib.findMaxMinListContinuous(Period.startTime, Period.endTime, 10.0, 24) { time =>
      val utc = time;
      val tdb = TimeLib.mjdutcToTdb(utc);
      val xyz_s = slowStarXyzFunc(tdb);
      val xyz_f = fastStarXyzFunc(tdb);
      VectorLib.angularDistance1(xyz_s, xyz_f);
    }.foreach { case (time, flag) =>
      if (flag > 0) {
        val ut1 = time; // 近似的
        val tdb = TimeLib.mjdutcToTdb(time);
        val xyz_f = fastStarXyzFunc(tdb);
        val xyz_s = slowStarXyzFunc(tdb);
        val distance = VectorLib.angularDistance(xyz_s, xyz_f);
        if (distance < distanceThres) {
          getConjunctionTweetTime(time, xyz_f) match {
            case Some(postTime) =>
              tweetsManager.putTweet(CloseStarsTweetContent(postTime, 24, slowStarName, fastStarName, distance, hashtags));
            case None =>
              // nop
          }
        }
      }
    }
  }

  val stars0: IndexedSeq[(String, Array[Double], List[String])] = IndexedSeq(
    ("03h47m", "+24°06′", "おうし座すばる", List("プレアデス星団")),
    ("04h36m", "+16°31′", "おうし座アルデバラン", Nil),
    ("07h45m", "+28°02′", "ふたご座ポルックス", Nil),
    ("10h08m", "+11°58′", "しし座レグルス", Nil),
    ("13h25m", "-11°09′", "おとめ座スピカ", Nil),
    ("16h29m", "-26°26′", "さそり座アンタレス", Nil),
  ).map { case (lngStr, latStr, name, hashtags) =>
    val lng = (lngStr.substring(0, 2).toInt.toDouble + lngStr.substring(3, 5).toInt.toDouble / 60) / 24 * Const.PI2;
    val lat = (latStr.substring(0, 3).toInt.toDouble + latStr.substring(4, 6).toInt.toDouble / 60) / 360 * Const.PI2;
    val c = Math.cos(lat);
    val x = c * Math.cos(lng);
    val y = c * Math.sin(lng);
    val z = Math.sin(lat);
    (name, Array(x, y, z), hashtags);
  }
  val stars1 = IndexedSeq(
    ("水星", JplData.Mercury),
    ("金星", JplData.Venus),
    ("火星", JplData.Mars),
    ("木星", JplData.Jupiter),
    ("土星", JplData.Saturn),
  );
  val stars2 = IndexedSeq(
    ("月", JplData.Moon),
  );
  stars2.foreach { star2 =>
    stars0.foreach { star0 =>
      calcClosestMoon(star0._1, star2._1, star0._3,
        { tdb: Double =>
          star0._2;
        },
        { tdb: Double =>
          JplData.calcPlanetFromEarth(tdb, star2._2);
        });
    }
  }
  stars2.foreach { star2 =>
    stars1.foreach { star1 =>
      calcClosestMoon(star1._1, star2._1, List(star1._1),
        { tdb: Double =>
          JplData.calcPlanetFromEarth(tdb, star1._2);
        },
        { tdb: Double =>
          JplData.calcPlanetFromEarth(tdb, star2._2);
        });
    }
  }
  stars1.foreach { star1 =>
    stars0.foreach { star0 =>
      calcClosest2(star0._1, star1._1, star1._1 :: star0._3,
        { tdb: Double =>
          star0._2;
        },
        { tdb: Double =>
          JplData.calcPlanetFromEarth(tdb, star1._2);
        });
    }
  }
  (0 until stars1.size).foreach { i =>
    ((i + 1) until stars1.size).foreach { j =>
      val star1 = stars1(i);
      val star0 = stars1(j);
      calcClosest2(star0._1, star1._1, List(star1._1, star0._1),
        { tdb: Double =>
          JplData.calcPlanetFromEarth(tdb, star0._2);
        },
        { tdb: Double =>
          JplData.calcPlanetFromEarth(tdb, star1._2);
        });
    }
  }
}

//==============================================================================
// 21時・23時の月
//==============================================================================

{
  val altThres = 10 / Const.PI57;
  (0 until Period.period).foreach { d =>
    val time = Period.startTime + d + 21.0 / 24.0 - 1.0 / (24 * 6);
    if (!tweetsManager.getTweets(time).tweets.filter(tc => isNightTime0(tc.time)).flatMap(_.starNames).contains("月")) {
      {
        val (azi, alt) = Acs.Horizontal.calcPlanetAziAlt(time, JplData.Moon, hcs);
        if (alt >= altThres) {
          Some((time, azi, alt));
        } else {
          val time = Period.startTime + d + 23.0 / 24.0;
          val (azi, alt) = Acs.Horizontal.calcPlanetAziAlt(time, JplData.Moon, hcs);
          if (alt >= altThres) {
            Some((time, azi, alt));
          } else {
            None;
          }
        }
      } match {
        case None => ;
        case Some((time, azi, alt)) =>
          val xyz = Acs.Icrs.calcPlanetXyz(time, JplData.Moon);
          val (cons, hashtags) = Constellations.icrsToConstellation(xyz);
          val hcsStr = Hcs.aziAltToNaturalString(azi, alt);
          if (azi >= Const.PI) {
            tweetsManager.putTweet(time, "月は%s、%sにいます。%s".format(hcsStr, cons, calcMoonPhaseString(time)) +
              hashtags.map(" #" + _).mkString);
          } else {
            tweetsManager.putTweet(time, "月は%s、%sにいます".format(hcsStr, cons) +
              hashtags.map(" #" + _).mkString);
          }
      }
    }
  }
}

//==============================================================================
// 21時・23時の惑星
//==============================================================================

OuterPlanets.nightTweets(hcs).foreach { tw =>
  tweetsManager.putTweet(tw);
}

// 惑星が見えない場合
{
  val altThres = 15 / Const.PI57;

  outerPlanets.zipWithIndex.foreach { case (planet, pi) =>
    (0 until Period.period).foreach { day =>
      if ((Period.startTime + day).toInt % 2 == 0) {
        val wday = TimeLib.wday(Period.startTime + day);
        val tweets = tweetsManager.getTweets(Period.startTime + day).tweets;
        if (wday == planet._3 && !tweets.exists(_.starNames.contains(planet._1))) {
          if (!(-6 to +6).exists { i =>
            tweetsManager.getTweets(Period.startTime + day + i).tweets.exists {
              case PlanetAstronomyTweetContent(_, _, planetName) if (planetName == planet._1) => true;
              case _ => false;
            }
          }) {
            val (azi, alt) = outerPlanetsNightAziAltList(pi)(day);
            if (alt < altThres) {
              val time = Period.startTime + day + 12.5 / 24;
              val (p1, p2) = calcNextSeason(outerPlanetsNightAziAltList(pi), day);
              tweetsManager.putTweet(NightNextPlanetTweetContent(time, planet._1, p1, p2));
            }
          }
        }
      }
    }
  }
}

//==============================================================================
// 月の出・月の入り
//==============================================================================

tweetMoonRiseSet();

//==============================================================================
// 手動のツイート
//==============================================================================

Diff.removeTweets1(tweetsManager);
Diff.putTweets1(tweetsManager);

//==============================================================================
// 星座
//==============================================================================

{
  import Ordering.Double.IeeeOrdering;
  case class StarTweetContent(time: Double, message: String,
    override val urlOpt: Option[String], hashtags2: List[String], starNames: List[String]) extends TweetContent {
    def hashtags: List[String] = hashtags2 ::: "星空" :: "星座" :: Nil;
  }

  // この時期21時ごろ見えやすい星座
  def putTweetConstellations21(day: Int): Unit = {
    val altThres = 30.0 / Const.PI57;
    val time = Period.startTime + day + 21.0 / 24;
    val tdb = TimeLib.mjdutcToTdb(time);
    val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
    val constellations = Constellations.constellationData.flatMap { constellation =>
      val xyz2 = VectorLib.multiplyMV(bpnMatrix, constellation.xyz);
      val (azi, alt) = hcs.trueEquatorialXyzToAziAltFromUtc(xyz2, time);
      if (alt >= altThres) {
        Some((azi, constellation.name));
      } else {
        None;
      }
    }
    val constellationsStr = constellations.sortBy(-_._1).map(_._2).mkString("、");
    val starNames = constellations.sortBy(-_._1).map(_._2).toList;
    val msg = "この時期21時ごろ見えやすい星座は、%sです #星空 #星座".format(constellationsStr);
    tweetsManager.putTweet(Period.startTime + day + (12.0 + 35.0 / 60) / 24, msg, starNames);
  }

  // この時期21時ごろ見える南の空低い星座
  def putTweetConstellationsSouth(day: Int): Boolean = {
    val decThres = - Const.PI5 + tokyoLat + 30.0 / Const.PI57;
    val time = Period.startTime + day + 21.0 / 24;
    val tdb = TimeLib.mjdutcToTdb(time);
    val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
    val constellations = Constellations.constellationData.flatMap { constellation =>
      if (constellation.dec < decThres) {
        val xyz2 = VectorLib.multiplyMV(bpnMatrix, constellation.xyz);
        val (azi, alt) = hcs.trueEquatorialXyzToAziAltFromUtc(xyz2, time);
        if (alt > 0.0) {
          Some((azi, constellation.name));
        } else {
          None;
        }
      } else {
        None;
      }
    }
    if (constellations.nonEmpty) {
      val constellationsStr = constellations.sortBy(-_._1).map(_._2).mkString("、");
      val starNames = constellations.sortBy(-_._1).map(_._2).toList;
      val msg = "この時期21時ごろ南の低い空に見える星座は、%sです #星空 #星座".format(constellationsStr);
      tweetsManager.putTweet(Period.startTime + day + (12.0 + 35.0 / 60) / 24, msg, starNames);
      true;
    } else {
      false;
    }
  }

  // この時期21時ごろ見える明るい星
  def putTweetBrightStars21(day: Int): Unit = {
    val altThres = 10.0 / Const.PI57;
    //val altThres = 0.0;
    val time = Period.startTime + day + 21.0 / 24;
    val tdb = TimeLib.mjdutcToTdb(time);
    val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
    val constellations = (Constellations.starData.flatMap { case (ra, xyz, name, hashtags) =>
      val xyz2 = VectorLib.multiplyMV(bpnMatrix, xyz);
      val (azi, alt) = hcs.trueEquatorialXyzToAziAltFromUtc(xyz2, time);
      if (alt >= altThres) {
        Some((azi, name));
      } else {
        None;
      }
    }) ++
    (List(("月", JplData.Moon), ("金星", JplData.Venus), ("火星", JplData.Mars), ("木星", JplData.Jupiter), ("土星", JplData.Saturn)).
      flatMap { case (name, targetPlanet) =>
      val (azi, alt) = Acs.Horizontal.calcPlanetAziAlt(time, targetPlanet, hcs);
      if (alt >= altThres) {
        Some((azi, name));
      } else {
        None;
      }
    }).toIndexedSeq;
    val constellationsStr = constellations.sortBy(-_._1).map(_._2).mkString("、");
    val starNames = constellations.sortBy(-_._1).map(_._2).toList;
    val msg = "この時期21時ごろ見える明るい星は、%sです #星空".format(constellationsStr);
    tweetsManager.putTweet(Period.startTime + day + (12.0 + 35.0 / 60) / 24, msg, starNames);
  }

  // この時期日没時の一番星となりうる明るい星
  def putTweetFirstStar(day: Int): Unit = {
    val altHor = -0.90 / Const.PI57;
    val altThres30 = 30.0 / Const.PI57;
    //val altThres = 0.0;
    val time = sunsetTimes(day);
    val tdb = TimeLib.mjdutcToTdb(time);
    val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
    val constellations0 = ((Constellations.starData2.flatMap { case (ra, xyz, magnitude, name, hashtags) =>
      val xyz2 = VectorLib.multiplyMV(bpnMatrix, xyz);
      val (azi, alt) = hcs.trueEquatorialXyzToAziAltFromUtc(xyz2, time);
      val aziAltStr = Hcs.aziAltToNaturalString(azi, alt);
      if (alt >= altHor) {
        Some((azi, alt, name, aziAltStr, magnitude));
      } else {
        None;
      }
    }) ++
    (List(("水星", JplData.Mercury, -2.5), ("金星", JplData.Venus, -4.9), ("火星", JplData.Mars, -2.8), ("木星", JplData.Jupiter, -2.9), ("土星", JplData.Saturn, -0.6)).
      flatMap { case (name, targetPlanet, magnitude) =>
      val (azi, alt) = Acs.Horizontal.calcPlanetAziAlt(time, targetPlanet, hcs);
      val aziAltStr = Hcs.aziAltToNaturalString(azi, alt);
      if (alt >= altHor) {
        Some((azi, alt, name, aziAltStr, magnitude));
      } else {
        None;
      }
    })).toIndexedSeq.map { case (azi, alt, name, aziAltStr, magnitude) =>
      val magnitude2 = if (alt >= altThres30) {
        magnitude;
      } else {
        magnitude + (altThres30 - alt) * Const.PI57 / 8.0;
      }
      (azi, alt, name, aziAltStr, magnitude, magnitude2);
    }.sortBy(_._6);
    val constellations = (
      constellations0.take(3).toSet ++
      constellations0.filter(_._2 >= altThres30).take(3).toSet
    ).toIndexedSeq.sortBy(_._6);
    if (constellations.nonEmpty) {
      val constellationsStr = constellations.map { case (azi, alt, name, aziAltStr, magnitude, magnitude2) =>
        "%s(%s)".format(name, aziAltStr);
      }.mkString("、");
      val starNames = constellations.map(_._3).toList;
      val msg = "この時期日没時の一番星となりうる明るい星は、%sです".format(constellationsStr);
      tweetsManager.putTweet(SunsetFirstStarTweetContent(day, msg, starNames));
    }
  }

  // この時期21時ごろ見えやすい黄道十二星座
  def putTweetConstellationsEcliptical(day: Int): Boolean = {
    val altThres = 30.0 / Const.PI57;
    val time = Period.startTime + day + 21.0 / 24;
    val tdb = TimeLib.mjdutcToTdb(time);
    val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
    val constellations = Constellations.constellationData.filter(_.eclipticalFlag).flatMap { constellation =>
      val xyz2 = VectorLib.multiplyMV(bpnMatrix, constellation.xyz);
      val (azi, alt) = hcs.trueEquatorialXyzToAziAltFromUtc(xyz2, time);
      if (alt >= altThres) {
        val aziAltStr = Hcs.aziAltToNaturalString(azi, alt);
        Some((aziAltStr, azi, constellation.name));
      } else {
        None;
      }
    }
    if (constellations.nonEmpty) {
      val constellationsStr = constellations.sortBy(-_._2).map { case (aziAltStr, azi, name) =>
        "%s(%s)".format(name, aziAltStr);
      }.mkString("、");
      val starNames = constellations.sortBy(-_._2).map(_._3).toList;
      val msg = "この時期21時ごろ見えやすい黄道十二星座は、%sです #星空 #星座".format(constellationsStr);
      tweetsManager.putTweet(Period.startTime + day + (12.0 + 35.0 / 60) / 24, msg, starNames);
      true;
    } else {
      false;
    }
  }

  // この時期21時ごろ見える天の川
  def putTweetConstellationsGalaxy(day: Int): Boolean = {
    val altThres = 30.0 / Const.PI57;
    val time = Period.startTime + day + 21.0 / 24;
    val tdb = TimeLib.mjdutcToTdb(time);
    val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
    val constellations = Constellations.constellationData.filter(_.galaxyFlag).flatMap { constellation =>
      val xyz2 = VectorLib.multiplyMV(bpnMatrix, constellation.xyz);
      val (azi, alt) = hcs.trueEquatorialXyzToAziAltFromUtc(xyz2, time);
      if (alt >= altThres) {
        val sid = hcs.siderealTimeFromUtc(time);
        val ra2 = constellation.ra - sid;
        val ra3 = if (ra2 < -Const.PI) ra2 + Const.PI2 else if (ra2 >= Const.PI) ra2 - Const.PI2 else ra2;
        Some((ra3, constellation.name));
      } else {
        None;
      }
    }
    if (constellations.nonEmpty) {
      val constellationsStr = constellations.sortBy(_._1).map(_._2).mkString("、");
      val starNames = constellations.sortBy(-_._1).map(_._2).toList;
      val msg = "この時期21時ごろ見える天の川は、%sを通っています #星空 #星座".format(constellationsStr);
      tweetsManager.putTweet(Period.startTime + day + (12.0 + 35.0 / 60) / 24, msg, starNames);
      true;
    } else {
      false;
    }
  }

  // 恒星の南中
  def putTweetCulminations(): Unit = {
    val altHor = -0.90 / Const.PI57;
    var index: Int = -1;
    val culminationContents = Constellations.culminationContents;
    (115 until Period.period).foreach { day => // PERIOD 2021/07/24
      if (index < 0) {
        val time = Period.startTime + day + 21.0 / 24.0; // PERIOD
        val sid = hcs.siderealTimeFromUtc(time);
        index = culminationContents.indexWhere(_._1 > sid);
        if (index < 0) {
          index = 0;
        }
      } else {
        val timeS = Period.startTime + day + 21.0 / 24.0;
        val timeE = Period.startTime + day + 23.0 / 24.0;
        val time1 = hcs.siderealTimeToUtc(culminationContents(index)._1, timeS);
        val time0 = if (time1 < timeS) {
          time1;
        } else if (tweetsManager.getTweets(timeS).tweets.map(_.time).filter(isNightTime3).isEmpty) {
          if (time1 < timeE) {
            time1;
          } else {
            0.0;
          }
        } else {
          0.0;
        }
        if (time0 > 0.0) {
          val msg = if (Acs.Horizontal.calcPlanetAlt(time0, JplData.Moon, hcs) >= altHor) {
            culminationContents(index)._2;
          } else {
            culminationContents(index)._2 + "。月明かりなし";
          }
          val urlOpt = culminationContents(index)._3;
          val hashtags = culminationContents(index)._4;
          val starNames = culminationContents(index)._5;
          tweetsManager.putTweet(StarTweetContent(time0, msg, urlOpt, hashtags, starNames));
          index += 1;
          if (index == culminationContents.size) {
            index = 0;
          }
        }
      }
    }
  }

  // 星座の解説など(かんむり座のパターン)
  def putTweetLunchTimeContents(): Unit = {
    val lunchTimeContents = Constellations.lunchTimeContents;
    var day1: Int = 123; // PERIOD 2021/08/01
    var day2: Int = day1;
    var index: Int = {
      val time = Period.startTime + day1 + 21.0 / 24.0;
      val sid = hcs.siderealTimeFromUtc(time);
      val index = lunchTimeContents.indexWhere(_._1 > sid);
      if (index < 0) {
        0;
      } else {
        index;
      }
    }
    while (day1 < Period.period) {
      val sid = hcs.siderealTimeFromUtc(Period.startTime + day2 + 21.0 / 24);
      if (MathLib.circleDiff1(sid, lunchTimeContents(index)._1) >= 0) {
        val msg = lunchTimeContents(index)._2;
        {
          val p = (day2 until (day1 + 14)).indexWhere { d =>
            !tweetsManager.getTweets(Period.startTime + d).tweets.exists(tc => DateTweets.isDayTime(tc.time));
          }
          day2 = if (p < 0) day2 else day2 + p;
        }
        val urlOpt = lunchTimeContents(index)._3;
        val hashtags = lunchTimeContents(index)._4;
        val starNames = Nil;
        tweetsManager.putTweet(StarTweetContent(Period.startTime + day2 + 12.0 / 24 + 5.0 / 60 / 24, msg, urlOpt, hashtags, starNames));
        index += 1;
        if (index == lunchTimeContents.size) {
          index = 0;
        }
        day2 += 1;
      }
      day1 += 1;
      if (day2 < day1) {
        day2 = day1;
      }
    }
  }

  putTweetCulminations();

  {
    var kind: Int = 0;
    (31 until Period.period).foreach { day => // PERIOD
      val wday = TimeLib.wday(Period.startTime + day);
      if (wday == 1) { // 月曜
        if ((Period.startTime + day).toInt % 2 == 0) {
          putTweetConstellations21(day);
        } else {
          var i: Int = 0;
          var f = false;
          while (!f && i < 3) {
            if (kind == 0) {
              f = putTweetConstellationsGalaxy(day);
              kind = 1;
            } else if (kind == 1) {
              f = putTweetConstellationsEcliptical(day);
              kind = 2;
            } else if (kind == 2) {
              f = putTweetConstellationsSouth(day);
              kind = 0;
            }
            i += 1;
          }
        }
      } else if (wday == 3) { // 水曜
        if ((Period.startTime + day).toInt % 2 == 1) {
          putTweetBrightStars21(day);
        }
      }
    }
  }
  (40 until Period.period).foreach { day => // PERIOD
    val wday = TimeLib.wday(Period.startTime + day);
    if (wday == 1) {
      if ((Period.startTime + day).toInt % 2 == 0) {
        if (tweetsManager.getTweets(Period.startTime + day).sunsetTweets.isEmpty) {
          //putTweetFirstStar(day);
        }
      }
    }
  }

  putTweetLunchTimeContents();

  (words.periodStart until Period.period).foreach { day =>
    if (tweetsManager.getTweets(Period.startTime + day).size < 3 || !tweetsManager.getTweets(Period.startTime + day).tweets.map(_.time).exists(isLunchTime)) {
      words.putTweetWord(day, hcs, tweetsManager);
    }
  }
}

//==============================================================================
// 手動のツイート
//==============================================================================

Diff.removeTweets2(tweetsManager);
Diff.putTweets2(tweetsManager);

//==============================================================================

  // TODO 上に書かれているコードを整理して少しずつmainの中に移動させるつもり

  def main(args: Array[String]): Unit = {

    words.saveHistory(wordHistoryPath);
    writeEmptyTweets();
    writeTweets();
  }

  // なにもツイートのない日付
  def writeEmptyTweets(): Unit = {
    (0 until Period.period).foreach { day =>
      val time = Period.startTime + day;
      if (tweetsManager.getTweets(time).isEmpty) {
        tweetsManager.putTweet(time, "##empty");
      } else {
        //if (tweetsManager.getTweets(time).tweets.map(_.time).filter(isNightTime0).isEmpty) {
        if (tweetsManager.getTweets(time).tweets.map(_.time).filter(isNightTime3).isEmpty) {
          tweetsManager.putTweet(time + 23.0 / 24, "##night empty");
        } else if (!tweetsManager.getTweets(time).tweets.exists(tc => DateTweets.isDayTime(tc.time))) {
          tweetsManager.putTweet(time + 9.0 / 24, "##daytime empty");
        }
      }
    }
  }

  // ツイート出力
  def writeTweets(): Unit = {
    scala.util.Using(new java.io.PrintWriter(new java.io.FileOutputStream(dataPath))) { writer =>
      (0 until Period.period).foreach { day =>
        tweetsManager.getTweets(Period.startTime + day).tweets.foreach { case tc =>
          val time = tc.time;
          if (time >= Period.startTime1 && time < Period.endTime) {
            val msg = tc.tweetContent;
            writer.println("%s %s".format(TimeLib.timeToDateTimeString(time), msg));
          }
        }
      }
    }
  }

}

