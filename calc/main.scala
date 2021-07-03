
object Main {

// PERIOD
val startTime = TimeLib.stringToModifiedJulianDay("2021-01-31T00:00:00+09:00");
val startTime1 = TimeLib.stringToModifiedJulianDay("2021-07-03T00:00:00+09:00");
val endTime1 = TimeLib.stringToModifiedJulianDay("2022-06-06T00:00:00+09:00");
val endTime = TimeLib.stringToModifiedJulianDay("2022-06-07T00:00:00+09:00");

val period = (endTime - startTime).toInt;


val jplDataPath = "../var/ssd.jpl.nasa.gov/pub/eph/planets/ascii/de430/ascp1950.430";
val nutLsDataPath = "../static/nut-ls.txt";
val nutPlDataPath = "../static/nut-pl.txt";
val constellationsDataPath = "constellations.txt";
val holidayDataPath = "holiday.txt";
val meteorDataPath = "meteor.txt";
val diffDataPath = "diff.txt";

val wordHistoryPath = "../etc/word-history.txt";

val dataPath = "../data.txt";

val PI = Math.PI;
val PI2 = Math.PI * 2.0;
val PI2R = 1.0 / PI2;
val PI5 = Math.PI * 0.5;
val PI57 = 180.0 / Math.PI;

val AU = 1.49597870700000000e+8;

val tokyoLng = 139.7 / PI57;
val tokyoLat = 35.7 / PI57;

System.err.println("Started calculating...");

object StringLib {
  def parseContent(content: String): (String, Option[String], List[String], List[String]) = {
    val UrlPattern = "(.+)\\s(https?://\\S+)".r;
    val HashtagsPattern = "(.+)\\s+#(\\S+)".r;
    val StarNamePattern = "(.+)\\s+\\$(\\S+)".r;
    var content0 = content.trim;
    var urlOpt: Option[String] = None;
    var starNames: List[String] = Nil;
    var hashtags: List[String] = Nil;
    var f: Boolean = true;
    while (f) {
      content0 match {
        case UrlPattern(c, s) =>
          content0 = c;
          urlOpt = Some(s);
        case HashtagsPattern(c, s) =>
          content0 = c;
          hashtags = s :: hashtags;
        case StarNamePattern(c, s) =>
          content0 = c;
          starNames = s :: starNames;
        case _ =>
          f = false;
      }
    }
    (content0, urlOpt, hashtags.reverse, starNames.reverse);
  }
}

val bpn = new Bpn(nutLsDataPath, nutPlDataPath);

class Hcs(lng: Double, lat: Double) {

  def siderealTimeFromUtc(utc: Double): Double = {
    val ut1 = utc; // 近似的
    siderealTimeFromUt1(ut1);
  }

  def siderealTimeFromUt1(ut1: Double): Double = {
    val gmst = TimeLib.mjdut1ToGmst(ut1);
    val gast = gmst; // 近似
    val s = gast + lng;
    if (s >= PI2) {
      s - PI2;
    } else {
      s;
    }
  }

  def siderealTimeToUtc(sid: Double, time0: Double): Double = {
    siderealTimeToUt1(sid, time0); // 近似
  }

  def siderealTimeToUt1(sid: Double, time0: Double): Double = {
    val s = sid - lng;
    val s2 = if (s < 0) {
      s + PI2;
    } else {
      s;
    }
    TimeLib.gmstToMjdut1(s2, time0);
  }

  def trueEquatorialXyzToAziAltFromUtc(xyz: Array[Double], utc: Double): (Double, Double) = {
    val ut1 = utc; // 近似的
    trueEquatorialXyzToAziAltFromUt1(xyz, ut1);
  }

  def trueEquatorialXyzToAziAltFromUt1(xyz: Array[Double], ut1: Double): (Double, Double) = {
    var r: Array[Double] = VectorLib.unitMatrix;
    r = VectorLib.rotateMatrixZ(PI5 - siderealTimeFromUt1(ut1), r);
    r = VectorLib.rotateMatrixX(PI5 - lat, r);
    // 地平座標系 X:西 Y:南 Z:天頂
    val hcs = VectorLib.multiplyMV(r, xyz);
    val x = hcs(0);
    val y = hcs(1);
    val z = hcs(2);
    val azi1 = Math.atan2(-x, -y);
    val azi = if (azi1 < 0) azi1 + PI2 else azi1; // 北0°、東90°、南180°、西270°
    val xy = Math.sqrt(x * x + y * y);
    val alt = Math.atan2(z, xy);
    (azi, alt);
  }

}

object Hcs {

  def aziAltToNaturalString(azi: Double, alt: Double): String = {
    val altHor = -0.90 / PI57;
    val alt360 = alt * PI57;
    if (alt360 >= 80) {
      "天頂付近";
    } else {
      val azi360 = azi * PI57;
      val s1 = if (azi360 <= 22.5 || azi360 >= 337.5) {
        "北の";
      } else if (azi360 < 67.5) {
        "北東の";
      } else if (azi360 > 292.5) {
        "北西の";
      } else if (azi360 <= 112.5) {
        "東の";
      } else if (azi360 >= 247.5) {
        "西の";
      } else if (azi360 < 157.5) {
        "南東の";
      } else if (azi360 > 202.5) {
        "南西の";
      } else {
        "南の";
      }
      val s2 = if (alt < altHor) {
        "地平線の下";
      } else if (alt360 < 10.0) {
        "空地平線近く";
      } else if (alt360 < 30.0) {
        "空低く";
      } else if (alt360 >= 60.0) {
        "空高く";
      } else {
        "空";
      }
      s1 + s2;
    }
  }

  def aziAltToNaturalStringC(azi: Double, alt: Double): (String, Double) = {
    val alt360 = alt * PI57;
    if (alt360 >= 80) {
      ("天頂", 0.0);
    } else if (alt360 < 30) {
      ("", 0.0);
    } else {
      val azi360 = azi * PI57;
      if (azi360 <= 45) {
        ("北の空", -alt);
      } else if (azi360 >= 315) {
        ("北の空", -alt);
      } else if (azi360 <= 135) {
        ("東の空", -alt);
      } else if (azi360 >= 225) {
        ("西の空", -alt);
      } else {
        ("南の空", -alt);
      }
    }
  }

}

object Constellations {

  case class Constellation (ra: Double, dec: Double, xyz: Array[Double], name: String, hashtags: List[String],
    eclipticalFlag: Boolean, galaxyFlag: Boolean);

  def siderealTimeToConstellation(siderealTime: Double, data: IndexedSeq[(String, String)]): String = {
    val lng5 = (siderealTime * PI57 / 5).toInt * 5;
    val key = "%02dh%02dm".format(lng5 / 15, lng5 % 15 * 4);
    val p = data.indexWhere(_._1 > key) - 1;
    val cons = if (p == -1) {
      data(data.size - 2)._2;
    } else if (p >= 0) {
      data(p)._2;
    } else {
      throw new Exception();
    }
    cons;
  }

}

class Constellations {

  val (culminationContents, lunchTimeContents, lunchTimeContents2, constellationData, starData, starData2, constellationsMap) = {
    import Ordering.Double.IeeeOrdering;
    val source = scala.io.Source.fromFile(constellationsDataPath);
    var culminationContents: List[(Double, String, Option[String], List[String], List[String])] = Nil;
    var lunchTimeContents:   List[(Double, String, Option[String], List[String])] = Nil;
    var lunchTimeContents2:   List[Words.ConstellationLunchTimeContent] = Nil;
    var constellationData: List[Constellations.Constellation] = Nil;
    var starData: List[(Double, Array[Double], String, List[String])] = Nil;
    var starData2: List[(Double, Array[Double], Double, String, List[String])] = Nil;
    var constellationsMap: Map[String, (String, List[String])] = Map.empty;
    source.getLines.foreach { line =>
      if (!line.startsWith("#") && line.length > 7) {
        val raStr = line.substring(0, 6);
        val (content0: String, urlOpt: Option[String], hashtags: List[String], starNames: List[String]) = StringLib.parseContent(line.substring(7));
        def calcRa(raStr: String): Double = {
          (raStr.substring(0, 2).toInt.toDouble + raStr.substring(3, 5).toInt.toDouble / 60) / 24 * PI2;
        }
        def calcDecXyz(raStr: String, decStr: String): (Double, Double, Array[Double]) = {
          val ra = calcRa(raStr);
          val dec = decStr.toDouble / PI57;
          val z = Math.sin(dec);
          val xy = Math.cos(dec);
          val x = xy * Math.cos(ra);
          val y = xy * Math.sin(ra);
          val xyz = Array(x, y, z);
          (ra, dec, xyz);
        }
        val ConstellationPattern = "([-+.0-9]+)\\s+Constellation\\s+(Ecliptical\\s+)?(Galaxy\\s+)?(.+)".r;
        val StarsPattern2 = "([-+.0-9]+)\\s+Stars\\s+Bright\\s+([-+.0-9]+)\\s+(.+)".r;
        val StarsPattern = "([-+.0-9]+)\\s+Stars\\s+Bright\\s+(.+)".r;
        val CulminationPattern = "Cul\\s+(.+)".r;
        val LunchTimeContentPattern = "Lunch\\s+(.+)".r;
        val LunchTimeContentPattern2 = "([-+.0-9]+)\\s+Lunch\\s+(.+)".r;
        val AreaPattern = "([-+.0-9]+)\\sArea\\s(.+)".r;
        content0 match {
          case ConstellationPattern(decStr, eclipticalStr, galaxyStr, content) =>
            val (ra, dec, xyz) = calcDecXyz(raStr, decStr);
            val eclipticalFlag = eclipticalStr != null;
            val galaxyFlag = galaxyStr != null;
            constellationData = Constellations.Constellation(ra, dec, xyz, content, hashtags,
              eclipticalFlag, galaxyFlag) :: constellationData;
          case StarsPattern2(decStr, magStr, content) =>
            val (ra, dec, xyz) = calcDecXyz(raStr, decStr);
            val mag = magStr.toDouble;
            starData = (ra, xyz, content, hashtags) :: starData;
            starData2 = (ra, xyz, mag, content, hashtags) :: starData2;
          case StarsPattern(decStr, content) =>
            val (ra, dec, xyz) = calcDecXyz(raStr, decStr);
            starData = (ra, xyz, content, hashtags) :: starData;
          case CulminationPattern(content) =>
            val ra = calcRa(raStr);
            culminationContents = (ra, content, urlOpt, hashtags, starNames) :: culminationContents;
          case LunchTimeContentPattern(content) =>
            val ra = calcRa(raStr);
            lunchTimeContents = (ra, content, urlOpt, hashtags) :: lunchTimeContents;
          case LunchTimeContentPattern2(decStr, content) =>
            val (ra, dec, xyz) = calcDecXyz(raStr, decStr);
            val p = content.indexOf("座");
            val word = content.substring(0, p + 1);
            val hashtags2 = hashtags ::: word :: "星座" :: Nil;
            lunchTimeContents2 = Words.ConstellationLunchTimeContent(word, content, urlOpt, hashtags2, ra, dec, xyz) :: lunchTimeContents2;
          case AreaPattern(decStr, content) =>
            val lngH = raStr.substring(0, 2).toInt;
            val lngM = raStr.substring(3, 5).toInt / 20 * 20;
            val lat5 = decStr.toInt;
            val key = "%2dh%02dm,%3d".format(lngH, lngM, lat5);
            constellationsMap = constellationsMap ++ Map(key -> (content, hashtags));
          case _ => // nothing
        }
      }
    }
    source.close();
    (
      culminationContents.reverse.toIndexedSeq.sortBy(_._1),
      lunchTimeContents.reverse.toIndexedSeq.sortBy(_._1),
      lunchTimeContents2.reverse.toIndexedSeq.sortBy(_.ra),
      constellationData.reverse.toIndexedSeq.sortBy(_.ra),
      starData.reverse.toIndexedSeq.sortBy(_._1),
      starData2.reverse.toIndexedSeq.sortBy(_._1),
      constellationsMap,
    );
  }

  private def icrsToConstellation(lng: Double, lat: Double): (String, List[String]) = {
    val lng5 = (lng * PI57 / 5).toInt * 5;
    val lat5 = ((lat * PI57 + 90) / 5).toInt * 5 - 90;
    val key = "%2dh%02dm,%3d".format(lng5 / 15, lng5 % 15 * 4, lat5);
    val (cons, hashtags) = constellationsMap.getOrElse(key, ("", Nil));
    if (cons == "") {
      ("(ERROR %s)".format(key), hashtags);
    } else {
      (cons, hashtags);
    }
  }

  def icrsToConstellation(xyz: Array[Double]): (String, List[String]) = {
    val lng = VectorLib.xyzToLng(xyz);
    val lat = VectorLib.xyzToLat(xyz);
    icrsToConstellation(lng, lat);
  }

}

object Words {

  trait LunchTimeContent {
    def word: String;
    def content(day: Int): String;
    def urlOpt: Option[String];
    def hashtags: List[String];
  }

  case class ConstellationLunchTimeContent (word: String, contentTemplate: String,
    urlOpt: Option[String], hashtags: List[String],
    ra: Double, dec: Double, xyz: Array[Double]) extends LunchTimeContent {
    def content(day: Int): String = {
      val altHor = -0.90 / PI57;
      val time = startTime + day + 21.0 / 24;
      val tdb = TimeLib.mjdutcToTdb(time);
      val bpnMatrix = bpn.icrsToTrueEquatorialMatrix(tdb);
      val xyz2 = VectorLib.multiplyMV(bpnMatrix, xyz);
      val (azi, alt) = hcs.trueEquatorialXyzToAziAltFromUtc(xyz2, time);
      if (alt < altHor) {
        "";
      } else {
        val hcsStr = Hcs.aziAltToNaturalString(azi, alt);
        replaceSouthPattern(contentTemplate, hcsStr);
      }
    }

    private[this] def replaceSouthPattern(template: String, hcsStr: String): String = {
      val template2 = if (hcsStr.startsWith("南の") || hcsStr.startsWith("北の空高く")) {
        template.replace("[", "").replace("]", "");
      } else {
        var f: Boolean = true;
        var template2 = template;
        while (f) {
          f = false;
          val p1 = template2.indexOf("[");
          if (p1 >= 0) {
            val p2 = template2.indexOf("]", p1);
            if (p2 > p1) {
              template2 = template2.substring(0, p1) + template2.substring(p2 + 1);
              f = true;
            }
          }
        }
        template2;
      }
      if (template2.indexOf("%") >= 0) {
        template2.replace("%", hcsStr);
      } else {
        template2;
      }
    }
  }

}

class Words(constellationData: Constellations) {

  private def fetchLunchTimeContents(word: String): IndexedSeq[Words.LunchTimeContent] = {
    constellationData.lunchTimeContents2.filter(_.word == word);
  }

  private val defaultLimit = 25;
  private var history: Map[String, (List[Int], Int)] = Map.empty;
  private var contents: Map[String, IndexedSeq[Words.LunchTimeContent]] = Map.empty;

  val periodStart: Int = 104; // PERIOD

  def loadHistory(path: String): Unit = {
    history = Map.empty;
    val source = scala.io.Source.fromFile(path);
    source.getLines.foreach { line =>
      val cols = line.split(":");
      val word = cols(0);
      val (days, limit) = history.getOrElse(word, (Nil, defaultLimit));
      if (cols(1).startsWith("limit ")) {
        val limit = cols(1).substring(6).toInt;
        history = history + (word -> (days, limit));
      } else {
        val day = TimeLib.dateStringToDay(cols(1));
        if (day < periodStart) {
          history = history + (word -> (day :: days, limit));
        }
      }
    }
    source.close();
  }

  def saveHistory(path: String): Unit = {
    scala.util.Using(new java.io.PrintWriter(new java.io.FileOutputStream(path))) { writer =>
      history.keySet.toSeq.sorted.foreach { word =>
        val (days, limit) = history(word);
        {
          val line = "%s:limit %d".format(word, limit);
          writer.println(line);
        }
        days.reverse.foreach { day =>
          val line = "%s:%s".format(word, TimeLib.timeToDateString(startTime + day));
          writer.println(line);
        }
      }
    }
  }

  private def tweetPoint(word: String, day: Int): (IndexedSeq[Words.LunchTimeContent], Double) = {
    val cs = contents.getOrElse(word, {
      val cs = fetchLunchTimeContents(word);
      contents = contents + (word -> cs);
      cs;
    });
    if (cs.isEmpty) {
      (IndexedSeq.empty, 0.0);
    } else {
      val (days, limit) = history.getOrElse(word, (Nil, defaultLimit));
      if (days.isEmpty) {
        (cs, 10000.0);
      } else {
        val d = day - days.head;
        if (d < limit) {
          (cs, 0.0);
        } else {
          (cs, d.toDouble / limit);
        }
      }
    }
  }

  private def findRecentWord(day: Int): Option[(String, IndexedSeq[Words.LunchTimeContent])] = {
    val time = startTime + day + 12.0 / 24;
    val lastTweets = (-13 to 0).flatMap(i => getTweets(startTime + day + i).tweets);
    val lst = lastTweets.flatMap(_.starNames);
    if (lst.nonEmpty) {
      import Ordering.Double.IeeeOrdering;
      val lst2 = lst.map { w => (w, tweetPoint(w, day)) }.sortBy(- _._2._2);
      lst2.find { case (word, (contents, point)) =>
        point >= 1.0 && contents.filter(_.content(day) != "").nonEmpty;
      } match {
        case Some((word, (contents, point))) => Some((word, contents));
        case _ => None;
      }
    } else {
      None;
    }
  }

  private def addHistory(word: String, day: Int): Unit = {
    val (days, limit) = history.getOrElse(word, (Nil, defaultLimit));
    history = history + (word -> (day :: days, limit));
  }

  def putTweetWord(day: Int): Unit = {
    findRecentWord(day) match {
      case None => ;
      case Some((word, contents)) =>
        val time1 = startTime + day + (12.0 + 50.0 / 60) / 24;
        var inc: Int = 0;
        contents.foreach { content =>
          val c = content.content(day + inc / 4);
          if (c != "") {
            val msg = c + content.hashtags.map(" #" + _).mkString("");
            putTweet(time1 + 1.0 * (inc % 4) / 24, msg, content.urlOpt);
            if (inc == 0) {
              addHistory(word, day + inc / 4);
            }
            inc += 1;
          }
        }
    }
  }

}

val jplData = new JplData(jplDataPath);
val hcs = new Hcs(tokyoLng, tokyoLat);
val constellationData = new Constellations();
val words = new Words(constellationData);
words.loadHistory(wordHistoryPath);

//==============================================================================
// イベント計算
//==============================================================================

val sunsetTimesData: IndexedSeq[(Double, Double, Array[Double])] = { // time, tdb, bpnMatrix
  val altHor = -0.90 / PI57;
  (0 until period).map { d =>
    val time = Lib2.findCrossingBoundaryTime(altHor, true, false, startTime + d + 16.0 / 24.0, 24 * 6, 4 * 6) { time =>
      val tdb = TimeLib.mjdutcToTdb(time);
      val sun = jplData.calcPlanetFromEarth(tdb, JplData.Sun);
      val bpnMatrix = bpn.icrsToTrueEquatorialMatrix(tdb);
      val sun2 = VectorLib.multiplyMV(bpnMatrix, sun);
      val (azi, alt) = hcs.trueEquatorialXyzToAziAltFromUtc(sun2, time);
      alt;
    }
    val tdb = TimeLib.mjdutcToTdb(time);
    val bpnMatrix = bpn.icrsToTrueEquatorialMatrix(tdb);
    (time, tdb, bpnMatrix);
  }
}

def sunsetTimes(day: Int): Double = sunsetTimesData(day)._1;

def calcPlanetOnSunsetTime(day: Int, targetPlanet: JplData.TargetPlanet): (Double, Double) = { // azi, alt
  val (time, tdb, bpnMatrix) = sunsetTimesData(day);
  val xyz = jplData.calcPlanetFromEarth(tdb, targetPlanet);
  val xyz2 = VectorLib.multiplyMV(bpnMatrix, xyz);
  hcs.trueEquatorialXyzToAziAltFromUtc(xyz2, time);
}

def calcPlanetXyzAziAlt(time: Double, targetPlanet: JplData.TargetPlanet): (Array[Double], Double, Double) = {
  val tdb = TimeLib.mjdutcToTdb(time);
  val ut1 = time; // 近似的
  val bpnMatrix = bpn.icrsToTrueEquatorialMatrix(tdb);
  val xyz = jplData.calcPlanetFromEarth(tdb, targetPlanet);
  val xyz2 = VectorLib.multiplyMV(bpnMatrix, xyz);
  val azialt = hcs.trueEquatorialXyzToAziAltFromUtc(xyz2, time);
  (xyz, azialt._1, azialt._2);
}

def calcRiseSetCulmination(targetPlanet: JplData.TargetPlanet, altHor: Double):
  IndexedSeq[(Double, Double, Double, Double)] = { // 出, 南中, 没, 南中高度
  val riseSet = MathLib.findMaxMinCrossingListContinuous(startTime, endTime, 0.25, 24 * 6) { time =>
    calcPlanetXyzAziAlt(time, targetPlanet)._3 - altHor;
  }.filter(_._2 % 2 == 0);
  var result: List[(Double, Double, Double, Double)] = Nil;
  (0 until (riseSet.size - 1)).foreach { i =>
    val (time0, flag0) = riseSet(i);
    if (flag0 == 0) {
      val (time2, flag2) = riseSet(i + 1);
      val time1 = Lib2.findCrossingBoundaryTime(0.0, false, false,
        time0, 24 * 6, ((time2 - time0) * 24 * 6).toInt) { time =>
        calcPlanetXyzAziAlt(time, JplData.Moon)._2 - PI;
      }
      val xyz = calcPlanetXyzAziAlt(time1, targetPlanet)._1;
      val x = xyz(0);
      val y = xyz(1);
      val z = xyz(2);
      val xy = Math.sqrt(x * x + y * y);
      val altm = Math.atan2(z, xy);
      val alt = altm + PI5 - tokyoLat;
      result = (time0, time1, time2, alt) :: result;
    }
  }
  result.reverse.toIndexedSeq;
}

val moonRiseSetTimesData: IndexedSeq[(Double, Double, Double, Double, Int)] = {
  val altHor = +0.4 / PI57; // だいたい
  val data: IndexedSeq[(Double, Double, Double, Double)] = calcRiseSetCulmination(JplData.Moon, altHor);
  val data2 = Array.fill[Int](data.size)(-1);
  (0 until 8).foreach { phase =>
    val q = phase.toDouble / 8 * PI2;
    MathLib.findMaxMinListDiscrete(0, data.size, 7) { p =>
      Math.abs(MathLib.circleDiff1(calcMoonPhase(data(p)._2) / 28 * PI2, q));
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
      putTweet(time2, riseset);
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
  val sun = jplData.calcPlanetFromEarth(tdb, JplData.Sun);
  val moon = jplData.calcPlanetFromEarth(tdb, JplData.Moon);
  val bpnMatrix = bpn.icrsToTrueEclipticMatrix(tdb);
  val sun2 = VectorLib.multiplyMV(bpnMatrix, sun);
  val moon2 = VectorLib.multiplyMV(bpnMatrix, moon);
  val sunLng = VectorLib.xyzToLng(sun2);
  val moonLng = VectorLib.xyzToLng(moon2);
  VectorLib.calcLngDiff(moonLng, sunLng);
}
def calcMoonPhase(time: Double): Double = {
  calcMoonLng(time) / PI2 * 28.0;
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
  MathLib.findCyclicPhaseListContinuous(8, startTime, endTime, 3.0, 24)(calcMoonLng);
}

val innerPlanets = IndexedSeq(
  ("金星", JplData.Venus, 5),
  ("水星", JplData.Mercury, 3)
);
val innerPlanetsSunsetAziAltList: IndexedSeq[IndexedSeq[(Double, Double)]] = innerPlanets.map { planet =>
  (0 until period).map { day =>
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
    period + 700;
  } else {
    period + 300;
  }
  (0 until period2).map { day =>
    val time = startTime + day + 21.0 / 24.0;
    val t = calcPlanetXyzAziAlt(time, planet._2);
    (t._2, t._3);
  }
}

//==============================================================================

def isLunchTime(time: Double): Boolean = {
  val day = (time - startTime).toInt;
  val s = sunsetTimes(day);
  //time - time.toInt >= 3.0 / 24 && time - time.toInt < 4.0 / 24; // 12時～13時
  time - time.toInt >= 1.0 / 24 && time - time.toInt < 7.0 / 24; // 10時～16時
}

def isNightTime0(time: Double): Boolean = {
  val day = (time - startTime).toInt;
  val s = sunsetTimes(day);
  //time - time.toInt >= 5.0 / 24 && time - time.toInt < 15.0 / 24;
  time > s && time - time.toInt < 14.0 / 24;
}

def isNightTime2(time: Double): Boolean = {
  val day = (time - startTime).toInt;
  val s = sunsetTimes(day);
  time >= s && time - time.toInt < 16.0 / 24;
}

def isNightTime3(time: Double): Boolean = {
  val day = (time - startTime).toInt;
  val s = sunsetTimes(day);
  time - time.toInt >= 11.0 / 24 && time - time.toInt < 14.1 / 24;
}

def getConjunctionTweetTime(time: Double, xyz: Array[Double]): Option[Double] = {
  val altThres0 = 10 / PI57;
  val sun_distanceThres = 20.0 / PI57;
  val utc = time;
  val ut1 = time; // 近似的
  val tdb = TimeLib.mjdutcToTdb(utc);
  val sun_xyz = jplData.calcPlanetFromEarth(tdb, JplData.Sun);
  val sun_distance = VectorLib.angularDistance(sun_xyz, xyz);
  if (sun_distance < sun_distanceThres) {
    None;
  } else {
    val time21 = time.toInt + 0.5;
    val tdb21 = TimeLib.mjdutcToTdb(time21);
    val bpnMatrix = bpn.icrsToTrueEquatorialMatrix(tdb21);
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
        val day = (time - startTime).toInt;
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

var _tweets: Map[String, DateTweets] = (0 until period).map { day =>
  val date = TimeLib.timeToDateString(startTime + day);
  (date, DateTweets());
}.toMap;

def getTweets(time: Double): DateTweets = {
  val date = TimeLib.timeToDateString(time);
  _tweets.getOrElse(date, DateTweets());
}

def putTweet(tc: TweetContent): Unit = {
  val date = TimeLib.timeToDateString(tc.time);
  if (_tweets.contains(date)) {
    val tw = _tweets(date);
    _tweets = _tweets.updated(date, tw.added(tc));
  }
}

def putTweet(time: Double, msg: String): Unit = {
  putTweet(LegacyTweetContent(time, msg, None, Nil));
}

def putTweet(time: Double, msg: String, urlOpt: Option[String]): Unit = {
  putTweet(LegacyTweetContent(time, msg, urlOpt, Nil));
}

def putTweet(time: Double, msg: String, starNames: List[String]): Unit = {
  putTweet(LegacyTweetContent(time, msg, None, starNames));
}

def removeTweet(time: Double, message: String): Unit = {
  val date = TimeLib.timeToDateString(time);
  if (_tweets.contains(date)) {
    val tw = _tweets(date);
    _tweets = _tweets.updated(date, tw.removed(time, message));
  }
}

//==============================================================================
// 祝日・記念日
//==============================================================================

{
  val yearStart = TimeLib.timeToDateTimeString(startTime).substring(0, 4).toInt;
  val yearEnd = TimeLib.timeToDateTimeString(endTime).substring(0, 4).toInt;
  val source = scala.io.Source.fromFile(holidayDataPath);
  source.getLines.foreach { line =>
    if (line != "" && !line.startsWith("#")) {
      val cols = line.split(" ", 2);
      if (cols(0).startsWith("YYYY-")) {
        val msg = cols(1);
        (yearStart to yearEnd).foreach { year =>
          val timeStr = year.toString + cols(0).substring(4);
          val time = TimeLib.stringToModifiedJulianDay(timeStr + ":00+09:00");
          putTweet(time, msg);
        }
      } else {
        val time = TimeLib.stringToModifiedJulianDay(cols(0) + ":00+09:00");
        val msg = cols(1);
        putTweet(time, msg);
      }
    }
  }
  source.close();
}

//==============================================================================
// 特定の時間に発生する天文の現象
//==============================================================================

def calcSunLng(time: Double): Double = {
  val utc = time;
  val tdb = TimeLib.mjdutcToTdb(utc);
  val sun = jplData.calcPlanetFromEarth(tdb, JplData.Sun);
  val bpnMatrix = bpn.icrsToTrueEclipticMatrix(tdb);
  val sun2 = VectorLib.multiplyMV(bpnMatrix, sun);
  val sunLng = VectorLib.xyzToLng(sun2);
  sunLng;
}
def calcSunLng2(time: Double): Double = {
  val utc = time;
  val tdb = TimeLib.mjdutcToTdb(utc);
  val sun = jplData.calcPlanetFromEarth(tdb, JplData.Sun);
  val bpnMatrix = bpn.icrsToMeanEclipticMatrix2000;
  val sun2 = VectorLib.multiplyMV(bpnMatrix, sun);
  val sunLng = VectorLib.xyzToLng(sun2);
  sunLng;
}

//--------------------------------------
// 24節気
//--------------------------------------

{
  val sunPhaseTerms = MathLib.findCyclicPhaseListContinuous(24, startTime, endTime, 3.0, 24) { time =>
    calcSunLng(time);
  }
  val termStrs = IndexedSeq(
    "春分", "清明", "穀雨", "立夏", "小満", "芒種", "夏至", "小暑", "大暑", "立秋", "処暑", "白露",
    "秋分", "寒露", "霜降", "立冬", "小雪", "大雪", "冬至", "小寒", "大寒", "立春", "雨水", "啓蟄",
  );
  sunPhaseTerms.foreach { case (time, term) =>
    if (term % 6 == 0) {
      putTweet(TimeLib.floor(time, 24) + 1.0 / (24 * 4), "%s。太陽の黄経が%d°です".format(termStrs(term), term * 15));
    } else {
      putTweet(TimeLib.floor(time, 24) + 1.0 / (24 * 4), "二十四節気の%s。太陽の黄経が%d°です".format(termStrs(term), term * 15));
    }
  }
}

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
          val lng = cols(0).toDouble / PI57;
          data1 = (lng, cols(1)) :: data1;
        }
      }
    }
    source.close();
    (data1.reverse.sortBy(_._1).toIndexedSeq, data2.reverse.toIndexedSeq);
  }

  var day: Int = 0;
  var index: Int = {
    val time = startTime + day + 0.5;
    val lng = calcSunLng2(time);
    val index = meteorShowerData._1.indexWhere(_._1 > lng);
    if (index < 0) {
      0;
    } else {
      index;
    }
  }
  while (day < period) {
    val time = startTime + day + 0.5;
    val lng = calcSunLng2(time);
    val d = MathLib.circleAdd(lng, -meteorShowerData._1(index)._1);
    if (d >= 0) {
      val time1 = time - 1.0 + 20.0 / 60 / 24;
      val name = meteorShowerData._1(index)._2;
      val msg = "%sは今夜が極大。%s".format(name, moonRiseSetForMeteorShower(time1));
      putTweet(time1, msg);
      index += 1;
      if (index == meteorShowerData._1.size) {
        index = 0;
      }
    } else {
      day += 1;
    }
  }

  meteorShowerData._2.foreach { case (time, msg) =>
    putTweet(time, msg);
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

  val altThres0 = 10 / PI57;

  def calcMoonConstellation(time: Double): Option[(String, List[String])] = {
    if (isNightTime2(time)) {
      val tdb = TimeLib.mjdutcToTdb(time);
      val xyz = jplData.calcPlanetFromEarth(tdb, JplData.Moon);
      val bpnMatrix = bpn.icrsToTrueEquatorialMatrix(tdb);
      val xyz2 = VectorLib.multiplyMV(bpnMatrix, xyz);
      val (azi, alt) = hcs.trueEquatorialXyzToAziAltFromUtc(xyz2, time);
      if (alt >= altThres0) {
        val (cons, hashtags) = constellationData.icrsToConstellation(xyz);
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
      val moon = jplData.calcPlanetFromEarth(tdb, JplData.Moon);
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
    }.foreach(putTweet);

    val midAltData = moonRiseSetTimesData.filter(_._5 == 4).map(t => (t._2, t._4));
    (0 until midAltData.size).foreach { i =>
      val (time, alt) = midAltData(i);
      if (i > 0 && i < midAltData.size - 1) {
        if (midAltData(i - 1)._2 <= alt && midAltData(i + 1)._2 < alt) {
          putTweet(time, "満月が南中(高度%.0f°)。冬の満月は空高く上り、今日はこの冬でもっとも天頂に近い満月です".format(alt * PI57));
        } else if (midAltData(i - 1)._2 >= alt && midAltData(i + 1)._2 > alt) {
          putTweet(time, "満月が南中(高度%.0f°)。夏の満月は空低く、今日はこの夏でもっとも低い満月です".format(alt * PI57));
        } else {
          putTweet(time, "満月が南中(高度%.0f°)".format(alt * PI57));
        }
      } else {
        putTweet(time, "満月が南中(高度%.0f°)".format(alt * PI57));
      }
    }
  }

  moonPhaseTerms.filter(_._2 != 4).filter(_._2 % 2 == 0).map { case (time, term) =>
    MoonPhaseTermTweetContent(time, term, 0, calcMoonConstellation(time));
  }.foreach(putTweet);
}

//--------------------------------------
// 近日点・遠日点
//--------------------------------------

MathLib.findMaxMinListContinuous(startTime, endTime, 30, 24) { time =>
  val utc = time;
  val tdb = TimeLib.mjdutcToTdb(utc);
  val sun = jplData.calcPlanetFromEarth(tdb, JplData.Sun);
  VectorLib.distance(sun);
}.foreach { case (time, flag) =>
  val s = if (flag < 0) "近日点" else "遠日点";
  putTweet(TimeLib.round(time, 24) - 1.0 / (24 * 4), "地球が%s通過".format(s));
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
    val sun = jplData.calcPlanetFromEarth(tdb, JplData.Sun);
    val planet = jplData.calcPlanetFromEarth(tdb, targetPlanet);
    val bpnMatrix = bpn.icrsToTrueEclipticMatrix(tdb);
    val sun2 = VectorLib.multiplyMV(bpnMatrix, sun);
    val planet2 = VectorLib.multiplyMV(bpnMatrix, planet);
    val sunLng = VectorLib.xyzToLng(sun2);
    val planetLng = VectorLib.xyzToLng(planet2);
    val d1 = VectorLib.calcLngDiff(planetLng, sunLng);
    if (d1 >= PI) d1 - PI2 else d1;
  }

  {
    val (planetName, targetPlanet) = ("水星", JplData.Mercury);
    MathLib.findMaxMinCrossingListContinuous(startTime, endTime, 10.0, 24) { time =>
      calcInnerPlanetLngEc(time, targetPlanet);
    }.map { case (time, term) =>
      if (term == 0) {
        putTweet(PlanetAstronomyTweetContent(TimeLib.floor(time, 24) + 1.0 / (24 * 4),
          "水星が外合(黄経基準)", planetName));
      } else if (term == 1) {
        putTweet(PlanetAstronomyTweetContent(TimeLib.round(time, 24) - 1.0 / (24 * 4),
          "水星が東方最大離角(黄経基準)\uD83C\uDF13", planetName));
      } else if (term == 2) {
        putTweet(PlanetAstronomyTweetContent(TimeLib.floor(time, 24) + 1.0 / (24 * 4),
          "水星が内合(黄経基準)", planetName));
      } else {
        putTweet(PlanetAstronomyTweetContent(TimeLib.round(time, 24) - 1.0 / (24 * 4),
          "水星が西方最大離角(黄経基準)\uD83C\uDF17", planetName));
      }
    }
  }

  {
    val (planetName, targetPlanet) = ("金星", JplData.Venus);
    MathLib.findMaxMinCrossingListContinuous(startTime, endTime, 10.0, 24) { time =>
      calcInnerPlanetLngEc(time, targetPlanet);
    }.map { case (time, term) =>
      if (term == 0) {
        putTweet(PlanetAstronomyTweetContent(TimeLib.floor(time, 24) + 1.0 / (24 * 4),
          "金星が外合(黄経基準)。数か月後に夕方の西の空に現れます", planetName));
      } else if (term == 1) {
        putTweet(PlanetAstronomyTweetContent(TimeLib.round(time, 24) - 1.0 / (24 * 4),
          "金星が東方最大離角(黄経基準)\uD83C\uDF13。宵の明星として夕方に西の空にいます", planetName));
      } else if (term == 2) {
        putTweet(PlanetAstronomyTweetContent(TimeLib.floor(time, 24) + 1.0 / (24 * 4),
          "金星が内合(黄経基準)。数週間後に明け方の東の空に現れます", planetName));
      } else {
        putTweet(PlanetAstronomyTweetContent(TimeLib.round(time, 24) - 1.0 / (24 * 4),
          "金星が西方最大離角(黄経基準)\uD83C\uDF17。明けの明星として明け方に東の空にいます", planetName));
      }
    }
  }
}

{
  def calcOuterPlanetLngEq(time: Double, targetPlanet: JplData.TargetPlanet): Double = {
    val utc = time;
    val tdb = TimeLib.mjdutcToTdb(utc);
    val sun = jplData.calcPlanetFromEarth(tdb, JplData.Sun);
    val planet = jplData.calcPlanetFromEarth(tdb, targetPlanet);
    val bpnMatrix = bpn.icrsToTrueEquatorialMatrix(tdb);
    val sun2 = VectorLib.multiplyMV(bpnMatrix, sun);
    val planet2 = VectorLib.multiplyMV(bpnMatrix, planet);
    val sunLng = VectorLib.xyzToLng(sun2);
    val planetLng = VectorLib.xyzToLng(planet2);
    val d1 = VectorLib.calcLngDiff(planetLng, sunLng);
    d1;
  }

  val planetPhases1: List[(Double, String, String, Boolean, Option[Array[Double]], Int)] = outerPlanets.zipWithIndex.toList.
  flatMap { case ((planetName, targetPlanet, _), pi) =>
    MathLib.findCyclicPhaseListContinuous(4, startTime, endTime, 30, 24) { time =>
      PI2 - calcOuterPlanetLngEq(time, targetPlanet);
    }.map { case (time, term) =>
      val utc = time;
      val tdb = TimeLib.mjdutcToTdb(utc);
      val xyz = jplData.calcPlanetFromEarth(tdb, targetPlanet);
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

  val altThres = 15 / PI57;

  planetPhases1.foreach { case (time, planetName, content, timeFlag, xyzOpt, pi) =>
    val time2 = if (timeFlag) {
      TimeLib.round(time, 24) - 1.0 / (24 * 4);
    } else {
      TimeLib.floor(time, 24) + 1.0 / (24 * 4);
    }
    val nextSeasonMsg = {
      val day = (time - startTime).toInt;
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
        val (cons, hashtags) = constellationData.icrsToConstellation(xyz);
        PlanetAstronomyTweetContent(time2, "%sが%s。%sにいます%s".format(planetName, content, cons, nextSeasonMsg) +
          hashtags.map(" #" + _).mkString, planetName);
    }
    putTweet(tweet);
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
  private[this] val azi360: Int = (azi * PI57 + 0.5).toInt;
  private[this] val alt360: Int = (alt * PI57 + 0.5).toInt;
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
  azi: Double, alt: Double, isIncreasing: Boolean, isMax: Boolean) extends OnSunsetTweetContent {
  def alt360: Int = (alt * PI57 + 0.5).toInt;
  def message: String = if (isMax) {
    "%sは日没時最大高度で西の空高度約%d°".format(planetName, alt360);
  } else if (isIncreasing) {
    "%sは日没時の高度を徐々に上げ、西の空高度約%d°にいます".format(planetName, alt360);
  } else {
    "%sは日没時の高度を徐々に下げ、西の空高度約%d°にいます".format(planetName, alt360);
  }
  def message2: String = if (isMax) {
    "%sは日没時最大高度で西の空高度約%d°です".format(planetName, alt360);
  } else if (isIncreasing) {
    "%sは日没時の高度を徐々に上げ、西の空高度約%d°にいます".format(planetName, alt360);
  } else {
    "%sは日没時の高度を徐々に下げ、西の空高度約%d°にいます".format(planetName, alt360);
  }
  def message3: String = if (isMax) {
    "%sは日没時最大高度で西の空高度約%d°です".format(planetName, alt360);
  } else if (isIncreasing) {
    "%sは日没時の高度を徐々に上げ、約%d°にいます".format(planetName, alt360);
  } else {
    "%sは日没時の高度を徐々に下げ、約%d°にいます".format(planetName, alt360);
  }
  def hashtags: List[String] = List(planetName);
  def starNames: List[String] = List(planetName);
}
case class SunsetStarTweetContent(day: Int, starName: String,
  azi: Double, alt: Double) extends OnSunsetTweetContent {
  def alt360: Int = (alt * PI57 + 0.5).toInt;
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
    MathLib.findMaxMinListDiscrete(0, period, 15) { day =>
      val (azi, alt) = calcPlanetOnSunsetTime(day, p._2);
      alt;
    }.foreach { case (day, flag) =>
      if (flag > 0) {
        val (azi, alt) = calcPlanetOnSunsetTime(day, p._2);
        putTweet(SunsetPlanetTweetContent(day, p._1, azi, alt, true, true));
      }
    }
  }
}

// 新月直後の月
{
  val aziThres0 = 200 / PI57;
  val aziThres1 = 315 / PI57;
  val altThres0 = 10 / PI57;

  (0 until moonPhaseTerms.size).foreach { i =>
    val (moonPhaseTime, term) = moonPhaseTerms(i);
    if (term == 0 && moonPhaseTime + 4 < endTime) {
      (0 until 4).foreach { d =>
        val day = (moonPhaseTime + d - startTime).toInt;
        val (azi, alt) = calcPlanetOnSunsetTime(day, JplData.Moon);
        if (azi >= aziThres0 && azi <= aziThres1 && alt >= altThres0) {
          val day = (moonPhaseTime + d - startTime).toInt;
          val (azi, alt) = calcPlanetOnSunsetTime(day, JplData.Moon);
          putTweet(SunsetMoonTweetContent(day, azi, alt));
        }
      }
    }
  }
}

// 水曜・金曜
{
  val aziThres0 = 200 / PI57;
  val aziThres1 = 315 / PI57;
  val altThres0 = 10 / PI57;

  val planets = IndexedSeq(("金星", JplData.Venus, 5), ("水星", JplData.Mercury, 3));
  innerPlanets.zipWithIndex.foreach { case (planet, pi) =>
    (1 until period).foreach { day =>
      val wday = TimeLib.wday(startTime + day);
      val sunsetTweets = getTweets(startTime + day).sunsetTweets;
      if (wday == planet._3 && !sunsetTweets.exists(_.starNames.contains(planet._1))) {
        val (azi, alt) = innerPlanetsSunsetAziAltList(pi)(day);
        if (alt >= altThres0) {
          val prevAlt = innerPlanetsSunsetAziAltList(pi)(day - 1)._2;
          val isIncreasing = (alt >= prevAlt);
          putTweet(SunsetPlanetTweetContent(day, planet._1, azi, alt, isIncreasing, false));
        }
      }
    }
  }
}

// 日没ツイートがある場合に他の天体のツイートも追加
{
  val aziThres0 = 200 / PI57;
  val aziThres1 = 315 / PI57;
  val altThres0 = 10 / PI57;

  val planets = IndexedSeq(("金星", JplData.Venus), ("水星", JplData.Mercury));
  val planets2 = IndexedSeq(("火星", JplData.Mars), ("木星", JplData.Jupiter), ("土星", JplData.Saturn));
  (1 until period).foreach { day =>
    val sunsetTweets = getTweets(startTime + day).sunsetTweets;
    if (sunsetTweets.nonEmpty) {
      planets.foreach { p =>
        if (!sunsetTweets.exists(_.starNames.contains(p._1))) {
          val (azi, alt) = calcPlanetOnSunsetTime(day, p._2);
          if (azi >= aziThres0 && azi <= aziThres1 && alt >= altThres0) {
            val (_, prevAlt) = calcPlanetOnSunsetTime(day - 1, p._2);
            val isIncreasing = (alt >= prevAlt);
            putTweet(SunsetPlanetTweetContent(day, p._1, azi, alt, isIncreasing, false));
          }
        }
      };
      {
        val p = ("月", JplData.Moon);
        if (!sunsetTweets.exists(_.starNames.contains(p._1))) {
          val (azi, alt) = calcPlanetOnSunsetTime(day, p._2);
          if (azi >= aziThres0 && azi <= aziThres1 && alt >= altThres0) {
            putTweet(SunsetMoonTweetContent(day, azi, alt));
          }
        }
      }
      planets2.foreach { p =>
          val (azi, alt) = calcPlanetOnSunsetTime(day, p._2);
          if (azi >= aziThres0 && azi <= aziThres1 && alt >= altThres0) {
            putTweet(SunsetStarTweetContent(day, p._1, azi, alt));
          }
      };
    }
  }
}

def calcNextSeason(aziAltList: IndexedSeq[(Double, Double)], day: Int): (Int, Int) = {
  val altThres = 15 / PI57;
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
  val altThres = 15 / PI57;

  innerPlanets.zipWithIndex.foreach { case (planet, pi) =>
    (1 until period).foreach { day =>
      if ((startTime + day).toInt % 2 == 0) {
        val wday = TimeLib.wday(startTime + day);
        val sunsetTweets = getTweets(startTime + day).sunsetTweets;
        if (wday == planet._3 && !sunsetTweets.exists(_.starNames.contains(planet._1))) {
          val (azi, alt) = innerPlanetsSunsetAziAltList(pi)(day);
          if (alt < altThres) {
            val time = startTime + day + 12.5 / 24;
            val (p1, p2) = calcNextSeason(innerPlanetsSunsetAziAltList(pi), day);
            putTweet(SunsetNextPlanetTweetContent(time, planet._1, p1, p2));
          }
        }
      }
    }
  }
}

// 日没時間
{
  MathLib.getMaxMinUpDownFlagListDiscrete(0, period, 90) { day =>
    sunsetTimes(day) - startTime - day;
  }.zipWithIndex.foreach { case (flag, day) =>
    val wday = TimeLib.wday(startTime + day);
    if (flag == 1 || flag == 3 || wday == 0) {
      putTweet(SunsetTweetContent(day, flag));
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
    val distance360 = distance * PI57;
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
  val altThres0 = 10 / PI57;
  val distanceThres = 3.0 / PI57;
  val distanceThresMoon = 6.0 / PI57;

  def calcClosestMoon(slowStarName: String, fastStarName: String, hashtags: List[String],
    slowStarXyzFunc: Double => Array[Double],
    fastStarXyzFunc: Double => Array[Double]): Unit = {
    MathLib.findMaxMinListContinuous(startTime, endTime, 10.0, 24 * 6) { time =>
      val utc = time;
      val tdb = TimeLib.mjdutcToTdb(utc);
      val xyz_s = slowStarXyzFunc(tdb);
      val xyz_f = fastStarXyzFunc(tdb);
      VectorLib.angularDistance1(xyz_s, xyz_f);
    }.foreach { case (time, flag) =>
      if (flag > 0 && isNightTime2(time)) {
        val tdb = TimeLib.mjdutcToTdb(time);
        val xyz_f = fastStarXyzFunc(tdb);
        val bpnMatrix = bpn.icrsToTrueEquatorialMatrix(tdb);
        val xyz_f2 = VectorLib.multiplyMV(bpnMatrix, xyz_f);
        val (azi, alt) = hcs.trueEquatorialXyzToAziAltFromUtc(xyz_f2, time);
        if (alt >= altThres0) {
          val xyz_s = slowStarXyzFunc(tdb);
          val distance = VectorLib.angularDistance(xyz_s, xyz_f);
          if (distance < distanceThresMoon) {
            putTweet(CloseStarsTweetContent(time, 24 * 6, slowStarName, fastStarName, distance, hashtags));
          }
        }
      }
    }
  }
  def calcClosest2(slowStarName: String, fastStarName: String, hashtags: List[String],
    slowStarXyzFunc: Double => Array[Double],
    fastStarXyzFunc: Double => Array[Double]): Unit = {
    MathLib.findMaxMinListContinuous(startTime, endTime, 10.0, 24) { time =>
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
              putTweet(CloseStarsTweetContent(postTime, 24, slowStarName, fastStarName, distance, hashtags));
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
    val lng = (lngStr.substring(0, 2).toInt.toDouble + lngStr.substring(3, 5).toInt.toDouble / 60) / 24 * PI2;
    val lat = (latStr.substring(0, 3).toInt.toDouble + latStr.substring(4, 6).toInt.toDouble / 60) / 360 * PI2;
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
          jplData.calcPlanetFromEarth(tdb, star2._2);
        });
    }
  }
  stars2.foreach { star2 =>
    stars1.foreach { star1 =>
      calcClosestMoon(star1._1, star2._1, List(star1._1),
        { tdb: Double =>
          jplData.calcPlanetFromEarth(tdb, star1._2);
        },
        { tdb: Double =>
          jplData.calcPlanetFromEarth(tdb, star2._2);
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
          jplData.calcPlanetFromEarth(tdb, star1._2);
        });
    }
  }
  (0 until stars1.size).foreach { i =>
    ((i + 1) until stars1.size).foreach { j =>
      val star1 = stars1(i);
      val star0 = stars1(j);
      calcClosest2(star0._1, star1._1, List(star1._1, star0._1),
        { tdb: Double =>
          jplData.calcPlanetFromEarth(tdb, star0._2);
        },
        { tdb: Double =>
          jplData.calcPlanetFromEarth(tdb, star1._2);
        });
    }
  }
}

//==============================================================================
// 21時・23時の月
//==============================================================================

{
  val altThres = 10 / PI57;
  (0 until period).foreach { d =>
    val time = startTime + d + 21.0 / 24.0 - 1.0 / (24 * 6);
    if (!getTweets(time).tweets.filter(tc => isNightTime0(tc.time)).flatMap(_.starNames).contains("月")) {
      {
        val (xyz, azi, alt) = calcPlanetXyzAziAlt(time, JplData.Moon);
        if (alt >= altThres) {
          Some((time, xyz, azi, alt));
        } else {
          val time = startTime + d + 23.0 / 24.0;
          val (xyz, azi, alt) = calcPlanetXyzAziAlt(time, JplData.Moon);
          if (alt >= altThres) {
            Some((time, xyz, azi, alt));
          } else {
            None;
          }
        }
      } match {
        case None => ;
        case Some((time, xyz, azi, alt)) =>
          val (cons, hashtags) = constellationData.icrsToConstellation(xyz);
          val hcsStr = Hcs.aziAltToNaturalString(azi, alt);
          if (azi >= PI) {
            putTweet(time, "月は%s、%sにいます。%s".format(hcsStr, cons, calcMoonPhaseString(time)) +
              hashtags.map(" #" + _).mkString);
          } else {
            putTweet(time, "月は%s、%sにいます".format(hcsStr, cons) +
              hashtags.map(" #" + _).mkString);
          }
      }
    }
  }
}

//==============================================================================
// 21時・23時の惑星
//==============================================================================

{
  val altThres = 15 / PI57;
  (0 until period).foreach { d =>
    val time = startTime + d + 21.0 / 24.0;
    val wday = TimeLib.wday(time);
    val (planetName, targetPlanet) = if (wday == 2) {
      ("火星", Some(JplData.Mars));
    } else if (wday == 4) {
      ("木星", Some(JplData.Jupiter));
    } else if (wday == 6) {
      ("土星", Some(JplData.Saturn));
    } else {
      ("", None);
    }
    targetPlanet match {
      case None => ;
      case Some(targetPlanet) =>
        {
          val (xyz, azi, alt) = calcPlanetXyzAziAlt(time, targetPlanet);
          if (alt >= altThres) {
            Some((time, xyz, azi, alt));
          } else {
            val time = startTime + d + 23.0 / 24.0;
            val (xyz, azi, alt) = calcPlanetXyzAziAlt(time, targetPlanet);
            if (alt >= altThres) {
              Some((time, xyz, azi, alt));
            } else {
              None;
            }
          }
        } match {
          case None => ;
          case Some((time, xyz, azi, alt)) =>
            val (cons, hashtags) = constellationData.icrsToConstellation(xyz);
            val hcsStr = Hcs.aziAltToNaturalString(azi, alt);
            putTweet(time - 1.0 / (24 * 4), "%sは%s、%sにいます #%s".format(planetName, hcsStr, cons, planetName) +
              hashtags.map(" #" + _).mkString);
        }
    }
  }
}

// 惑星が見えない場合
{
  val altThres = 15 / PI57;

  outerPlanets.zipWithIndex.foreach { case (planet, pi) =>
    (0 until period).foreach { day =>
      if ((startTime + day).toInt % 2 == 0) {
        val wday = TimeLib.wday(startTime + day);
        val tweets = getTweets(startTime + day).tweets;
        if (wday == planet._3 && !tweets.exists(_.starNames.contains(planet._1))) {
          if (!(-6 to +6).exists { i =>
            getTweets(startTime + day + i).tweets.exists {
              case PlanetAstronomyTweetContent(_, _, planetName) if (planetName == planet._1) => true;
              case _ => false;
            }
          }) {
            val (azi, alt) = outerPlanetsNightAziAltList(pi)(day);
            if (alt < altThres) {
              val time = startTime + day + 12.5 / 24;
              val (p1, p2) = calcNextSeason(outerPlanetsNightAziAltList(pi), day);
              putTweet(NightNextPlanetTweetContent(time, planet._1, p1, p2));
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

{
  val diffData: List[(Double, Int, String)] = {
    var diffData: List[(Double, Int, String)] = Nil;
    val source = scala.io.Source.fromFile(diffDataPath);
    source.getLines.foreach { line =>
      if (line != "" && !line.startsWith("#")) {
        val cols = line.split(" ", 2);
        val time = TimeLib.stringToModifiedJulianDay(cols(0) + ":00+09:00");
        val sign = cols(1).charAt(0);
        val content = cols(1).substring(1);
        if (sign == '-') {
          diffData = (time, -1, content) :: diffData;
        } else if (sign == '+') {
          diffData = (time, +1, content) :: diffData;
        }
      }
    }
    source.close();
    import Ordering.Double.IeeeOrdering;
    diffData.reverse.sortBy(_._1);
  }
  diffData.filter(_._2 < 0).foreach { d =>
    removeTweet(d._1, d._3);
  }
  diffData.filter(_._2 > 0).foreach { d =>
    putTweet(d._1, d._3);
  }
}

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
    val altThres = 30.0 / PI57;
    val time = startTime + day + 21.0 / 24;
    val tdb = TimeLib.mjdutcToTdb(time);
    val bpnMatrix = bpn.icrsToTrueEquatorialMatrix(tdb);
    val constellations = constellationData.constellationData.flatMap { constellation =>
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
    putTweet(startTime + day + (12.0 + 35.0 / 60) / 24, msg, starNames);
  }

  // この時期21時ごろ見える南の空低い星座
  def putTweetConstellationsSouth(day: Int): Boolean = {
    val decThres = - PI5 + tokyoLat + 30.0 / PI57;
    val time = startTime + day + 21.0 / 24;
    val tdb = TimeLib.mjdutcToTdb(time);
    val bpnMatrix = bpn.icrsToTrueEquatorialMatrix(tdb);
    val constellations = constellationData.constellationData.flatMap { constellation =>
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
      putTweet(startTime + day + (12.0 + 35.0 / 60) / 24, msg, starNames);
      true;
    } else {
      false;
    }
  }

  // この時期21時ごろ見える明るい星
  def putTweetBrightStars21(day: Int): Unit = {
    val altThres = 10.0 / PI57;
    //val altThres = 0.0;
    val time = startTime + day + 21.0 / 24;
    val tdb = TimeLib.mjdutcToTdb(time);
    val bpnMatrix = bpn.icrsToTrueEquatorialMatrix(tdb);
    val constellations = (constellationData.starData.flatMap { case (ra, xyz, name, hashtags) =>
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
      val (xyz, azi, alt) = calcPlanetXyzAziAlt(time, targetPlanet);
      if (alt >= altThres) {
        Some((azi, name));
      } else {
        None;
      }
    }).toIndexedSeq;
    val constellationsStr = constellations.sortBy(-_._1).map(_._2).mkString("、");
    val starNames = constellations.sortBy(-_._1).map(_._2).toList;
    val msg = "この時期21時ごろ見える明るい星は、%sです #星空".format(constellationsStr);
    putTweet(startTime + day + (12.0 + 35.0 / 60) / 24, msg, starNames);
  }

  // この時期日没時の一番星となりうる明るい星
  def putTweetFirstStar(day: Int): Unit = {
    val altHor = -0.90 / PI57;
    val altThres30 = 30.0 / PI57;
    //val altThres = 0.0;
    val time = sunsetTimes(day);
    val tdb = TimeLib.mjdutcToTdb(time);
    val bpnMatrix = bpn.icrsToTrueEquatorialMatrix(tdb);
    val constellations0 = ((constellationData.starData2.flatMap { case (ra, xyz, magnitude, name, hashtags) =>
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
      val (xyz, azi, alt) = calcPlanetXyzAziAlt(time, targetPlanet);
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
        magnitude + (altThres30 - alt) * PI57 / 8.0;
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
      putTweet(SunsetFirstStarTweetContent(day, msg, starNames));
    }
  }

  // この時期21時ごろ見えやすい黄道十二星座
  def putTweetConstellationsEcliptical(day: Int): Boolean = {
    val altThres = 30.0 / PI57;
    val time = startTime + day + 21.0 / 24;
    val tdb = TimeLib.mjdutcToTdb(time);
    val bpnMatrix = bpn.icrsToTrueEquatorialMatrix(tdb);
    val constellations = constellationData.constellationData.filter(_.eclipticalFlag).flatMap { constellation =>
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
      putTweet(startTime + day + (12.0 + 35.0 / 60) / 24, msg, starNames);
      true;
    } else {
      false;
    }
  }

  // この時期21時ごろ見える天の川
  def putTweetConstellationsGalaxy(day: Int): Boolean = {
    val altThres = 30.0 / PI57;
    val time = startTime + day + 21.0 / 24;
    val tdb = TimeLib.mjdutcToTdb(time);
    val bpnMatrix = bpn.icrsToTrueEquatorialMatrix(tdb);
    val constellations = constellationData.constellationData.filter(_.galaxyFlag).flatMap { constellation =>
      val xyz2 = VectorLib.multiplyMV(bpnMatrix, constellation.xyz);
      val (azi, alt) = hcs.trueEquatorialXyzToAziAltFromUtc(xyz2, time);
      if (alt >= altThres) {
        val sid = hcs.siderealTimeFromUtc(time);
        val ra2 = constellation.ra - sid;
        val ra3 = if (ra2 < -PI) ra2 + PI2 else if (ra2 >= PI) ra2 - PI2 else ra2;
        Some((ra3, constellation.name));
      } else {
        None;
      }
    }
    if (constellations.nonEmpty) {
      val constellationsStr = constellations.sortBy(_._1).map(_._2).mkString("、");
      val starNames = constellations.sortBy(-_._1).map(_._2).toList;
      val msg = "この時期21時ごろ見える天の川は、%sを通っています #星空 #星座".format(constellationsStr);
      putTweet(startTime + day + (12.0 + 35.0 / 60) / 24, msg, starNames);
      true;
    } else {
      false;
    }
  }

  // 恒星の南中
  def putTweetCulminations(): Unit = {
    val altHor = -0.90 / PI57;
    var index: Int = -1;
    val culminationContents = constellationData.culminationContents;
    (87 until period).foreach { day => // PERIOD
      if (index < 0) {
        val time = startTime + day + 20.9 / 24.0; // PERIOD
        val sid = hcs.siderealTimeFromUtc(time);
        index = culminationContents.indexWhere(_._1 > sid);
        if (index < 0) {
          index = 0;
        }
      } else {
        val timeS = startTime + day + 21.0 / 24.0;
        val timeE = startTime + day + 23.0 / 24.0;
        val time1 = hcs.siderealTimeToUtc(culminationContents(index)._1, timeS);
        val time0 = if (time1 < timeS) {
          time1;
        } else if (getTweets(timeS).tweets.map(_.time).filter(isNightTime3).isEmpty) {
          if (time1 < timeE) {
            time1;
          } else {
            0.0;
          }
        } else {
          0.0;
        }
        if (time0 > 0.0) {
          val msg = if (calcPlanetXyzAziAlt(time0, JplData.Moon)._3 >= altHor) {
            culminationContents(index)._2;
          } else {
            culminationContents(index)._2 + "。月明かりなし";
          }
          val urlOpt = culminationContents(index)._3;
          val hashtags = culminationContents(index)._4;
          val starNames = culminationContents(index)._5;
          putTweet(StarTweetContent(time0, msg, urlOpt, hashtags, starNames));
          index += 1;
          if (index == culminationContents.size) {
            index = 0;
          }
        }
      }
    }
  }

  // 星座の解説など
  def putTweetLunchTimeContents(): Unit = {
    val lunchTimeContents = constellationData.lunchTimeContents;
    var day1: Int = 62; // PERIOD
    var day2: Int = 62; // PERIOD
    var index: Int = {
      val time = startTime + day1 + 21.0 / 24.0; // PERIOD
      val sid = hcs.siderealTimeFromUtc(time);
      val index = lunchTimeContents.indexWhere(_._1 > sid);
      if (index < 0) {
        0;
      } else {
        index;
      }
    }
    while (day1 < period) {
      val sid = hcs.siderealTimeFromUtc(startTime + day1 + 21.0 / 24);
      if (MathLib.circleAdd(sid, -lunchTimeContents(index)._1) >= 0) {
        val msg = lunchTimeContents(index)._2;
        {
          val p = (day1 until (day2 + 14)).indexWhere { d =>
            !getTweets(startTime + d).tweets.exists(tc => DateTweets.isDayTime(tc.time));
          }
          day1 = if (p < 0) day1 else day1 + p;
        }
        val urlOpt = lunchTimeContents(index)._3;
        val hashtags = lunchTimeContents(index)._4;
        val starNames = Nil;
        putTweet(StarTweetContent(startTime + day1 + 12.0 / 24 + 5.0 / 60 / 24, msg, urlOpt, hashtags, starNames));
        index += 1;
        if (index == lunchTimeContents.size) {
          index = 0;
        }
      }
      day1 += 1;
      day2 += 1;
    }
  }

  putTweetCulminations();

  {
    var kind: Int = 0;
    (90 until period).foreach { day => // PERIOD
      val wday = TimeLib.wday(startTime + day);
      if (wday == 1) { // 月曜
        if ((startTime + day).toInt % 2 == 0) {
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
        if ((startTime + day).toInt % 2 == 1) {
          putTweetBrightStars21(day);
        }
      }
    }
  }
  (99 until period).foreach { day => // PERIOD
    val wday = TimeLib.wday(startTime + day);
    if (wday == 1) {
      if ((startTime + day).toInt % 2 == 0) {
        if (getTweets(startTime + day).sunsetTweets.isEmpty) {
          //putTweetFirstStar(day);
        }
      }
    }
  }

  putTweetLunchTimeContents();

  (words.periodStart until period).foreach { day =>
    if (getTweets(startTime + day).size < 3 || !getTweets(startTime + day).tweets.map(_.time).exists(isLunchTime)) {
      words.putTweetWord(day);
    }
  }
}

//==============================================================================
// なにもツイートのない日付

{
  (0 until period).foreach { day =>
    val time = startTime + day;
    if (getTweets(time).isEmpty) {
      putTweet(time, "#empty");
    } else {
      //if (getTweets(time).tweets.map(_.time).filter(isNightTime0).isEmpty) {
      if (getTweets(time).tweets.map(_.time).filter(isNightTime3).isEmpty) {
        putTweet(time + 23.0 / 24, "#night empty");
      } else if (!getTweets(time).tweets.exists(tc => DateTweets.isDayTime(tc.time))) {
        putTweet(time + 9.0 / 24, "#daytime empty");
      }
    }
  }
}

//==============================================================================
// ツイート出力

words.saveHistory(wordHistoryPath);

scala.util.Using(new java.io.PrintWriter(new java.io.FileOutputStream(dataPath))) { writer =>
  (0 until period).foreach { day =>
    getTweets(startTime + day).tweets.foreach { case tc =>
      val time = tc.time;
      if (time >= startTime1 && time < endTime1) {
        val msg = tc.tweetContent;
        writer.println("%s %s".format(TimeLib.timeToDateTimeString(time), msg));
      }
    }
  }
}

//==============================================================================

  def main(args: Array[String]): Unit = {
  }
}

