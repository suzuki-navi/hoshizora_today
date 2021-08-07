
object Constellations {

  val PI57 = 180.0 / Math.PI;
  val PI2 = Math.PI * 2.0;

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

  val (culminationContents, lunchTimeContents, lunchTimeContents2, constellationData, starData, starData2, constellationsMap) = {
    import Ordering.Double.IeeeOrdering;

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

    val constellationsDataPath = "constellations.txt";
    val source = scala.io.Source.fromInputStream(
      getClass.getClassLoader.getResourceAsStream(constellationsDataPath));
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
        val (content0: String, urlOpt: Option[String], hashtags: List[String], starNames: List[String]) = parseContent(line.substring(7));
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

  val PI57 = 180.0 / Math.PI;

  trait LunchTimeContent {
    def word: String;
    def content(day: Int, hcs: Hcs): String;
    def urlOpt: Option[String];
    def hashtags: List[String];
  }

  case class ConstellationLunchTimeContent (word: String, contentTemplate: String,
    urlOpt: Option[String], hashtags: List[String],
    ra: Double, dec: Double, xyz: Array[Double]) extends LunchTimeContent {
    def content(day: Int, hcs: Hcs): String = {
      val altHor = -0.90 / PI57;
      val time = Period.startTime + day + 21.0 / 24;
      val tdb = TimeLib.mjdutcToTdb(time);
      val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
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

class Words() {

  private def fetchLunchTimeContents(word: String): IndexedSeq[Words.LunchTimeContent] = {
    Constellations.lunchTimeContents2.filter(_.word == word);
  }

  private val defaultLimit = 25;
  private var history: Map[String, (List[Int], Int)] = Map.empty;
  private var contents: Map[String, IndexedSeq[Words.LunchTimeContent]] = Map.empty;

  val periodStart: Int = 109; // PERIOD 2021/07/18

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
          val line = "%s:%s".format(word, TimeLib.timeToDateString(Period.startTime + day));
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

  private def findRecentWord(day: Int, hcs: Hcs, tweets: Tweets): Option[(String, IndexedSeq[Words.LunchTimeContent])] = {
    val time = Period.startTime + day + 12.0 / 24;
    val lastTweets = (-13 to 0).flatMap(i => tweets.getTweets(Period.startTime + day + i).tweets);
    val lst = lastTweets.flatMap(_.starNames);
    if (lst.nonEmpty) {
      import Ordering.Double.IeeeOrdering;
      val lst2 = lst.map { w => (w, tweetPoint(w, day)) }.sortBy(- _._2._2);
      lst2.find { case (word, (contents, point)) =>
        point >= 1.0 && contents.filter(_.content(day, hcs) != "").nonEmpty;
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

  def putTweetWord(day: Int, hcs: Hcs, tweets: Tweets): Unit = {
    findRecentWord(day, hcs, tweets) match {
      case None => ;
      case Some((word, contents)) =>
        val time1 = Period.startTime + day + (12.0 + 50.0 / 60) / 24;
        var inc: Int = 0;
        contents.foreach { content =>
          val c = content.content(day + inc / 4, hcs);
          if (c != "") {
            val msg = c + content.hashtags.map(" #" + _).mkString("");
            tweets.putTweet(time1 + 1.0 * (inc % 4) / 24, msg, content.urlOpt);
            if (inc == 0) {
              addHistory(word, day + inc / 4);
            }
            inc += 1;
          }
        }
    }
  }

}

