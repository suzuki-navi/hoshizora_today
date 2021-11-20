
object Diff {

  private val diffDataPath = "diff.txt";

  private lazy val diffData = load();

  private def load(): List[(Double, Int, String)] = {
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

  def removeTweets1(tweetsManager: Tweets): Unit = {
    diffData.filter(_._2 < 0).foreach { d =>
      tweetsManager.removeTweet(d._1, d._3);
    }
  }

  def putTweets1(tweetsManager: Tweets): Unit = {
    diffData.filter(_._2 > 0).foreach { d =>
      tweetsManager.putTweet(d._1, d._3);
    }
  }

}
