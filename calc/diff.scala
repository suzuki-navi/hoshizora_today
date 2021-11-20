
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
        val s = cols(1);
        if (s.startsWith("++")) {
          diffData = (time, +2, s.substring(2)) :: diffData;
        } else if (s.startsWith("--")) {
          diffData = (time, -2, s.substring(2)) :: diffData;
        } else if (s.startsWith("+")) {
          diffData = (time, +1, s.substring(1)) :: diffData;
        } else if (s.startsWith("-")) {
          diffData = (time, -1, s.substring(1)) :: diffData;
        }
      }
    }
    source.close();
    import Ordering.Double.IeeeOrdering;
    diffData.reverse.sortBy(_._1);
  }

  def removeTweets1(tweetsManager: Tweets): Unit = {
    diffData.filter(_._2 == -1).foreach { d =>
      tweetsManager.removeTweet(d._1, d._3, true);
    }
  }

  def putTweets1(tweetsManager: Tweets): Unit = {
    diffData.filter(_._2 == +1).foreach { d =>
      tweetsManager.putTweet(d._1, d._3);
    }
  }

  def removeTweets2(tweetsManager: Tweets): Unit = {
    diffData.filter(_._2 == -2).foreach { d =>
      tweetsManager.removeTweet(d._1, d._3, true);
    }
  }

  def putTweets2(tweetsManager: Tweets): Unit = {
    diffData.filter(_._2 == +2).foreach { d =>
      tweetsManager.putTweet(d._1, d._3);
    }
  }

}
