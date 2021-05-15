
trait TweetContent {
  def time: Double;
  def message: String;
  def urlOpt: Option[String] = None;
  def hashtags: List[String];
  def starNames: List[String];

  def date: String = TimeLib.timeToDateString(time);
  def tweetContent: String = {
    val (msg, length) = urlOpt match {
      case None      =>
        val msg = message + hashtags.map(" #" + _).mkString;
        val length = msg.length;
        (msg, length);
      case Some(url) =>
        val msg = message + " " + url + hashtags.map(" #" + _).mkString;
        val length = msg.length - url.length - 1 + 12;
        (msg, length);
    }
    if (msg.indexOf("ERROR") >= 0) {
      "#" + msg;
    } else if (length > 140) {
      "#TOO_LONG(%d) %s".format(length, msg);
    } else {
      msg;
    }
  }
}

case class LegacyTweetContent(time: Double, message: String,
  override val urlOpt: Option[String], starNames: List[String]) extends  TweetContent {
  def hashtags: List[String] = Nil;
}

case class DateTweets(otherTweets: List[TweetContent], daytimeTweets: List[TweetContent],
  sunsetTweets: List[Main.OnSunsetTweetContent],
  nightTweets: List[TweetContent]) {
  def isEmpty: Boolean = otherTweets.isEmpty && sunsetTweets.isEmpty && daytimeTweets.isEmpty && nightTweets.isEmpty;
  def size: Int = otherTweets.size + daytimeTweets.size + nightTweets.size + (if (sunsetTweets.isEmpty) 0 else 1);
  def added(tc: TweetContent): DateTweets = {
    tc match {
      case tc: Main.OnSunsetTweetContent =>
        this.copy(sunsetTweets = tc :: this.sunsetTweets);
      case _ =>
        if (DateTweets.isDayTime(tc.time)) {
          this.copy(daytimeTweets = tc :: this.daytimeTweets);
        } else if (Main.isNightTime0(tc.time)) {
          this.copy(nightTweets = tc :: this.nightTweets);
        } else {
          this.copy(otherTweets = tc :: this.otherTweets);
        }
    }
  }
  def removed(time: Double, message: String): DateTweets = {
    def p(tc: TweetContent): Boolean = {
      TimeLib.timeToDateTimeString(tc.time) == TimeLib.timeToDateTimeString(time) &&
      tc.message == message;
    }
    if ((otherTweets ::: daytimeTweets ::: sunsetTweets ::: nightTweets).exists(p)) {
      DateTweets(otherTweets.filter(!p(_)), daytimeTweets.filter(!p(_)),
        sunsetTweets.filter(!p(_)), nightTweets.filter(!p(_)));
    } else {
      added(LegacyTweetContent(time, "#-" + message, None, Nil));
    }
  }
  def tweets: List[TweetContent] = {
    import Ordering.Double.IeeeOrdering;
    ((if (sunsetTweets.isEmpty) {
      Nil;
    } else if (sunsetTweets.tail.isEmpty) {
      sunsetTweets.head :: Nil;
    } else {
      Main.MultiSunsetTweetContent(sunsetTweets.head.day, sunsetTweets.reverse) :: Nil;
    }) ::: otherTweets ::: daytimeTweets ::: nightTweets).sortBy(_.time);
  }
}

object DateTweets {

  def isDayTime(time: Double): Boolean = {
    val day = (time - Main.startTime).toInt;
    val s = Main.sunsetTimes(day);
    time - time.toInt >= 0.0 / 24 && time < s; // 9時～日没
  }

}

