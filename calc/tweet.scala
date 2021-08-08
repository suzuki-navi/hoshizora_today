
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

case class DateTweets(otherTweets: List[TweetContent],
  sunsetTweets: List[Main.OnSunsetTweetContent]) {
  def isEmpty: Boolean = otherTweets.isEmpty && sunsetTweets.isEmpty;
  def size: Int = otherTweets.size + (if (sunsetTweets.isEmpty) 0 else 1);
  def added(tc: TweetContent): DateTweets = {
    tc match {
      case tc: Main.OnSunsetTweetContent =>
        this.copy(sunsetTweets = tc :: this.sunsetTweets);
      case _ =>
        this.copy(otherTweets = tc :: this.otherTweets);
    }
  }
  def removed(time: Double, message: String): DateTweets = {
    def p(tc: TweetContent): Boolean = {
      TimeLib.timeToDateTimeString(tc.time) == TimeLib.timeToDateTimeString(time) &&
      tc.message == message || tc.tweetContent == message;
    }
    if (sunsetTweet.nonEmpty && p(sunsetTweet.get)) {
      this.copy(sunsetTweets = Nil);
    } else if (otherTweets.exists(p)) {
      this.copy(otherTweets = otherTweets.filter(!p(_)));
    } else {
      added(LegacyTweetContent(time, "#-" + message, None, Nil));
    }
  }
  val sunsetTweet: Option[TweetContent] = {
    if (sunsetTweets.isEmpty) {
      None;
    } else if (sunsetTweets.tail.isEmpty) {
      Some(sunsetTweets.head);
    } else {
      Some(Main.MultiSunsetTweetContent(sunsetTweets.head.day, sunsetTweets.reverse));
    }
  }
  def tweets: List[TweetContent] = {
    import Ordering.Double.IeeeOrdering;
    (sunsetTweet match {
      case Some(sunsetTweet) => sunsetTweet :: otherTweets;
      case None => otherTweets;
    }).sortBy(_.time);
  }
}

object DateTweets {

  def apply(): DateTweets = DateTweets(Nil, Nil);

  def isDayTime(time: Double): Boolean = {
    val day = (time - Period.startTime).toInt;
    val s = Main.sunsetTimes(day);
    time - time.toInt >= 0.0 / 24 && time < s; // 9時～日没
  }

}

class Tweets() {

  var _tweets: Map[String, DateTweets] = (0 until Period.period).map { day =>
    val date = TimeLib.timeToDateString(Period.startTime + day);
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

}

