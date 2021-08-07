
object SolarTerms {

  def tweets: Seq[TweetContent] = {
    val termStrs = IndexedSeq(
      "春分", "清明", "穀雨", "立夏", "小満", "芒種", "夏至", "小暑", "大暑", "立秋", "処暑", "白露",
      "秋分", "寒露", "霜降", "立冬", "小雪", "大雪", "冬至", "小寒", "大寒", "立春", "雨水", "啓蟄",
    );
    val sunPhaseTerms = MathLib.findCyclicPhaseListContinuous(24, Period.startTime, Period.endTime, 3.0, 24) { time =>
      JplData.calcPlanetLngTrueEcliptic(time, JplData.Sun);
    }
    sunPhaseTerms.map { case (time, term) =>
      val tweetTime = TimeLib.floor(time, 24) + 1.0 / (24 * 4);
      val msg = if (term % 6 == 0) {
        "%s。太陽の黄経が%d°です".format(termStrs(term), term * 15);
      } else {
        "二十四節気の%s。太陽の黄経が%d°です".format(termStrs(term), term * 15);
      }
      LegacyTweetContent(tweetTime, msg, None, Nil);
    }
  }

}

