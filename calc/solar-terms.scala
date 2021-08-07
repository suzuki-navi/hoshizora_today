
object SolarTerms {

  def tweets: Seq[TweetContent] = {
    val termStrs: IndexedSeq[(String, String)] = IndexedSeq(
      ("春分", ""),
      ("清明", ""),
      ("穀雨", ""),
      ("立夏", "夏の気配が感じられるころです"),
      ("小満", ""),
      ("芒種", ""),
      ("夏至", "日中の太陽がもっとも空高い時期です"),
      ("小暑", ""),
      ("大暑", "暑さがもっとも厳しいころです"),
      ("立秋", "秋の気配が感じられるはずのころです"),
      ("処暑", "暑さがおさまるころです"),
      ("白露", ""),
      ("秋分", ""),
      ("寒露", ""),
      ("霜降", ""),
      ("立冬", "冬の気配が感じられるころです"),
      ("小雪", ""),
      ("大雪", ""),
      ("冬至", "日中の太陽がもっとも空低い時期です"),
      ("小寒", ""),
      ("大寒", "寒さがもっとも厳しいころです"),
      ("立春", "春の気配が感じられるころです"),
      ("雨水", ""),
      ("啓蟄", ""),
    );
    val sunPhaseTerms = MathLib.findCyclicPhaseListContinuous(24, Period.startTime, Period.endTime, 3.0, 24) { time =>
      JplData.calcPlanetLngTrueEcliptic(time, JplData.Sun);
    }
    sunPhaseTerms.map { case (time, term) =>
      val tweetTime = TimeLib.floor(time, 24) + 1.0 / (24 * 4);
      val msg1 = if (term % 6 == 0) {
        "";
      } else {
        "二十四節気の";
      }
      val msg2 = if (term % 6 == 0) {
        "。二十四節気のひとつで、";
      } else {
        "。";
      }
      val msg3 = if (termStrs(term)._2 == "") {
        "";
      } else {
        "。" + termStrs(term)._2;
      }
      val msg = msg1 + termStrs(term)._1 + msg2 + "太陽の黄経が%d°です".format(term * 15) + msg3;
      LegacyTweetContent(tweetTime, msg, None, Nil);
    }
  }

}

