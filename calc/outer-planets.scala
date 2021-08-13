
object OuterPlanets {

  def nightTweets(hcs: Hcs): Seq[TweetContent] = {
    (0 until Period.period).flatMap { day =>
      val time = Period.startTime + day;
      val wday = TimeLib.wday(time);
      (if (wday == 2) {
        Some(("火星", JplData.Mars));
      } else if (wday == 4) {
        Some(("木星", JplData.Jupiter));
      } else if (wday == 6) {
        Some(("土星", JplData.Saturn));
      } else {
        None;
      }).flatMap { case (planetName, targetPlanet) =>
        nightTweet(targetPlanet, planetName, hcs, day);
      }
    }
  }

  private def nightTweet(targetPlanet: JplData.TargetPlanet, planetName: String,
    hcs: Hcs, day: Int): Option[TweetContent] = {
    val altThres = 15 / Const.PI57;
    {
      val time = Period.startTime + day + 21.0 / 24.0;
      val (azi, alt) = Acs.Horizontal.calcPlanetAziAlt(time, targetPlanet, hcs);
      if (alt >= altThres) {
        Some((time, azi, alt));
      } else {
        val time = Period.startTime + day + 23.0 / 24.0;
        val (azi, alt) = Acs.Horizontal.calcPlanetAziAlt(time, targetPlanet, hcs);
        if (alt >= altThres) {
          Some((time, azi, alt));
        } else {
          None;
        }
      }
    }.map { case (time, azi, alt) =>
      val xyz = Acs.Icrs.calcPlanetXyz(time, targetPlanet);
      val (cons, hashtags) = Constellations.icrsToConstellation(xyz);
      val hcsStr = Hcs.aziAltToNaturalString(azi, alt);
      val msg = "%sは%s、%sにいます".format(planetName, hcsStr, cons);
      if (day <= 136) { // PERIOD 2021/08/14 互換性対策
        LegacyTweetContent(time - 1.0 / (24 * 4), msg + " #%s".format(planetName) +
          hashtags.map(" #" + _).mkString, None, Nil);
      } else {
        BasicTweetContent(time - 1.0 / (24 * 4), msg, None, planetName :: hashtags, planetName :: Nil);
      }
    }
  }

}

