
// 星座のツイート

object TweetConstellations {

  // この時期21時ごろ見えやすい星座
  def constellations21Tweet(day: Int, hcs: Hcs): TweetContent = {
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
    LegacyTweetContent(Period.startTime + day + (12.0 + 35.0 / 60) / 24, msg, None, starNames);
  }

}

