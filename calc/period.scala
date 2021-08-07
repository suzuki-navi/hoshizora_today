
object Period {

  // PERIOD

  val startTime = TimeLib.stringToModifiedJulianDay ("2021-01-31T00:00:00+09:00");
  val startTime1 = TimeLib.stringToModifiedJulianDay("2021-07-18T00:00:00+09:00");
  val endTime = TimeLib.stringToModifiedJulianDay   ("2022-07-18T00:00:00+09:00");

  val period = (endTime - startTime).toInt;

}
