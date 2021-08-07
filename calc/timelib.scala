import java.time.OffsetDateTime;

object TimeLib {

  // 文字列から修正ユリウス日に変換
  def stringToModifiedJulianDay(timeStr: String): Double = {
    val offsetDateTime = OffsetDateTime.parse(timeStr)
    val second = offsetDateTime.toEpochSecond; // 1970-01-01T00:00:00Zからの秒数
    second.toDouble / 86400.0 + 40587.0;
  }

  def modifiedJulianDayToString(time: Double): String = {
    java.time.Instant.ofEpochSecond(((time - 40587.0) * 86400.0).toLong).toString;
  }

  def dateStringToDay(str: String): Int = {
    (stringToModifiedJulianDay(str + "T00:00:00Z") - Period.startTime).toInt;
  }

  // 2021-04-01
  def timeToDateString(time: Double): String = {
    java.time.Instant.ofEpochSecond(((time - 40587.0) * 86400.0 + 0.5).toLong + 9 * 3600).toString.substring(0, 10);
  }

  // 2021-04-01T12:00
  def timeToDateTimeString(time: Double): String = {
    java.time.Instant.ofEpochSecond(((time - 40587.0) * 86400.0 + 0.5).toLong + 9 * 3600).toString.substring(0, 16);
  }

  // 12時00分
  def timeToTimeNaturalString(time: Double): String = {
    val str = java.time.Instant.ofEpochSecond(((time - 40587.0) * 86400.0 + 0.5).toLong + 9 * 3600).toString;
    "%s時%s分".format(str.substring(11, 13), str.substring(14, 16));
  }

  def wday(time: Double): Int = {
    ((time + (9.0 / 24.0 + 0.5 / 86400.0)).toInt + 3) % 7;
  }

  private def offsetDateTime(time: Double): OffsetDateTime = {
    java.time.Instant.ofEpochSecond(((time - 40587.0) * 86400.0 + 0.5).toLong).atOffset(java.time.ZoneOffset.ofHours(9));
  }

  def dayOfMonth(time: Double): Int = {
    offsetDateTime(time).getDayOfMonth;
  }

  def dayOfYear(time: Double): Int = {
    offsetDateTime(time).getDayOfYear;
  }

  def month(time: Double): Int = {
    offsetDateTime(time).getMonthValue;
  }

  def floor(time: Double, stepCountPerDay: Int): Double = {
    val t1 = time.toInt;
    val t2 = time - t1;
    if (stepCountPerDay == 1) {
      if (t2 < 15.0 / 24) {
        t1.toDouble - 9.0 / 24;
      } else {
        t1.toDouble + 15.0 / 24;
      }
    } else {
      t1 + (t2 * stepCountPerDay).toInt.toDouble / stepCountPerDay;
    }
  }

  def round(time: Double, stepCountPerDay: Int): Double = {
    val t1 = time.toInt;
    val t2 = time - t1;
    if (stepCountPerDay == 1) {
      if (t2 < 15.0 / 24) {
        t1.toDouble - 9.0 / 24;
      } else {
        t1.toDouble + 15.0 / 24;
      }
    } else {
      t1 + (t2 * stepCountPerDay + 0.5).toInt.toDouble / stepCountPerDay;
    }
  }

  // "来年1月"
  def monthString(time0: Double, day1: Int, day2: Int): String = {
    def sub(time: Double): (Int, Int) = {
      val str = java.time.Instant.ofEpochSecond(((time - 40587.0) * 86400.0 + 0.5).toLong + 9 * 3600).toString;
      val year = str.substring(0, 4);
      val month = str.substring(5, 7);
      (year.toInt, month.toInt);
    }
    val (y0, m0) = sub(time0);
    val (y1, m1) = sub(Period.startTime + day1);
    val (y2, m2) = sub(Period.startTime + day2);
    val s1 = if (y0 == y1 && m0 == m1) {
      "今月";
    } else if (y0 == y1 && m0 + 1 == m1 || y0 + 1 == y1 && m0 - 11 == m1) {
      "来月"
    } else {
      (if (y0 == y1) {
        "";
      } else if (y0 + 1 == y1) {
        "来年";
      } else {
        "%d年".format(y1);
      }) + "%d月".format(m1);
    }
    if (y1 == y2 && m1 == m2) {
      s1;
    } else {
      s1 + "から" + (if (y0 == y2 || y1 == y2) {
        "";
      } else if (y0 + 1 == y2) {
        "来年";
      } else {
        "%d年".format(y2);
      }) + "%d月".format(m2) + "にかけて";
    }
  }

  // UTCからUT1に変換
  def mjdutcToUt1(time: Double): Double = {
    time; // 近似的に同一とする
  }

  // UTCからTDBに変換
  def mjdutcToTdb(time: Double): Double = {
    // 近似的
    //if (time >= 57754.0) {
      time + (37.0 + 32.184) / 86400.0;
    //}
  }

  // UT1から平均恒星時に変換
  def mjdut1ToGmst(time: Double): Double = {
    val r = 0.280191 + gmstK * (time - 59215.0);
    Const.PI2 * (if (r < 0) {
      r + 1 - r.toInt;
    } else {
      r - r.toInt;
    });
  }

  def gmstToMjdut1(sid: Double, time0: Double): Double = {
    val d = MathLib.circleAdd(sid, -mjdut1ToGmst(time0));
    time0 + d / Const.PI2 / gmstK;
  }

  private val gmstK = 1.0027379094;

}

