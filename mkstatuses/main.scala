import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.time.DayOfWeek;
import java.time.Instant;
import java.time.OffsetDateTime;
import java.time.ZoneOffset;

val planetsDataPath = "../var/jpl.dat";

val startTime = OffsetDateTime.parse(args(0) + "T00:00:00+09:00").toInstant;
val endTime   = OffsetDateTime.parse(args(1) + "T00:00:00+09:00").toInstant;

val jst = ZoneOffset.ofHours(+9);
val pi57 = 180.0 / Math.PI;
val pi = Math.PI;
val pi2 = 2.0 * Math.PI;
val pi5 = 0.5 * Math.PI;

println("# time: %s - %s".format(startTime, endTime));

val PLANET_COUNT = 7;
val ELEMENT_COUNT = 6;

val SUN_OFFSET     = 0 * ELEMENT_COUNT;
val MOON_OFFSET    = 1 * ELEMENT_COUNT;
val MERCURY_OFFSET = 2 * ELEMENT_COUNT;
val VENUS_OFFSET   = 3 * ELEMENT_COUNT;
val MARS_OFFSET    = 4 * ELEMENT_COUNT;
val JUPITER_OFFSET = 5 * ELEMENT_COUNT;
val SATURN_OFFSET  = 6 * ELEMENT_COUNT;
val DISTANCE_IDX       = 0;
val J2000_LNG_IDX      = 1;
val J2000_LAT_IDX      = 2;
val T_ECLIPTIC_LNG_IDX = 3;
val HCS_AZI_IDX        = 4;
val HCS_ALT_IDX        = 5;

def readJplData(): IndexedSeq[(Instant, OffsetDateTime, IndexedSeq[Double])] = {
  val input = new DataInputStream(new BufferedInputStream(new FileInputStream(planetsDataPath)));
  val durationCount = (endTime.getEpochSecond() - startTime.getEpochSecond()).toInt / 600;
  val data = (0 until durationCount).map { i =>
    val time = startTime.plusSeconds(i * 600);
    val data = (0 until (PLANET_COUNT * ELEMENT_COUNT)).map { _ =>
      input.readDouble();
    }
    val timeJst = OffsetDateTime.ofInstant(time, jst);
    (time, timeJst, data);
  }
  input.close();
  data;
}

def putMessage(timeJst: OffsetDateTime, message: String): Unit = {
  val timeStr = timeJst.toString.substring(0, 16);
  println("%s %s".format(timeStr, message));
}

println("# reading ...");

val jplData = readJplData();
val durationCount = jplData.size;

val durationCount6 = durationCount / 6;
val jplData6 = (0 until durationCount6).map(i => jplData(i * 6));

val durationCount24 = durationCount6 / 24;

val sunriseSunsetTimes: IndexedSeq[(Double, Double)] = {
  val altHor = -0.90 / pi57;
  (0 until durationCount24).map { day =>
    val sunrise = {
      val s = day * 144 + (4 * 6);
      val alts = (s to s+24).map { i =>
        jplData(i)._3(SUN_OFFSET + HCS_ALT_IDX);
      }
      val p = alts.indexWhere(_ > altHor);
      val alt1 = alts(p-1);
      val alt2 = alts(p);
      s.toDouble + p - 1.0 - (alt1 - altHor) / (alt2 - alt1);
    }
    val sunset = {
      val s = day * 144 + (16 * 6);
      val alts = (s to s+24).map { i =>
        jplData(i)._3(SUN_OFFSET + HCS_ALT_IDX);
      }
      val p = alts.indexWhere(_ < altHor);
      val alt1 = alts(p-1);
      val alt2 = alts(p);
      s.toDouble + p - 1.0 + (alt1 - altHor) / (alt1 - alt2);
    }
    (sunrise, sunset);
  }
}

def jplDataTime(time: Int): OffsetDateTime = {
  jplData(time)._2;
}

def jplDataTime(time: Double): OffsetDateTime = {
  val time1 = time.toInt;
  val time2 = time - time1;
  jplData(time1)._2.plusSeconds((time2 * 600.0).toInt);
}

def jplDataPlanet(time: Int, planetOffset: Int): IndexedSeq[Double] = {
  (0 until ELEMENT_COUNT).map { i =>
    jplData(time)._3(planetOffset + i);
  }
}

def jplDataPlanet(time: Double, planetOffset: Int): IndexedSeq[Double] = {
  val time1 = time.toInt;
  val time2 = time - time1;
  if (time2 == 0.0) {
    jplDataPlanet(time1, planetOffset);
  } else {
    (0 until ELEMENT_COUNT).map { i =>
      val v1 = jplData(time1)._3(planetOffset + i);
      val v2 = jplData(time1 + 1)._3(planetOffset + i);
      val v2b = if (i == DISTANCE_IDX) {
        v2;
      } else if (v2 - v1 > pi) {
        v2 - pi2;
      } else if (v1 - v2 > pi) {
        v2 + pi2;
      } else {
        v2;
      }
      val v3 = v1 + (v2b - v1) * time2;
      val v3b = if (i == J2000_LNG_IDX || i == T_ECLIPTIC_LNG_IDX || i == HCS_AZI_IDX) {
        if (v3 >= pi2) {
          v3 - pi2;
        } else if (v3 < 0) {
          v3 + pi2;
        } else {
          v3;
        }
      } else if (i == J2000_LAT_IDX || i == HCS_ALT_IDX) {
        if (v3 >= pi) {
          v3 - pi2;
        } else if (v3 < -pi) {
          v3 + pi2;
        } else {
          v3;
        }
      } else {
        v3;
      }
      v3b;
    }
  }
}

def calcLngDiff(lng1: Double, lng2: Double): Double = {
  val d = lng1 - lng2;
  if (d < 0) d + pi2 else d;
}

val constellationsMap = Map (
  " 0h00m, -5" -> "うお座の西側の魚(ペガスス座の南)の南",

  " 0h20m, -5" -> "くじら座のしっぽ付近",

  " 1h00m,  0" -> "うお座の西側の魚(ペガスス座の南)のしっぽ付近",

  " 1h20m, 10" -> "うお座の東側の魚(アンドロメダ座の南)のしっぽ付近",

  " 1h40m,  5" -> "うお座の2匹の魚の付け根付近",
  " 1h40m, 10" -> "うお座の東側の魚(アンドロメダ座の南)のしっぽ付近",

  " 2h00m,  5" -> "くじら座の頭付近",
  " 2h00m, 10" -> "おひつじ座の頭とくじら座の頭の間",
  " 2h00m, 15" -> "おひつじ座の頭付近",

  " 2h20m, 10" -> "おひつじ座とくじら座の頭の間",
  " 2h20m, 15" -> "おひつじ座の南",

  " 2h40m, 10" -> "おひつじ座とくじら座の頭の間",
  " 2h40m, 15" -> "おひつじ座の南",

  " 3h00m, 15" -> "おひつじ座のしっぽ側",

  " 3h20m, 15" -> "おうし座の西",
  " 3h20m, 20" -> "おうし座すばる付近",

  " 3h40m, 15" -> "おうし座すばるの南",
  " 3h40m, 20" -> "おうし座すばる付近",

  " 4h00m, 15" -> "おうし座ヒアデスの西",
  " 4h00m, 20" -> "おうし座すばる付近",

  " 4h20m, 20" -> "おうし座とペルセウス座とぎょしゃ座の間",

  " 4h40m, 20" -> "おうし座ヒアデスの北東",

  " 5h00m, 20" -> "おうし座の角付近",

  " 5h20m, 20" -> "おうし座の角とぎょしゃ座付近",

  " 5h40m, 20" -> "おうし座の角とオリオン座の腕付近",

  " 6h00m, 20" -> "ふたご座の西側の子の足元とオリオン座の腕付近",

  " 6h20m, 20" -> "ふたご座の西側の子の足元付近",
  " 6h20m, 25" -> "ふたご座の西側の子の腰付近",

  " 6h40m, 20" -> "ふたご座の2人の足元付近",
  " 6h40m, 25" -> "ふたご座の西側の子の胴体付近",

  " 7h00m, 20" -> "ふたご座の東側の子の腰付近",
  " 7h00m, 25" -> "ふたご座の西側の子の胴体付近",

  " 7h20m, 20" -> "ふたご座の東側の子の胴体付近",
  " 7h20m, 25" -> "ふたご座の東側の子の胴体付近",

  " 7h40m, 20" -> "ふたご座ポルックスの南でふたご座の東",

  " 8h00m, 20" -> "かに座の西",

  " 8h20m, 20" -> "かに座",

  " 8h40m, 20" -> "かに座",

  " 9h00m, 20" -> "かに座の東",

  " 9h20m, 15" -> "しし座とかに座の間",

  " 9h40m, 15" -> "しし座の西",

  "10h00m, 15" -> "しし座の肩付近",

  "10h20m, 10" -> "しし座レグルスの東",
  "10h20m, 15" -> "しし座の中央",

  "10h40m, 10" -> "しし座の腹付近",

  "11h00m, 10" -> "しし座の後ろ足付近",

  "11h20m,  5" -> "しし座の後ろ足とおとめ座の間",

  "11h40m,  5" -> "おとめ座の頭付近",

  "12h00m,  0" -> "おとめ座の四角形の西",
  "12h00m,  5" -> "おとめ座の頭付近",

  "12h20m,  0" -> "おとめ座の四角形の西",

  "12h40m,  0" -> "おとめ座の四角形",
  "12h40m, -5" -> "おとめ座の四角形",

  "13h00m, -5" -> "おとめ座の四角形",

  "13h20m, -5" -> "おとめ座スピカの北",
  "13h20m,-10" -> "おとめ座スピカの北",

  "13h40m,-10" -> "おとめ座スピカの東",

  "14h00m,-10" -> "おとめ座スピカの東",
  "14h00m,-15" -> "おとめ座スピカとてんびん座の間",

  "14h20m,-15" -> "てんびん座の西おとめ座との間",

  "14h40m,-15" -> "てんびん座の西",

  "15h00m,-20" -> "てんびん座",

  "15h20m,-20" -> "てんびん座",

  "15h40m,-20" -> "てんびん座の東さそり座との間",

  "16h00m,-25" -> "さそり座アンタレスの北西",

  "16h20m,-25" -> "へびつかい座でさそり座アンタレスの北",

  "16h40m,-25" -> "へびつかい座でさそり座アンタレスの北東",

  "17h00m,-25" -> "へびつかい座でさそり座アンタレスの東",

  "17h20m,-25" -> "へびつかい座でさそり座アンタレスの東",

  "17h40m,-30" -> "いて座さそり座へびつかい座の間",

  "18h00m,-30" -> "いて座の西",

  "18h20m,-30" -> "いて座南斗六星",

  "18h40m,-30" -> "いて座南斗六星",

  "19h00m,-30" -> "いて座南斗六星の先",

  "19h20m,-30" -> "いて座南斗六星の東",

  "20h00m,-25" -> "やぎ座の西部",

  "20h20m,-20" -> "やぎ座の西部",

  "20h40m,-20" -> "やぎ座の中央",
  "20h40m,-25" -> "やぎ座の中央",

  "21h00m,-20" -> "やぎ座の中央",
  "21h00m,-25" -> "やぎ座の中央",

  "21h40m,-15" -> "やぎ座とみずがめ座の間",
  "21h40m,-20" -> "やぎ座の東部",

  "22h00m,-20" -> "みずがめ座の南部でやぎ座の東",
  "22h00m,-15" -> "みずがめ座の南部でやぎ座の東",

  "22h20m,-15" -> "みずがめ座の南部でやぎ座の東",

  "22h40m,-15" -> "みずがめ座の南部でうお座の西側の魚の頭の南西",

  "23h00m,-15" -> "みずがめ座の南部でうお座の西側の魚の頭の南",

  "23h20m,-10" -> "みずがめ座とうお座の西側の魚の頭の間",

  "23h40m,-10" -> "みずがめ座とうお座の西側の魚の頭の間",
);

def j2000ToConstellations(lng: Double, lat: Double): (String, String) = {
  val lng5 = (lng * pi57 / 5).toInt * 5;
  val lat5 = ((lat * pi57 + 90) / 5).toInt * 5 - 90;
  val key = "%2dh%02dm,%3d".format(lng5 / 15, lng5 % 15 * 4, lat5);
  val cons = constellationsMap.getOrElse(key, "-");
  if (cons == "") {
    ("#", "(%s)".format(key));
  } else if (cons == "-") {
    ("##", "(%s)".format(key));
  } else {
    ("", cons);
  }
}

def procSun1(): Unit = {
  val sunLng360 = (0 until durationCount6).map { i =>
    val d = jplData6(i)._3(SUN_OFFSET + T_ECLIPTIC_LNG_IDX);
    d * pi57;
  }

  val sunTermStr = IndexedSeq(
    "春分", "清明", "穀雨", "立夏", "小満", "芒種", "夏至", "小暑", "大暑", "立秋", "処暑", "白露",
    "秋分", "寒露", "霜降", "立冬", "小雪", "大雪", "冬至", "小寒", "大寒", "立春", "雨水", "啓蟄",
  );
  (1 until durationCount6).map { i =>
    val v1 = (sunLng360(i-1) / 15).toInt;
    val v2 = (sunLng360(i) / 15).toInt;
    if (v1 != v2) {
      if (v2 % 6 != 0) {
        (i, "二十四節気の%s。太陽の黄経が%d°です".format(sunTermStr(v2), v2 * 15));
      } else {
        (i, "%s。太陽の黄経が%d°です".format(sunTermStr(v2), v2 * 15));
      }
    } else {
      (i, "");
    }
  }.filter(_._2 != "").foreach { case (i, msg) =>
    putMessage(jplData6(i)._2.minusSeconds(2700), msg);
  }
}

def procSun2(): Unit = {
  (1 until durationCount6 - 1).map { i =>
    val d0 = jplData6(i - 1)._3(SUN_OFFSET + DISTANCE_IDX);
    val d1 = jplData6(i    )._3(SUN_OFFSET + DISTANCE_IDX);
    val d2 = jplData6(i + 1)._3(SUN_OFFSET + DISTANCE_IDX);
    if (d1 < d0 && d1 <= d2) {
      (i, "地球が近日点通過");
    } else if (d1 > d0 && d1 >= d2) {
      (i, "地球が遠日点通過");
    } else {
      (i, "");
    }
  }.filter(_._2 != "").foreach { case (i, msg) =>
    putMessage(jplData6(i)._2.minusSeconds(900), msg);
  }
}

def procSun3(): Unit = {
  (1 until durationCount24 - 1).flatMap { day =>
    val dayOfWeek = jplDataTime(day * 144).getDayOfWeek();
    val d0 = sunriseSunsetTimes(day - 1)._2 - (day - 1) * 144;
    val d1 = sunriseSunsetTimes(day    )._2 - (day    ) * 144;
    val d2 = sunriseSunsetTimes(day + 1)._2 - (day + 1) * 144;
    if (d1 < d0 && d1 <= d2) {
      val time = jplDataTime(d1 + day * 144);
      val hour = time.getHour();
      val minute = time.getMinute();
      Some((time, "日没が最も早く、%d時%02d分ごろです".format(hour, minute)));
    } else if (d1 > d0 && d1 >= d2) {
      val time = jplDataTime(d1 + day * 144);
      val hour = time.getHour();
      val minute = time.getMinute();
      Some((time, "日没が最も遅く、%d時%02d分ごろです".format(hour, minute)));
    } else {
      None;
    }
  }.foreach { case (time, msg) =>
    putMessage(time.minusSeconds(600), msg);
  }
}

def procSun4(): Unit = {
  (0 until durationCount24).flatMap { day =>
    val dayOfWeek = jplDataTime(day * 144).getDayOfWeek();
    if (dayOfWeek == DayOfWeek.SUNDAY) {
      val time = jplDataTime(sunriseSunsetTimes(day)._2);
      val hour = time.getHour();
      val minute = time.getMinute();
      val msg = "日没は%d時%02d分ごろ".format(hour, minute);
      Some((time, msg));
    } else {
      None;
    }
  }.foreach { case (time, msg) =>
    putMessage(time.minusSeconds(600), msg);
  }
}

def procMoon(): Unit = {
  val moonSunLng360 = (0 until durationCount6).map { i =>
    val d = jplDataPlanet(i * 6, MOON_OFFSET)(T_ECLIPTIC_LNG_IDX) - jplDataPlanet(i * 6, SUN_OFFSET)(T_ECLIPTIC_LNG_IDX);
    val d360 = d * pi57;
    if (d360 < 0) d360 + 360 else d360;
  }
  val moonTermStr = IndexedSeq("新月", "上弦の月", "満月", "下弦の月");

  val fullMoonDistanceData: (Map[Int, Int], IndexedSeq[Double]) = {
    val moonSunLng = (0 until durationCount).map { i =>
      val d1 = jplDataPlanet(i, MOON_OFFSET)(T_ECLIPTIC_LNG_IDX);
      val d0 = jplDataPlanet(i, SUN_OFFSET)(T_ECLIPTIC_LNG_IDX);
      calcLngDiff(d1, d0);
    }
    val timeDistance = (1 until durationCount).flatMap { i =>
      if (moonSunLng(i) >= pi && moonSunLng(i - 1) < pi) {
        val time = i.toDouble - 1.0 + (pi - moonSunLng(i - 1)) / (moonSunLng(i) - moonSunLng(i - 1));
        val distance = jplDataPlanet(time, MOON_OFFSET)(DISTANCE_IDX);
        Some((time, distance));
      } else {
        None;
      }
    }
    val data1 = timeDistance.zipWithIndex.map { case ((time, _), idx) =>
      ((time / 6).toInt, idx);
    }.toMap;
    val data2 = timeDistance.map(_._2);
    (data1, data2);
  }

  (1 until durationCount6).map { i =>
    val v1 = (moonSunLng360(i-1) / 90).toInt;
    val v2 = (moonSunLng360(i) / 90).toInt;
    if (v1 != v2) {
      if (v2 == 2) {
        val idx = fullMoonDistanceData._1.getOrElse(i - 1, 0);
        if (idx >= 1 && idx < fullMoonDistanceData._2.size - 1) {
          val d0 = fullMoonDistanceData._2(idx - 1);
          val d1 = fullMoonDistanceData._2(idx    );
          val d2 = fullMoonDistanceData._2(idx + 1);
          if (d1 < d0 && d1 <= d2) {
            (i * 6, "満月。月が地球に近いため、もっとも大きい満月です");
          } else if (d1 > d0 && d1 >= d2) {
            (i * 6, "満月。月が地球に遠いため、もっとも小さい満月です");
          } else {
            (i * 6, moonTermStr(v2));
          }
        } else {
          (i * 6, moonTermStr(v2));
        }
      } else {
        (i * 6, moonTermStr(v2));
      }
    } else {
      (i * 6, "");
    }
  }.filter(_._2 != "").foreach { case (i, msg) =>
    putMessage(jplDataTime(i).minusSeconds(2700), msg);
  }
}

def procPlanets1(): Unit = {
  val planets = IndexedSeq(
    (MERCURY_OFFSET, "水星"),
    (VENUS_OFFSET,   "金星"),
  );
  planets.foreach { case (offset, planetName) =>
    val data = (0 until durationCount6).map { i =>
      val d1 = jplDataPlanet(i * 6, SUN_OFFSET);
      val d2 = jplDataPlanet(i * 6, offset);
      calcLngDiff(d2(T_ECLIPTIC_LNG_IDX), d1(T_ECLIPTIC_LNG_IDX));
      //calcLngDiff(d2(J2000_LNG_IDX), d1(J2000_LNG_IDX));
    }
    (1 until durationCount6 - 1).foreach { i =>
      val d0 = data(i - 1);
      val d1 = data(i    );
      val d2 = data(i + 1);
      if (d0 > pi && d1 < pi) {
        val msg = "%sが外合(黄経基準)".format(planetName);
        putMessage(jplDataTime(i * 6).minusSeconds(2700), msg);
      } else if (d0 < pi && d1 > pi) {
        val msg = "%sが内合(黄経基準)".format(planetName);
        putMessage(jplDataTime(i * 6).minusSeconds(2700), msg);
      } else if (d0 < pi && d1 < pi && d2 < pi && d1 > d0 && d1 >= d2) {
        val msg = "%sが東方最大離角(黄経基準)".format(planetName);
        putMessage(jplDataTime(i * 6).minusSeconds(900), msg);
      } else if (d0 > pi && d1 > pi && d2 > pi && d1 < d0 && d1 <= d2) {
        val msg = "%sが西方最大離角(黄経基準)".format(planetName);
        putMessage(jplDataTime(i * 6).minusSeconds(900), msg);
      }
    }
  }
}

def procPlanets2(): Unit = {
  val planets = IndexedSeq(
    (MARS_OFFSET,    "火星"),
    (JUPITER_OFFSET, "木星"),
    (SATURN_OFFSET,  "土星"),
  );
  val termStr = IndexedSeq("合", "東矩", "衝", "西矩");
  planets.foreach { case (offset, planetName) =>
    val data = (0 until durationCount6).map { i =>
      val d1 = jplDataPlanet(i * 6, SUN_OFFSET);
      val d2 = jplDataPlanet(i * 6, offset);
      //val d = calcLngDiff(d2(T_ECLIPTIC_LNG_IDX), d1(T_ECLIPTIC_LNG_IDX));
      val d = calcLngDiff(d2(J2000_LNG_IDX), d1(J2000_LNG_IDX));
      ((d * pi57 / 90).toInt, d2(J2000_LNG_IDX), d2(J2000_LAT_IDX));
    }
    (1 until durationCount6).foreach { i =>
      if (data(i)._1 != data(i - 1)._1) {
        val v = data(i-1)._1;
        val msg = if (v == 0) {
          "%sが%s(赤経基準)".format(planetName, termStr(v));
        } else {
          val (flag, cons) = j2000ToConstellations(data(i)._2, data(i)._3);
          "%s%sが%s(赤経基準)。%sにいます".format(flag, planetName, termStr(v), cons);
        }
        putMessage(jplDataTime(i * 6).minusSeconds(2700), msg);
      }
    }
  }
}

def procPlanets3(): Unit = {
  val planets = IndexedSeq(
    (MERCURY_OFFSET, "水星"),
    (VENUS_OFFSET,   "金星"),
  );
  val data = planets.map { case (planetOffset, _) =>
    (0 until durationCount24).map { day =>
      val time = sunriseSunsetTimes(day)._2;
      val d = jplDataPlanet(time, planetOffset);
      d(HCS_ALT_IDX);
    }
  }
  (planets zip data).foreach { case ((_, planetName), data) =>
    (1 until durationCount24 - 1).flatMap { day =>
      val d0 = data(day - 1);
      val d1 = data(day    );
      val d2 = data(day + 1);
      if (d1 > d0 && d1 >= d2) {
        val alt = (d1 * pi57 + 0.5).toInt;
        Some((sunriseSunsetTimes(day)._2, "%sは日没時最大高度で西の空高度約%d°".format(planetName, alt)));
      } else {
        None;
      }
    }.foreach { case (time, msg) =>
      val time2 = (time * 2).toInt * 0.5;
      putMessage(jplDataTime(time2), msg);
    }
  }
}

def procPlanets4(): Unit = {
  val planets = IndexedSeq(
    (MERCURY_OFFSET, "水星", List(DayOfWeek.WEDNESDAY)),
    (VENUS_OFFSET,   "金星", List(DayOfWeek.FRIDAY)),
  );
  val altThres = 10 / pi57;
  (1 until durationCount24).foreach { day =>
    val dayOfWeek = jplDataTime(day * 144).getDayOfWeek();
    val time = sunriseSunsetTimes(day)._2;
    val timePrev = sunriseSunsetTimes(day - 1)._2;
    planets.filter(_._3.contains(dayOfWeek)).
      foreach { case (offset, planetName, _) =>
      val d0 = jplDataPlanet(timePrev, offset)(HCS_ALT_IDX);
      val d1 = jplDataPlanet(time, offset)(HCS_ALT_IDX);
      if (d1 >= altThres) {
        val alt = (d1 * pi57 + 0.5).toInt;
        val msg = if (d1 > d0) {
          "%sは日没時の高度を徐々に上げ、西の空高度約%d°にいます".format(planetName, alt);
        } else {
          "%sは日没時の高度を徐々に下げ、西の空高度約%d°にいます".format(planetName, alt);
        }
        val time2 = (time * 2).toInt * 0.5;
        putMessage(jplDataTime(time2), msg);
      }
    }
  }
}

def procPlanets5(): Unit = {
  val planets = IndexedSeq(
    (MOON_OFFSET, "月", List(DayOfWeek.MONDAY, DayOfWeek.TUESDAY, DayOfWeek.WEDNESDAY,
      DayOfWeek.THURSDAY, DayOfWeek.FRIDAY, DayOfWeek.SATURDAY, DayOfWeek.SUNDAY)),
    (MARS_OFFSET,    "火星", List(DayOfWeek.TUESDAY)),
    (JUPITER_OFFSET, "木星", List(DayOfWeek.THURSDAY)),
    (SATURN_OFFSET,  "土星", List(DayOfWeek.SATURDAY)),
  );
  val altThres = 10 / pi57;
  (0 until durationCount24 - 1).foreach { day =>
    val dayOfWeek = jplDataTime(day * 144).getDayOfWeek();
    val timeList = IndexedSeq(day * 144 + 21 * 6, day * 144 + 23 * 6);
    planets.filter(_._3.contains(dayOfWeek)).
      foreach { case (offset, planetName, _) =>
      val d = timeList.map { time =>
        val d = jplDataPlanet(time, offset);
        (d(J2000_LNG_IDX), d(J2000_LAT_IDX), d(HCS_AZI_IDX), d(HCS_ALT_IDX));
      }
      val timeDelta = if (offset == MOON_OFFSET) 600 else 900;
      if (d(0)._4 >= altThres) {
        val (flag, cons) = j2000ToConstellations(d(0)._1, d(0)._2);
        val msg = "%s%sは%sにいます".format(flag, planetName, cons);
        putMessage(jplDataTime(day * 144 + 21 * 6).minusSeconds(timeDelta), msg);
      } else if (d(1)._4 >= altThres && offset != MOON_OFFSET) {
        val (flag, cons) = j2000ToConstellations(d(1)._1, d(1)._2);
        val msg = "%s%sは%sにいます".format(flag, planetName, cons);
        putMessage(jplDataTime(day * 144 + 23 * 6).minusSeconds(timeDelta), msg);
      }
    }
  }
}

procSun1();
procSun2();
procSun3();
procSun4();
procMoon();

procPlanets1();
procPlanets2();
procPlanets3();
procPlanets4();
procPlanets5();

