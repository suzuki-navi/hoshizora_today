import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.time.DayOfWeek;
import java.time.Instant;
import java.time.OffsetDateTime;
import java.time.ZoneOffset;

val planetsDataPath = "../var/planets.dat";

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
  ( 25, 10) -> "うお座の東側の魚(アンドロメダ座の南)のしっぽ付近",
  ( 30, 10) -> "おひつじ座の頭とくじら座の頭の間",
  ( 35, 10) -> "おひつじ座とくじら座の頭の間",
  ( 40, 10) -> "おひつじ座とくじら座の頭の間",
  ( 45, 15) -> "おひつじ座のしっぽ側",
  ( 60, 15) -> "おうし座ヒアデスの西",
  ( 65, 20) -> "おうし座とペルセウス座とぎょしゃ座の間",
  ( 70, 20) -> "おうし座ヒアデスの北東",
  ( 75, 20) -> "おうし座の角付近",
  ( 80, 20) -> "おうし座の角とぎょしゃ座付近",
  ( 85, 20) -> "おうし座の角とオリオン座の腕付近",
  ( 90, 20) -> "ふたご座の西側の子の足元とオリオン座の腕付近",
  ( 95, 20) -> "ふたご座の西側の子の足元付近",
  (100, 20) -> "ふたご座の2人の足元付近",
  (100, 25) -> "ふたご座の西側の子の胴体付近",
  (105, 20) -> "ふたご座の東側の子の腰付近",
  (110, 20) -> "ふたご座の東側の子の胴体付近",
  (115, 20) -> "ふたご座ポルックスの南でふたご座の東",
  (120, 20) -> "かに座の西",
  (125, 20) -> "かに座",
  (130, 20) -> "かに座",
  (135, 20) -> "かに座の東",
  (140, 15) -> "しし座とかに座の間",
  (145, 15) -> "しし座の西",
  (150, 15) -> "しし座の肩付近",
  (155, 15) -> "しし座の中央",
  (160, 10) -> "しし座の腹付近",
  (165, 10) -> "しし座の後ろ足付近",
  (175,  5) -> "おとめ座の頭付近",
  (190,  0) -> "おとめ座の胴体付近",
  (200,-10) -> "おとめ座スピカの北",
);

def j2000ToConstellations(lng: Double, lat: Double): (String, String) = {
  val lng5 = (lng * pi57 / 5).toInt * 5;
  val lat5 = ((lat * pi57 + 90) / 5).toInt * 5 - 90;
  val key = (lng5, lat5);
  val cons = constellationsMap.getOrElse(key, "");
  if (cons == "") {
    ("#", "(%dh%02dm(%d), %d)".format(lng5 / 15, lng5 % 15 * 4, lng5, lat5));
  } else {
    ("", cons);
  }
}

def procSun(): Unit = {
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

def procMoon(): Unit = {
  val moonSunLng360 = (0 until durationCount6).map { i =>
    val d = jplData6(i)._3(MOON_OFFSET + T_ECLIPTIC_LNG_IDX) - jplData6(i)._3(SUN_OFFSET + T_ECLIPTIC_LNG_IDX);
    val d360 = d * pi57;
    if (d360 < 0) d360 + 360 else d360;
  }

  val moonTermStr = IndexedSeq("新月", "上弦の月", "満月", "下弦の月");
  (1 until durationCount6).map { i =>
    val v1 = (moonSunLng360(i-1) / 90).toInt;
    val v2 = (moonSunLng360(i) / 90).toInt;
    if (v1 != v2) {
      (i, moonTermStr(v2));
    } else {
      (i, "");
    }
  }.filter(_._2 != "").foreach { case (i, msg) =>
    putMessage(jplData6(i)._2.minusSeconds(2700), msg);
  }
}

procSun();
procMoon();

def procPlanets1(): Unit = {
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
        val (flag, cons) = j2000ToConstellations(data(i)._2, data(i)._3);
        val msg = "%s%sが%s。%sにいます".format(flag, planetName, termStr(data(i-1)._1), cons);
        putMessage(jplDataTime(i * 6).minusSeconds(2700), msg);
      }
    }
  }
}

def procPlanets2(): Unit = {
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
    val timeList = IndexedSeq(sunriseSunsetTimes(day)._2, day * 144.0 + 21 * 6, day * 144.0 + 23 * 6);
    planets.filter(_._3.contains(dayOfWeek)).
      foreach { case (offset, planetName, _) =>
      val d = timeList.map { time =>
        val d = jplDataPlanet(time, offset);
        (d(J2000_LNG_IDX), d(J2000_LAT_IDX), d(HCS_AZI_IDX), d(HCS_ALT_IDX));
      }
      (d, planetName);
      if (d(1)._4 >= altThres) {
        val (flag, cons) = j2000ToConstellations(d(1)._1, d(1)._2);
        val msg = "%s%sは%sにいます".format(flag, planetName, cons);
        putMessage(jplDataTime(day * 144 + 21 * 6).minusSeconds(600), msg);
      } else if (d(2)._4 >= altThres && offset != MOON_OFFSET) {
        val (flag, cons) = j2000ToConstellations(d(1)._1, d(1)._2);
        val msg = "%s%sは%sにいます".format(flag, planetName, cons);
        putMessage(jplDataTime(day * 144 + 23 * 6).minusSeconds(600), msg);
      }
    }
  }
}

procPlanets1();
procPlanets2();

