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

println("# time: %s - %s".format(startTime, endTime));

val PLANET_COUNT = 3;
val ELEMENT_COUNT = 6;

val SUN_OFFSET  = 0 * ELEMENT_COUNT;
val MOON_OFFSET = 1 * ELEMENT_COUNT;
val MARS_OFFSET = 2 * ELEMENT_COUNT;
val DISTANCE_IDX       = 0;
val J2000_LNG_IDX      = 1;
val J2000_LAT_IDX      = 2;
val T_ECLIPTIC_LNG_IDX = 3;
val HCS_AZI_IDX        = 4;
val HCS_ALT_IDX        = 5;

def readData(): IndexedSeq[(Instant, OffsetDateTime, IndexedSeq[Double])] = {
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

val data = readData();
val durationCount = data.size;

val durationCount6 = durationCount / 6;
val data6 = (0 until durationCount6).map(i => data(i * 6));

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

def j2000ToConstellations(lng: Double, lat: Double): String = {
  val lng5 = (lng * pi57 / 5).toInt * 5;
  val lat5 = ((lat * pi57 + 90) / 5).toInt * 5 - 90;
  val key = (lng5, lat5);
  val cons = constellationsMap.getOrElse(key, "");
  if (cons == "") {
    println("# ERROR %d(%dh+%dm), %d".format(lng5, lng5 / 15, lng5 % 15 * 4, lat5));
  }
  cons;
}

def calcSun(): Unit = {
  val sunLng360 = (0 until durationCount6).map { i =>
    val d = data6(i)._3(SUN_OFFSET + T_ECLIPTIC_LNG_IDX);
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
    putMessage(data6(i)._2.minusSeconds(2700), msg);
  }

  (1 until durationCount6 - 1).map { i =>
    val d0 = data6(i - 1)._3(SUN_OFFSET + DISTANCE_IDX);
    val d1 = data6(i    )._3(SUN_OFFSET + DISTANCE_IDX);
    val d2 = data6(i + 1)._3(SUN_OFFSET + DISTANCE_IDX);
    if (d1 < d0 && d1 <= d2) {
      (i, "地球が近日点通過");
    } else if (d1 > d0 && d1 >= d2) {
      (i, "地球が遠日点通過");
    } else {
      (i, "");
    }
  }.filter(_._2 != "").foreach { case (i, msg) =>
    putMessage(data6(i)._2.minusSeconds(900), msg);
  }
}

def calcMoon(): Unit = {
  val moonSunLng360 = (0 until durationCount6).map { i =>
    val d = data6(i)._3(MOON_OFFSET + T_ECLIPTIC_LNG_IDX) - data6(i)._3(SUN_OFFSET + T_ECLIPTIC_LNG_IDX);
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
    putMessage(data6(i)._2.minusSeconds(2700), msg);
  }

  (0 until durationCount6).map { i =>
    if (data6(i)._2.getHour() == 21) {
      if (data6(i)._3(MOON_OFFSET + HCS_ALT_IDX) * pi57 >= 20) {
      //val d = moonSunLng360(i);
      //if (d >= 90 && d < 225) {
        val lng = data6(i)._3(MOON_OFFSET + J2000_LNG_IDX);
        val lat = data6(i)._3(MOON_OFFSET + J2000_LAT_IDX);
        val cons = j2000ToConstellations(lng, lat);
        if (cons != "") {
          (i, "月は%sにいます".format(cons));
        } else {
          (i, "");
        }
      } else {
        (i, "");
      }
    } else {
      (i, "");
    }
  }.filter(_._2 != "").foreach { case (i, msg) =>
    putMessage(data6(i)._2.minusSeconds(600), msg);
  }
}

def calcMars(): Unit = {
  val marsSunLng360 = (0 until durationCount6).map { i =>
    val d = data6(i)._3(MARS_OFFSET + T_ECLIPTIC_LNG_IDX) - data6(i)._3(SUN_OFFSET + T_ECLIPTIC_LNG_IDX);
    val d360 = d * pi57;
    if (d360 < 0) d360 + 360 else d360;
  }

  val marsTermStr = IndexedSeq("合", "東矩", "衝", "西矩");
  (1 until durationCount6).map { i =>
    val v1 = (marsSunLng360(i-1) / 90).toInt;
    val v2 = (marsSunLng360(i) / 90).toInt;
    if (v1 != v2) {
      (i, "火星が%s".format(marsTermStr(v2)));
    } else {
      (i, "");
    }
  }.filter(_._2 != "").foreach { case (i, msg) =>
    putMessage(data6(i)._2.minusSeconds(2700), msg);
  }

  (0 until durationCount6).map { i =>
    val localTime = data6(i)._2;
    if (localTime.getHour() == 21 && localTime.getDayOfWeek() == DayOfWeek.TUESDAY) {
      val d = marsSunLng360(i);
      if (d >= 90 && d < 225) {
        val lng = data6(i)._3(MARS_OFFSET + J2000_LNG_IDX);
        val lat = data6(i)._3(MARS_OFFSET + J2000_LAT_IDX);
        val cons = j2000ToConstellations(lng, lat);
        if (cons != "") {
          (i, "火星は%sにいます".format(cons));
        } else {
          (i, "");
        }
      } else {
        (i, "");
      }
    } else {
      (i, "");
    }
  }.filter(_._2 != "").foreach { case (i, msg) =>
    putMessage(data6(i)._2.minusSeconds(600), msg);
  }
}

calcSun();
calcMoon();
calcMars();

