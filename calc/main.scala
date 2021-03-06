import java.time.OffsetDateTime;

val startTime = TimeLib.stringToModifiedJulianDay(args(0) + "T00:00:00+09:00");
val endTime = TimeLib.stringToModifiedJulianDay(args(1) + "T00:00:00+09:00");
val period = (endTime - startTime).toInt;


val jplDataPath = "../var/ssd.jpl.nasa.gov/pub/eph/planets/ascii/de430/ascp1950.430";
val nutLsDataPath = "nut-ls.txt";
val nutPlDataPath = "nut-pl.txt";

val PI = Math.PI;
val PI2 = Math.PI * 2.0;
val PI2R = 1.0 / PI2;
val PI5 = Math.PI * 0.5;
val PI57 = 180.0 / Math.PI;
val PI_AS2R = Math.PI / (3600 * 180);

val AU = 1.49597870700000000e+8;

val tokyoLng = 139.7 / PI57;
val tokyoLat = 35.7 / PI57;

val DEBUG_FLAG = true;

System.err.println("Started calculating...");

object JplData {

  val EMRAT = 81.3005690741906200; // 地球と月の質量比
  val EMRAT1 = 1.0 / (1.0 + EMRAT);

  sealed trait JplPlanet {
    def dataOffset: Int;
    def subPeriodCount: Int;
    def coefficientCount: Int;
  }
  object MERCURY_JPL extends JplPlanet {
    val dataOffset = 2;
    val subPeriodCount = 4;
    val coefficientCount = 14;
  }
  object VENUS_JPL extends JplPlanet {
    val dataOffset = 170;
    val subPeriodCount = 2;
    val coefficientCount = 10;
  }
  object EARTH_MOON_BARYCENTER_JPL extends JplPlanet {
    val dataOffset = 230;
    val subPeriodCount = 2;
    val coefficientCount = 13;
  }
  object MARS_JPL extends JplPlanet {
    val dataOffset = 308;
    val subPeriodCount = 1;
    val coefficientCount = 11;
  }
  object JUPITER_JPL extends JplPlanet {
    val dataOffset = 341;
    val subPeriodCount = 1;
    val coefficientCount = 8;
  }
  object SATURN_JPL extends JplPlanet {
    val dataOffset = 365;
    val subPeriodCount = 1;
    val coefficientCount = 7;
  }
  object MOON_JPL extends JplPlanet {
    val dataOffset = 440;
    val subPeriodCount = 8;
    val coefficientCount = 13;
  }
  object SUN_JPL extends JplPlanet {
    val dataOffset = 752;
    val subPeriodCount = 2;
    val coefficientCount = 11;
  }

  sealed trait TargetPlanet;
  object Mercury extends TargetPlanet;
  object Venus   extends TargetPlanet;
  object Mars    extends TargetPlanet;
  object Jupiter extends TargetPlanet;
  object Saturn  extends TargetPlanet;

}

class JplData(dataPath: String) {

  private[this] val jplData: IndexedSeq[IndexedSeq[Double]] = loadJplData(dataPath);

  private[this] def loadJplData(dataPath: String): IndexedSeq[IndexedSeq[Double]] = {
    val data1 = {
      val sc = new java.util.Scanner(new java.io.BufferedInputStream(new java.io.FileInputStream(dataPath)));
      val data1 = scala.collection.mutable.ArrayBuffer[String]();
      while (sc.hasNext) {
        data1 += sc.next();
      }
      sc.close();
      data1;
    }
    val periodCount = data1.size / 1022;
    val data2: IndexedSeq[IndexedSeq[Double]] = (0 until periodCount).map { i =>
      (0 until 1018).map { j =>
        data1(i * 1022 + 2 + j).replace("D", "E").toDouble;
      }
    }
    data2;
  }

  private[this] def getPeriodData(time: Double): (IndexedSeq[Double], Double) = {
    val time2 = (time - (jplData(0)(0) - 2400000.5)) / 32.0;
    val periodIndex = time2.toInt;
    (jplData(periodIndex), time2 - periodIndex);
  }

  private[this] def calcChebyshevArr(count: Int, x: Double): Array[Double] = {
    val arr = new Array[Double](count);
    val x2 = 2.0 * x;
    arr(0) = 1.0;
    arr(1) = x;
    var i: Int = 2;
    while (i < count) {
      arr(i) = x2 * arr(i-1) - arr(i-2);
      i = i + 1;
    }
    arr;
  }

  private[this] def calcChebyshev(periodData: IndexedSeq[Double], offset: Int, count: Int, arr: Array[Double]): Double = {
    var result: Double = 0.0;
    var i: Int = 0;
    while (i < count) {
      result += periodData(offset + i) * arr(i);
      i = i + 1;
    }
    result;
  }

  def calcPlanetPosition(time: Double, planet: JplData.JplPlanet): Array[Double] = {
    val (periodData, time2) = getPeriodData(time);
    val time3 = time2 * planet.subPeriodCount;
    val subPeriodIndex = time3.toInt;
    val time4 = 2.0 * (time3 - subPeriodIndex) - 1.0;
    val coefficientCount = planet.coefficientCount;
    val offset = 3 * coefficientCount * subPeriodIndex + planet.dataOffset;
    val charr = calcChebyshevArr(coefficientCount, time4);
    val ret = new Array[Double](3);
    ret(0) = calcChebyshev(periodData, offset                       , coefficientCount, charr);
    ret(1) = calcChebyshev(periodData, offset + 1 * coefficientCount, coefficientCount, charr);
    ret(2) = calcChebyshev(periodData, offset + 2 * coefficientCount, coefficientCount, charr);
    ret;
  }

  def calcEarthPosition(time: Double): Array[Double] = {
    val em = calcPlanetPosition(time, JplData.EARTH_MOON_BARYCENTER_JPL);
    val moon = calcPlanetPosition(time, JplData.MOON_JPL);
    val emrat1 = JplData.EMRAT1;
    val ret1 = new Array[Double](3);
    ret1(0) = em(0) - moon(0) * emrat1;
    ret1(1) = em(1) - moon(1) * emrat1;
    ret1(2) = em(2) - moon(2) * emrat1;
    ret1;
  }

  def calcEarthAndMoonPosition(time: Double): (Array[Double], Array[Double]) = {
    val em = calcPlanetPosition(time, JplData.EARTH_MOON_BARYCENTER_JPL);
    val moon = calcPlanetPosition(time, JplData.MOON_JPL);
    val emrat1 = JplData.EMRAT1;
    val ret1 = new Array[Double](3);
    ret1(0) = em(0) - moon(0) * emrat1;
    ret1(1) = em(1) - moon(1) * emrat1;
    ret1(2) = em(2) - moon(2) * emrat1;
    val ret2 = new Array[Double](3);
    ret2(0) = ret1(0) + moon(0);
    ret2(1) = ret1(1) + moon(1);
    ret2(2) = ret1(2) + moon(2);
    (ret1, ret2);
  }

  def calcSunFromEarth(time: Double): Array[Double] = {
    val earth = calcEarthPosition(time);
    val sun = calcPlanetPosition(time, JplData.SUN_JPL);
    VectorLib.minus(sun, earth);
  }

  def calcMoonFromEarth(time: Double): Array[Double] = {
    calcPlanetPosition(time, JplData.MOON_JPL);
  }

  def calcPlanetFromEarth(time: Double, targetPlanet: JplData.TargetPlanet): Array[Double] = {
    val earth = calcEarthPosition(time);
    val target = targetPlanet match {
      case JplData.Mercury => calcPlanetPosition(time, JplData.MERCURY_JPL);
      case JplData.Venus   => calcPlanetPosition(time, JplData.VENUS_JPL);
      case JplData.Mars    => calcPlanetPosition(time, JplData.MARS_JPL);
      case JplData.Jupiter => calcPlanetPosition(time, JplData.JUPITER_JPL);
      case JplData.Saturn  => calcPlanetPosition(time, JplData.SATURN_JPL);
    }
    VectorLib.minus(target, earth);
  }

}

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

  def modifiedJulianDayToStringJST(time: Double): String = {
    java.time.Instant.ofEpochSecond(((time - 40587.0) * 86400.0 + 0.5).toLong + 9 * 3600).toString.substring(0, 16);
  }

  def modifiedJulianDayToStringJSTNaturalTime(time: Double): String = {
    val str = java.time.Instant.ofEpochSecond(((time - 40587.0) * 86400.0 + 0.5).toLong + 9 * 3600).toString;
    "%s時%s分".format(str.substring(11, 13), str.substring(14, 16));
  }

  def modifiedJulianDayToStringJSTDate(time: Double): String = {
    val str = java.time.Instant.ofEpochSecond(((time - 40587.0) * 86400.0 + 0.5).toLong + 9 * 3600).toString;
    str.substring(0, 10);
  }

  def wday(time: Double): Int = {
    val wday1 = (time + (9.0 / 24.0 + 0.5 / 86400.0)).toInt % 7;
    if (wday1 >= 4) {
      wday1 - 4;
    } else {
      wday1 + 3;
    }
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
    val r = 0.280191 + 1.0027379094 * (time - 59215.0);
    PI2 * (if (r < 0) {
      r + 1 - r.toInt;
    } else {
      r - r.toInt;
    });
  }

}

object VectorLib {

  def minus(a: Array[Double], b: Array[Double]): Array[Double] = {
    val ret = new Array[Double](3);
    ret(0) = a(0) - b(0);
    ret(1) = a(1) - b(1);
    ret(2) = a(2) - b(2);
    ret;
  }

  def multiplyMV(m: Array[Double], v: Array[Double]): Array[Double] = {
    val ret = new Array[Double](3);
    ret(0) = m(0) * v(0) + m(1) * v(1) + m(2) * v(2);
    ret(1) = m(3) * v(0) + m(4) * v(1) + m(5) * v(2);
    ret(2) = m(6) * v(0) + m(7) * v(1) + m(8) * v(2);
    ret;
  }

  def distance(v: Array[Double]): Double = {
    val x = v(0);
    val y = v(1);
    val z = v(2);
    Math.sqrt(x * x + y * y + z * z);
  }

  val unitMatrix: Array[Double] = {
    val m = new Array[Double](9);
    m(0) = 1.0;
    m(4) = 1.0;
    m(8) = 1.0;
    m;
  }

  def rotateMatrixX(th: Double, r: Array[Double]): Array[Double] = {
    val cos = Math.cos(th);
    val sin = Math.sin(th);
    val ret = new Array[Double](9);
    ret(0) = r(0);
    ret(1) = r(1);
    ret(2) = r(2);
    ret(3) =  cos * r(3) - sin * r(6);
    ret(4) =  cos * r(4) - sin * r(7);
    ret(5) =  cos * r(5) - sin * r(8);
    ret(6) =  sin * r(3) + cos * r(6);
    ret(7) =  sin * r(4) + cos * r(7);
    ret(8) =  sin * r(5) + cos * r(8);
    ret;
  }

  def rotateMatrixY(th: Double, r: Array[Double]): Array[Double] = {
    val cos = Math.cos(th);
    val sin = Math.sin(th);
    val ret = new Array[Double](9);
    ret(0) =  cos * r(0) + sin * r(6);
    ret(1) =  cos * r(1) + sin * r(7);
    ret(2) =  cos * r(2) + sin * r(8);
    ret(3) =  r(3);
    ret(4) =  r(4);
    ret(5) =  r(5);
    ret(6) = -sin * r(0) + sin * r(6);
    ret(7) = -sin * r(1) + sin * r(7);
    ret(8) = -sin * r(2) + sin * r(8);
    ret;
  }

  def rotateMatrixZ(th: Double, r: Array[Double]): Array[Double] = {
    val cos = Math.cos(th);
    val sin = Math.sin(th);
    val ret = new Array[Double](9);
    ret(0) =  cos * r(0) - sin * r(3);
    ret(1) =  cos * r(1) - sin * r(4);
    ret(2) =  cos * r(2) - sin * r(5);
    ret(3) =  sin * r(0) + cos * r(3);
    ret(4) =  sin * r(1) + cos * r(4);
    ret(5) =  sin * r(2) + cos * r(5);
    ret(6) =  r(6);
    ret(7) =  r(7);
    ret(8) =  r(8);
    ret;
  }

  def xyzToLng(xyz: Array[Double]): Double = {
    val x = xyz(0);
    val y = xyz(1);
    val lng = Math.atan2(y, x);
    if (lng < 0) lng + PI2 else lng;
  }

  def xyzToLat(xyz: Array[Double]): Double = {
    val x = xyz(0);
    val y = xyz(1);
    val z = xyz(2);
    val r = Math.sqrt(x * x + y * y);
    val lat = Math.atan2(z, r);
    lat;
  }

  def calcLngDiff(lng1: Double, lng2: Double): Double = {
    val d = lng1 - lng2;
    if (d < 0) d + PI2 else d;
  }

}

object Bpn {

  private def calcGamma(time: Double): Double = {
    val jc = (time - 51544.5) / 36525.0;
    (-0.052928    +
    (10.556378    +
    ( 0.4932044   +
    (-0.00031238  +
    (-0.000002788 +
    ( 0.0000000260)
    * jc) * jc) * jc) * jc) * jc) * PI_AS2R;
  }

  private def calcPhi(time: Double): Double = {
    val jc = (time - 51544.5) / 36525.0;
    (84381.412819    +
    (  -46.811016    +
    (    0.0511268   +
    (    0.00053289  +
    (   -0.000000440 +
    (   -0.0000000176)
    * jc) * jc) * jc) * jc) * jc) * PI_AS2R;
  }

  private def calcPsi(time: Double): Double = {
    val jc = (time - 51544.5) / 36525.0;
    (  -0.041775    +
    (5038.481484    +
    (   1.5584175   +
    (  -0.00018522  +
    (  -0.000026452 +
    (  -0.0000000148)
    * jc) * jc) * jc) * jc) * jc) * PI_AS2R;
  }

  private def calcEps(time: Double): Double = {
    val jc = (time - 51544.5) / 36525.0;
    (84381.406        +
    (  -46.836769     +
    (   -0.0001831    +
    (    0.00200340   +
    (   -5.76 * 10e-7 +
    (   -4.34 * 10e-8)
    * jc) * jc) * jc) * jc) * jc) * PI_AS2R;
  }

  // ICRS座標系から平均黄道座標系に変換
  def icrsToMeanEclipticMatrix(time: Double): Array[Double] = {
    val gamma = calcGamma(time);
    val phi = calcPhi(time);
    val psi = calcPsi(time);
    var r: Array[Double] = VectorLib.unitMatrix;
    r = VectorLib.rotateMatrixZ(-gamma, r);
    r = VectorLib.rotateMatrixX(-phi, r);
    r = VectorLib.rotateMatrixZ(psi, r);
    r;
  }

  // ICRS座標系から平均赤道座標系に変換
  def icrsToMeanEquatorialMatrix(time: Double): Array[Double] = {
    val eps = calcEps(time);
    var r: Array[Double] = icrsToMeanEclipticMatrix(time);
    r = VectorLib.rotateMatrixX(eps, r);
    r;
  }

  // ICRS座標系から真黄道座標系に変換
  private def icrsToTrueEclipticMatrix(time: Double, nutation: (Double, Double)): Array[Double] = {
    val gamma = calcGamma(time);
    val phi = calcPhi(time);
    val psi = calcPsi(time);
    val (dpsi, deps) = nutation;
    var r: Array[Double] = VectorLib.unitMatrix;
    r = VectorLib.rotateMatrixZ(-gamma, r);
    r = VectorLib.rotateMatrixX(-phi, r);
    r = VectorLib.rotateMatrixZ(psi + dpsi, r);
    r;
  }

  // ICRS座標系から真黄道座標系に変換
  def icrsToTrueEclipticMatrix(time: Double): Array[Double] = {
    val nutation = calcNutation(time);
    icrsToTrueEclipticMatrix(time, nutation);
  }

  // ICRS座標系から真赤道座標系に変換
  def icrsToTrueEquatorialMatrix(time: Double): Array[Double] = {
    val nutation = calcNutation(time);
    val (dpsi, deps) = nutation;
    val eps = calcEps(time);
    var r: Array[Double] = icrsToTrueEclipticMatrix(time, nutation);
    r = VectorLib.rotateMatrixX(eps + deps, r);
    r;
  }

  // 章動
  private def calcNutation(time: Double): (Double, Double) = {
    val jc = (time - 51544.5) / 36525.0;
    val fj2 = -2.7774e-6 * jc;
    val nl = calcNutationLuniSolar(jc);
    val np = calcNutationPlanetary(jc);
    val dpsi = (nl._1 + np._1) * (1.0 + 0.4697e-6 * fj2);
    val deps = (nl._2 + np._2) * (1.0 + fj2);
    (dpsi, deps);
  }

  // 日月章動
  private def calcNutationLuniSolar(jc: Double): (Double, Double) = {
    val l_iers2003 =
      (    485868.249036  +
      (1717915923.2178    +
      (        31.8792    +
      (         0.051635  +
      (        -0.00024470)
      * jc) * jc) * jc) * jc) * PI_AS2R % PI2;
    val lp_mhb2000 =
      (  1287104.79305   +
      (129596581.0481    +
      (       -0.5532    +
      (        0.000136  +
      (       -0.00001149)
      * jc) * jc) * jc) * jc) * PI_AS2R % PI2;
    val f_iers2003 =
      (     335779.526232 +
      (1739527262.8478    +
      (       -12.7512    +
      (        -0.001037  +
      (         0.00000417)
      * jc) * jc) * jc) * jc) * PI_AS2R % PI2;
    val d_mhb2000 =
      (   1072260.70369   +
      (1602961601.2090    +
      (        -6.3706    +
      (         0.006593  +
      (        -0.00003169)
      * jc) * jc) * jc) * jc) * PI_AS2R % PI2;
    val om_iers2003 =
      (  450160.398036  +
      (-6962890.5431    +
      (       7.4722    +
      (       0.007702  +
      (      -0.00005939)
      * jc) * jc) * jc) * jc) * PI_AS2R % PI2;

    var dpsi: Double = 0.0;
    var deps: Double = 0.0;
    var i: Int = nutLsData.size - 1;
    while (i >= 0) {
      var d = nutLsData(i);
      val arg = (d(0) * l_iers2003 +
                 d(1) * lp_mhb2000 +
                 d(2) * f_iers2003 +
                 d(3) * d_mhb2000 +
                 d(4) * om_iers2003) % PI2;
      val sin = Math.sin(arg);
      val cos = Math.cos(arg);
      dpsi = dpsi + (d(5) + d(6) * jc) * sin + d(7) * cos;
      deps = deps + (d(8) + d(9) * jc) * cos + d(10) * sin;
      i = i - 1;
    }
    (dpsi, deps);
  }

  // 惑星章動
  private def calcNutationPlanetary(jc: Double): (Double, Double) = {
    val l_mhb2000 =
      (   2.35555598  +
      (8328.6914269554)
      * jc) % PI2;
    val f_mhb2000 =
      (   1.627905234 +
      (8433.466158131)
      * jc) % PI2;
    val d_mhb2000_2 =
      (   5.198466741 +
      (7771.3771468121)
      * jc) % PI2;
    val om_mhb2000 =
      (  2.18243920 +
      (-33.757045)
      * jc) % PI2;
    val pa_iers2003 =
      (0.024381750 +
      (0.00000538691)
      * jc) * jc % PI2;
    val lme_iers2003 =
      (   4.402608842 +
      (2608.7903141574)
      * jc) % PI2;
    val lve_iers2003 =
      (   3.176146697 +
      (1021.3285546211)
      * jc) % PI2;
    val lea_iers2003 =
      (  1.753470314 +
      (628.3075849991)
      * jc) % PI2;
    val lma_iers2003 =
      (  6.203480913 +
      (334.0612426700)
      * jc) % PI2;
    val lju_iers2003 =
      ( 0.599546497 +
      (52.9690962641)
      * jc) % PI2;
    val lsa_iers2003 =
      ( 0.874016757 +
      (21.3299104960)
      * jc) % PI2;
    val lur_iers2003 =
      (5.481293872 +
      (7.4781598567)
      * jc) % PI2;
    val lne_mhb2000 =
      (5.321159000 +
      (3.8127774000)
      * jc) % PI2;


    var dpsi: Double = 0.0;
    var deps: Double = 0.0;
    var i: Int = nutPlData.size - 1;
    while (i >= 0) {
      var d = nutPlData(i);
      val arg = (d(0) * l_mhb2000 +
                 d(2) * f_mhb2000 +
                 d(3) * d_mhb2000_2 +
                 d(4) * om_mhb2000 +
                 d(13) * pa_iers2003 +
                 d(5) * lme_iers2003 +
                 d(6) * lve_iers2003 +
                 d(7) * lea_iers2003 +
                 d(8) * lma_iers2003 +
                 d(9) * lju_iers2003 +
                 d(10) * lsa_iers2003 +
                 d(11) * lur_iers2003 +
                 d(12) * lne_mhb2000) % PI2;
      val sin = Math.sin(arg);
      val cos = Math.cos(arg);
      dpsi = dpsi + d(14) * sin + d(15) * cos;
      deps = deps + d(17) * cos + d(16) * sin;
      i = i - 1;
    }
    (dpsi, deps);
  }

  private def loadNutLs(): IndexedSeq[IndexedSeq[Double]] = {
    // L L' F D Om   PS PST PC EC ECT ES
    val data1 = {
      val sc = new java.util.Scanner(new java.io.BufferedInputStream(new java.io.FileInputStream(nutLsDataPath)));
      val data1 = scala.collection.mutable.ArrayBuffer[String]();
      while (sc.hasNext) {
        data1 += sc.next();
      }
      sc.close();
      data1;
    }
    val count = data1.size / 11;
    val data2: IndexedSeq[IndexedSeq[Double]] = (0 until count).map { i =>
      (0 until 11).map { j =>
        if (j < 5) {
          data1(i * 11 + j).toDouble;
        } else {
          data1(i * 11 + j).toDouble * 0.001 * PI_AS2R;
        }
      }
    }
    data2;
  }

  private def loadNutPl(): IndexedSeq[IndexedSeq[Double]] = {
    // L L' F D Om Lm Lv Le LM Lj Ls Lu Ln Pa   PS PC ES EC
    val data1 = {
      val sc = new java.util.Scanner(new java.io.BufferedInputStream(new java.io.FileInputStream(nutPlDataPath)));
      val data1 = scala.collection.mutable.ArrayBuffer[String]();
      while (sc.hasNext) {
        data1 += sc.next();
      }
      sc.close();
      data1;
    }
    val count = data1.size / 18;
    val data2: IndexedSeq[IndexedSeq[Double]] = (0 until count).map { i =>
      (0 until 18).map { j =>
        if (j < 14) {
          data1(i * 18 + j).toDouble;
        } else {
          data1(i * 18 + j).toDouble * 0.001 * PI_AS2R;
        }
      }
    }
    data2;
  }

  private val nutLsData: IndexedSeq[IndexedSeq[Double]] = loadNutLs();
  private val nutPlData: IndexedSeq[IndexedSeq[Double]] = loadNutPl();

}

class Hcs(lng: Double, lat: Double) {

  def siderealTime(ut1Time: Double): Double = {
    val gmst = TimeLib.mjdut1ToGmst(ut1Time);
    val gast = gmst; // 近似
    val s = gast + lng;
    if (s >= PI2) {
      s - PI2;
    } else {
      s;
    }
  }

  def trueEquatorialXyzToAziAlt(xyz: Array[Double], ut1Time: Double): (Double, Double) = {
    var r: Array[Double] = VectorLib.unitMatrix;
    r = VectorLib.rotateMatrixZ(PI5 - siderealTime(ut1Time), r);
    r = VectorLib.rotateMatrixX(PI5 - lat, r);
    // 地平座標系 X:西 Y:南 Z:天頂
    val hcs = VectorLib.multiplyMV(r, xyz);
    val x = hcs(0);
    val y = hcs(1);
    val z = hcs(2);
    val azi1 = Math.atan2(-x, -y);
    val azi = if (azi1 < 0) azi1 + PI2 else azi1; // 北0°、東90°、南180°、西270°
    val xy = Math.sqrt(x * x + y * y);
    val alt = Math.atan2(z, xy);
    (azi, alt);
  }

}

object Constellations {
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

  def icrsToConstellation(lng: Double, lat: Double): (String, String) = {
    val lng5 = (lng * PI57 / 5).toInt * 5;
    val lat5 = ((lat * PI57 + 90) / 5).toInt * 5 - 90;
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

  def icrsToConstellation(xyz: Array[Double]): (String, String) = {
    val lng = VectorLib.xyzToLng(xyz);
    val lat = VectorLib.xyzToLat(xyz);
    icrsToConstellation(lng, lat);
  }

  val ecliptical = IndexedSeq(
    ("02h00m", "この時期に南の空にいる黄道十二星座はおひつじ座です。西から東に順番に、みずがめ座、うお座、おひつじ座、おうし座、ふたご座、かに座の順に並んでいます"),
    ("03h00m", "この時期に南の空高くにいる黄道十二星座はおうし座です。西から東に順番に、うお座、おひつじ座、おうし座、ふたご座、かに座、の順に並んでいます"),
    ("06h00m", "この時期に南の空高くにいる黄道十二星座はふたご座です。西から東に順番に、おひつじ座、おうし座、ふたご座、かに座、しし座の順に並んでいます"),
    ("08h00m", "この時期に南の空高くにいる黄道十二星座はかに座です。西から東に順番に、おうし座、ふたご座、かに座、しし座、おとめ座の順に並んでいます"),
    ("09h20m", "この時期に南の空にいる黄道十二星座はしし座です。西から東に順番に、ふたご座、かに座、しし座、おとめ座の順に並んでいます"),
    ("12h00m", "この時期に南の空にいる黄道十二星座はおとめ座です。西から東に順番に、ふたご座、かに座、しし座、おとめ座、てんびん座、さそり座の順に並んでいます"),
    ("14h00m", "この時期に南の空にいる黄道十二星座はてんびん座です。西から東に順番に、しし座、おとめ座、てんびん座、さそり座、いて座の順に並んでいます"),
    ("16h00m", "この時期に南の空にいる黄道十二星座はさそり座です。西から東に順番に、おとめ座、てんびん座、さそり座、いて座、やぎ座の順に並んでいます"),
    ("17h40m", "この時期に南の空にいる黄道十二星座はいて座です。西から東に順番に、てんびん座、さそり座、いて座、やぎ座、みずがめ座の順に並んでいます"),
    ("20h00m", "この時期に南の空にいる黄道十二星座はやぎ座です。西から東に順番に、いて座、やぎ座、みずがめ座、うお座の順に並んでいます"),
    ("22h00m", "この時期に南の空にいる黄道十二星座はみずがめ座です。西から東に順番に、やぎ座、みずがめ座、うお座、おひつじ座、おうし座の順に並んでいます"),
    ("23h20m", "この時期に南の空にいる黄道十二星座はうお座です。西から東に順番に、やぎ座、みずがめ座、うお座、おひつじ座、おうし座、ふたご座の順に並んでいます"),
    ("24h00m", ""),
  );

  val northern = IndexedSeq(
    ("00h00m", "この時期の北の空は、北極星の上にカシオペア座がいます"),
    ("02h00m", "この時期の北の空は、北極星の11時の方向(左上)にカシオペア座がいます"),
    ("04h00m", "この時期の北の空は、北極星の10時の方向(左上)にカシオペア座がいます"),
    ("05h00m", "この時期の北の空は、北極星の10時の方向(左上)にカシオペア座、右に北斗七星がいます"),
    ("06h00m", "この時期の北の空は、北極星の左にカシオペア座、右に北斗七星がいます"),
    ("07h00m", "この時期の北の空は、北極星の2時の方向(右上)に北斗七星がいます"),
    ("09h00m", "この時期の北の空は、北極星の1時の方向(右上)に北斗七星がいます"),
    ("11h00m", "この時期の北の空は、北極星の上に北斗七星がいます"),
    ("13h00m", "この時期の北の空は、北極星の11時の方向(左上)に北斗七星がいます"),
    ("15h00m", "この時期の北の空は、北極星の10時の方向(左上)に北斗七星がいます"),
    ("17h00m", "この時期の北の空は、北極星の左に北斗七星がいます"),
    ("18h00m", "この時期の北の空は、北極星の左に北斗七星、右にカシオペア座がいます"),
    ("20h00m", "この時期の北の空は、北極星の2時の方向(右上)にカシオペア座がいます"),
    ("22h00m", "この時期の北の空は、北極星の1時の方向(右上)にカシオペア座がいます"),
    ("24h00m", ""),
  );

  val winter = IndexedSeq(
    ("01h00m", "冬のダイヤモンドは東の空から昇ってきています"),
    ("03h00m", "冬のダイヤモンドは南東の空にいます"),
    ("04h40m", "冬のダイヤモンドは南の空高くにいます"),
    ("07h20m", "冬のダイヤモンドは南西の空にいます"),
    ("09h00m", "冬のダイヤモンドは西の空に沈みかけています"),
    ("11h00m", ""),
    ("24h00m", ""),
  );

  val summer = IndexedSeq(
    ("16h20m", "夏の大三角形は東の空にいます"),
    ("17h20m", "夏の大三角形は東の空高くにいます"),
    ("18h00m", "夏の大三角形は頭上空高くにいます"),
    ("21h00m", "夏の大三角形は西の空高くにいます"),
    ("22h00m", "夏の大三角形は西の空にいます"),
    ("23h20m", ""),
    ("24h00m", ""),
  );

  def siderealTimeToConstellation(siderealTime: Double, data: IndexedSeq[(String, String)]): String = {
    val lng5 = (siderealTime * PI57 / 5).toInt * 5;
    val key = "%02dh%02dm".format(lng5 / 15, lng5 % 15 * 4);
    val p = data.indexWhere(_._1 > key) - 1;
    val cons = if (p == -1) {
      data(data.size - 2)._2;
    } else if (p >= 0) {
      data(p)._2;
    } else {
      throw new Exception();
    }
    cons;
  }
}

object Lib {

  // 外惑星の衝合矩、月の相、月と恒星の合などの探索
  // 単調増加
  // 0～PI2
  // cycle: 1, 2, 4, 8, 24
  def findEvents1(cycle: Int, startTime: Double, endTime: Double, firstStep: Double, stepCountPerDay: Int,
    isDebug: Boolean = false)(f: Double => Double): IndexedSeq[(Double, Int)] = {
    var result: List[(Double, Int)] = Nil;
    var prevTime = startTime;
    var prevValue: Double = f(prevTime) * PI2R * cycle;
    while (prevTime < endTime) {
      if (isDebug) {
        System.err.println("DEBUG %f %f".format(prevTime, prevValue));
      }
      val time = prevTime + firstStep;
      val value = f(time) * PI2R * cycle;
      if (value.toInt - prevValue.toInt != 0) {
        val t = findCrossingBoundaryTime(value.toInt.toDouble * PI2 / cycle, false, true, prevTime, stepCountPerDay, (stepCountPerDay * firstStep).toInt, isDebug)(f);
        result = (t, value.toInt) :: result;
      }
      prevValue = value;
      prevTime = time;
    }
    return result.reverse.toIndexedSeq;
  }

  // 距離などの最大値最小値の探索
  // 2つ目の返値 +1: 最大値
  // 2つ目の返値 -1: 最小値
  def findEvents2(startTime: Double, endTime: Double, firstStep: Double, stepCountPerDay: Int,
    isDebug: Boolean = false)(f: Double => Double): IndexedSeq[(Double, Int)] = {
    findEvents3Sub(startTime, endTime, firstStep, stepCountPerDay, false, isDebug)(f).map { t =>
      (t._1, 2 - t._2);
    }
  }

  // 内惑星の合、最大離角の探索
  // 2つ目の返値 0: 外合
  // 2つ目の返値 1: 東方最大離角
  // 2つ目の返値 2: 内合
  // 2つ目の返値 3: 西方最大離角
  def findEvents3(startTime: Double, endTime: Double, firstStep: Double, stepCountPerDay: Int,
    isDebug: Boolean = false)(f: Double => Double): IndexedSeq[(Double, Int)] = {
    findEvents3Sub(startTime, endTime, firstStep, stepCountPerDay, true, isDebug)(f);
  }

  private def findEvents3Sub(startTime: Double, endTime: Double, firstStep: Double, stepCountPerDay: Int,
    needCrossing: Boolean, isDebug: Boolean = false)(f: Double => Double): IndexedSeq[(Double, Int)] = {
    var result: List[(Double, Int)] = Nil;
    var prevIsIncrease: Int = 0;
    var prevTime = startTime;
    var prevValue: Double = f(prevTime)
    while (prevTime < endTime) {
      val time = prevTime + firstStep;
      val value = f(time);
      val isIncrease = if (value > prevValue) +1 else if (value < prevValue) -1 else 0;
      if (isIncrease * prevIsIncrease < 0) {
        val t = findMaxMinTime(isIncrease > 0, false, prevTime - firstStep, stepCountPerDay, (2 * firstStep * stepCountPerDay).toInt, isDebug)(f);
        result = (t._1, 2 - prevIsIncrease) :: result;
        prevIsIncrease = isIncrease;
      } else if (prevIsIncrease == 0) {
        prevIsIncrease = isIncrease;
      }
      if (needCrossing && prevValue * value < 0.0) {
        val t = findCrossingBoundaryTime(0.0, value < 0.0, false, prevTime, stepCountPerDay, (firstStep * stepCountPerDay).toInt, isDebug)(f);
        if (value > 0.0) {
          result = (t, 0) :: result;
        } else {
          result = (t, 2) :: result;
        }
      }
      prevTime = time;
      prevValue = value;
    }
    return result.reverse.toIndexedSeq;
  }

  def findCrossingBoundaryTime(boundary: Double, isDecrease: Boolean, isCircle: Boolean,
    startTime: Double, stepCountPerDay: Int, maxStep: Int,
    isDebug: Boolean = false)(f: Double => Double): Double = {
    val k = if (isDecrease) -1.0 else 1.0;
    def f2(time: Double): Double = {
      val ret = k * (f(time) - boundary);
      if (isCircle) {
        val z = ret;
        if (z >= PI) {
          ret - PI2;
        } else if (z < -PI) {
          ret + PI2;
        } else {
          ret;
        }
      } else {
        ret;
      }
    }
    if (isDebug) {
      System.err.println("DEBUG Start");
    }
    var minIndex = 0;
    var minValue = f2(startTime);
    var maxIndex = maxStep - 1;
    var maxValue = f2(startTime + maxIndex.toDouble / stepCountPerDay);
    if (isDebug) {
      System.err.println("DEBUG %d %d %f %f".format(minIndex, maxIndex, minValue, maxValue));
    }
    if (minValue > 0 || maxValue <= 0) {
      return -1.0;
    }
    while (minIndex + 1 < maxIndex) {
      val i = (minIndex + maxIndex) / 2;
      val value = f2(startTime + i.toDouble / stepCountPerDay);
      if (value <= 0) {
        minIndex = i;
        minValue = value;
      } else {
        maxIndex = i;
        maxValue = value;
      }
      if (isDebug) {
        System.err.println("DEBUG %d %d %f %f".format(minIndex, maxIndex, minValue, maxValue));
      }
    }
    startTime + (minIndex.toDouble - minValue / (maxValue - minValue)) / stepCountPerDay;
  }

  def findMaxMinTime(isMin: Boolean, isCircle: Boolean,
    startTime: Double, stepCountPerDay: Int, maxStep: Int,
    isDebug: Boolean = false)(f: Double => Double): (Double, Double) = {
    val k = if (isMin) -1.0 else 1.0;
    def f2(time: Double, other: Double): Double = {
      val ret = k * f(time);
      if (isCircle) {
        val z = ret - other;
        if (z >= PI) {
          ret - PI2;
        } else if (z < -PI) {
          ret + PI2;
        } else {
          ret;
        }
      } else {
        ret;
      }
    }
    def findMax(s: Double, t: Double, u: Double): (Double, Double) = {
      val a: Double = (u + s) * 0.5 - t;
      val b: Double = (u - s) * 0.5;
      val maxTime = -0.5 * b / a;
      val maxValue = ((a * maxTime) + b) * maxTime + t;
      (maxTime, maxValue);
    }
    var tIndex = (maxStep - 1) / 2;
    var tValue = f2(startTime + tIndex.toDouble / stepCountPerDay, 0.0);
    var sIndex = 0;
    var sValue = f2(startTime, tValue);
    var uIndex = maxStep - 1;
    var uValue = f2(startTime + uIndex.toDouble / stepCountPerDay, tValue);
    if (isDebug) {
      System.err.println("DEBUG Start");
    }
    while (true) {
      if (isDebug) {
        System.err.println("DEBUG %d %d %d %f %f %f".format(sIndex, tIndex, uIndex, sValue, tValue, uValue));
      }
      if (sValue > tValue) {
        return (-1.0, 0.0);
      } else if (uValue > tValue) {
        return (-1.0, 0.0);
      }
      if (sIndex + 1 == tIndex && tIndex + 1 == uIndex) {
        val (t, v) = findMax(sValue, tValue, uValue);
        val t2 = startTime + (tIndex + t) / stepCountPerDay;
        val v2 = if (isCircle) {
          if (v >= PI2) {
            v - PI2;
          } else if (v < 0.0) {
            v + PI2;
          } else {
            v;
          }
        } else {
          v;
        }
        val v3 = if (isMin) {
          if (isCircle) {
            PI2 - v2;
          } else {
            -v2;
          }
        } else {
          v2;
        }
        return (t2, v3);
      }
      val f = if (sIndex + uIndex < 2 * tIndex) {
        false;
      } else {
        true;
      }
      if (f) {
        val cIndex = (uIndex + tIndex) / 2;
        val cValue = f2(startTime + cIndex.toDouble / stepCountPerDay, tValue);
        if (cValue >= tValue) {
          sIndex = tIndex;
          sValue = tValue;
          tIndex = cIndex;
          tValue = cValue;
        } else {
          uIndex = cIndex;
          uValue = cValue;
        }
      } else {
        val cIndex = (sIndex + tIndex) / 2;
        val cValue = f2(startTime + cIndex.toDouble / stepCountPerDay, tValue);
        if (cValue >= tValue) {
          uIndex = tIndex;
          uValue = tValue;
          tIndex = cIndex;
          tValue = cValue;
        } else {
          sIndex = cIndex;
          sValue = cValue;
        }
      }
    }
    throw new Exception();
  }

}

val jplData = new JplData(jplDataPath);
val hcs = new Hcs(tokyoLng, tokyoLat);

//==============================================================================
// イベント計算
//==============================================================================

def calcSunsetTimes(): IndexedSeq[Double] = {
  val altHor = -0.90 / PI57;
  (0 until period).map { d =>
    Lib.findCrossingBoundaryTime(altHor, true, false, startTime + d + 16.0 / 24.0, 24 * 6, 4 * 6) { time =>
      val utc = time;
      val ut1 = utc; // 近似的
      val tdb = TimeLib.mjdutcToTdb(utc);
      val sun = jplData.calcSunFromEarth(tdb);
      val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
      val sun2 = VectorLib.multiplyMV(bpnMatrix, sun);
      val (azi, alt) = hcs.trueEquatorialXyzToAziAlt(sun2, ut1);
      alt;
    }
  }
}

val moonPhaseTerms: IndexedSeq[(Double, Int)] = { // time, term
  def calcMoonLng(time: Double): Double = {
    val utc = time;
    val tdb = TimeLib.mjdutcToTdb(utc);
    val sun = jplData.calcSunFromEarth(tdb);
    val moon = jplData.calcMoonFromEarth(tdb);
    val bpnMatrix = Bpn.icrsToTrueEclipticMatrix(tdb);
    val sun2 = VectorLib.multiplyMV(bpnMatrix, sun);
    val moon2 = VectorLib.multiplyMV(bpnMatrix, moon);
    val sunLng = VectorLib.xyzToLng(sun2);
    val moonLng = VectorLib.xyzToLng(moon2);
    VectorLib.calcLngDiff(moonLng, sunLng);
  }
  Lib.findEvents1(8, startTime, endTime, 3.0, 24)(calcMoonLng);
}

// TODO 廃止したい
def calcMoonPhaseTerms(): IndexedSeq[(Double, Int, Double)] = { // time, phase, distance
  def calcMoonPhase(time: Double): Double = {
    val utc = time;
    val tdb = TimeLib.mjdutcToTdb(utc);
    val sun = jplData.calcSunFromEarth(tdb);
    val moon = jplData.calcMoonFromEarth(tdb);
    val bpnMatrix = Bpn.icrsToTrueEclipticMatrix(tdb);
    val sun2 = VectorLib.multiplyMV(bpnMatrix, sun);
    val moon2 = VectorLib.multiplyMV(bpnMatrix, moon);
    val sunLng = VectorLib.xyzToLng(sun2);
    val moonLng = VectorLib.xyzToLng(moon2);
    VectorLib.calcLngDiff(moonLng, sunLng);
  }
  def calcMoonDistance(time: Double): Double = {
    val utc = time;
    val tdb = TimeLib.mjdutcToTdb(utc);
    val sun = jplData.calcSunFromEarth(tdb);
    val moon = jplData.calcMoonFromEarth(tdb);
    VectorLib.distance(moon);
  }
  var moonPhaseTerms: List[(Double, Int, Double)] = Nil;
  var currTerm: Int = (calcMoonPhase(startTime) * 8 / PI2).toInt;
  var nextStartTime: Double = startTime - 30.0;
  var nextRange: Int = 5;
  while (nextStartTime >= 0.0 && nextStartTime < endTime + 30.0) {
    val nextTerm = if (currTerm == 7) 0 else currTerm + 1;
    val nextTime = Lib.findCrossingBoundaryTime(nextTerm.toDouble * PI2 / 8, false, true,
      nextStartTime, 24, nextRange * 24) { time =>
      calcMoonPhase(time);
    }
    if (nextTime >= 0.0) {
      val distance = calcMoonDistance(nextTime);
      moonPhaseTerms = (nextTime, nextTerm, distance) :: moonPhaseTerms;
      currTerm = nextTerm;
      nextStartTime = TimeLib.floor(nextTime, 1) + 3.0;
      nextRange = 3;
    } else {
      nextStartTime = -1.0;
    }
  }
  moonPhaseTerms.reverse.toIndexedSeq;
}

val sunsetTimes: IndexedSeq[Double] = calcSunsetTimes();

// TODO 廃止したい
val moonPhaseTerms2: IndexedSeq[(Double, Int, Double)] = calcMoonPhaseTerms();

sealed trait TweetContent {
  def time: Double;
  def message: String;
  def hashtags: List[String];

  def date: String = TimeLib.modifiedJulianDayToStringJSTDate(time);
}

var tweets: Map[String, List[TweetContent]] = Map.empty;

sealed trait OnSunsetTweetContent extends TweetContent {
  def day: Int;
  def message: String;
  def hashtags: List[String];

  def time: Double = sunsetTimes(day);
}

def putTweet(tc: TweetContent): Unit = {
  val date = TimeLib.modifiedJulianDayToStringJSTDate(tc.time);
  val list = tweets.getOrElse(date, Nil);
  tweets = tweets.updated(date, tc :: list);
}

// 日没時の西の空
{
  sealed trait OnSunsetTweetContent2 extends OnSunsetTweetContent {
    def day: Int;
    def message: String;
    def message2: String;
    def hashtags: List[String];
  }
  case class SunsetTweetContent(day: Int) extends OnSunsetTweetContent2 {
    def message: String = "日没は%sごろ".format(TimeLib.modifiedJulianDayToStringJSTNaturalTime(time));
    def message2: String = "日没は%sごろです".format(TimeLib.modifiedJulianDayToStringJSTNaturalTime(time));
    def hashtags: List[String] = Nil;
  }
  case class SunsetMoonTweetContent(day: Int, azi: Double, alt: Double) extends OnSunsetTweetContent2 {
    def azi360: Int = (azi * PI57 + 0.5).toInt;
    def alt360: Int = (alt * PI57 + 0.5).toInt;
    def message: String = "新月後の細い月は日没時に西の空高度約%d°".format(alt360);
    def message2: String = "新月後の細い月は日没時に西の空高度約%d°にいます".format(alt360);
    def hashtags: List[String] = Nil;
  }
  case class SunsetPlanetTweetContent(day: Int, planetName: String,
    azi: Double, alt: Double, isIncreasing: Boolean, isMax: Boolean) extends OnSunsetTweetContent2 {
    def alt360: Int = (alt * PI57 + 0.5).toInt;
    def message: String = if (isMax) {
      "%sは日没時最大高度で西の空高度約%d°".format(planetName, alt360);
    } else if (isIncreasing) {
      "%sは日没時の高度を徐々に上げ、西の空高度約%d°にいます".format(planetName, alt360);
    } else {
      "%sは日没時の高度を徐々に下げ、西の空高度約%d°にいます".format(planetName, alt360);
    }
    def message2: String = if (isMax) {
      "%sは日没時最大高度で西の空高度約%d°です".format(planetName, alt360);
    } else if (isIncreasing) {
      "%sは日没時の高度を徐々に上げ、西の空高度約%d°にいます".format(planetName, alt360);
    } else {
      "%sは日没時の高度を徐々に下げ、西の空高度約%d°にいます".format(planetName, alt360);
    }
    def hashtags: List[String] = if (planetName == "水星" || planetName == "金星") { List(planetName) } else { Nil };
  }
  case class MultiSunsetTweetContent(day: Int, tc: List[OnSunsetTweetContent2]) extends OnSunsetTweetContent {
    def message: String = tc.map(_.message2).mkString("。");
    def hashtags: List[String] = tc.flatMap(_.hashtags);
  }

  val aziThres0 = 200 / PI57;
  val aziThres1 = 315 / PI57;
  val altThres0 = 10 / PI57;
  val Moon = 0;
  val Mercury = 1;
  val Venus = 2;
  val data = new Array[Array[Double]](period);
  val moonFlag = new Array[Boolean](period);
  (0 until period).foreach { day =>
    data(day) = new Array[Double](6);
  }

  val planets = IndexedSeq(("金星", JplData.Venus, 5, Venus), ("水星", JplData.Mercury, 3, Mercury));
  var altData = (0 until period).foreach { day =>
    val time = sunsetTimes(day);
    val ut1 = time; // 近似的
    val tdb = TimeLib.mjdutcToTdb(time);
    val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
    {
      val xyz = jplData.calcMoonFromEarth(tdb);
      val xyz2 = VectorLib.multiplyMV(bpnMatrix, xyz);
      val (azi, alt) = hcs.trueEquatorialXyzToAziAlt(xyz2, ut1);
      data(day)(2 * Moon + 0) = azi;
      data(day)(2 * Moon + 1) = alt;
    }
    (0 until planets.size).foreach { i =>
      val xyz = jplData.calcPlanetFromEarth(tdb, planets(i)._2);
      val xyz2 = VectorLib.multiplyMV(bpnMatrix, xyz);
      val (azi, alt) = hcs.trueEquatorialXyzToAziAlt(xyz2, ut1);
      data(day)(2 * planets(i)._4 + 0) = azi;
      data(day)(2 * planets(i)._4 + 1) = alt;
    }
  }

  (0 until moonPhaseTerms.size).foreach { i =>
    val (moonPhaseTime, term) = moonPhaseTerms(i);
    if (term == 0) {
      val d = (0 until 4).indexWhere { d =>
        val day = (moonPhaseTime + d - startTime).toInt;
        val alt = data(day)(2 * Moon + 1);
        alt >= altThres0;
      }
      val day = (moonPhaseTime + d - startTime).toInt;
      moonFlag(day) = true;
    }
  }

  // 最大高度
  (1 until (period-1)).foreach { day =>
    val dayData = data(day);
    val wday = TimeLib.wday(startTime + day);
    var tc: List[OnSunsetTweetContent2] = Nil;
    planets.foreach { p =>
      val azi = dayData(2 * p._4 + 0);
      val alt = dayData(2 * p._4 + 1);
      val isIncreasing = (alt >= data(day-1)(2 * p._4 + 1));
      val isMax = (isIncreasing && (alt > data(day+1)(2 * p._4 + 1)));
      if (isMax) {
        tc = SunsetPlanetTweetContent(day, p._1, azi, alt, isIncreasing, isMax) :: tc;
      }
    }

    // 新月直後の月
    if (moonFlag(day)) {
      val azi = dayData(2 * Moon + 0);
      val alt = dayData(2 * Moon + 1);
      tc = SunsetMoonTweetContent(day, dayData(2 * Moon + 0), dayData(2 * Moon + 1)) :: tc;
    }

    // 水曜・金曜
    planets.foreach { p =>
      val azi = dayData(2 * p._4 + 0);
      val alt = dayData(2 * p._4 + 1);
      if (azi >= aziThres0 && azi <= aziThres1 && alt >= altThres0 && wday == p._3) {
        if (!tc.exists {
          case SunsetPlanetTweetContent(_, p._1, _, _, _, _) => true;
          case _ => false;
        }) {
          val isIncreasing = (alt >= data(day-1)(2 * p._4 + 1));
          tc = SunsetPlanetTweetContent(day, p._1, azi, alt, isIncreasing, false) :: tc;
        }
      }
    }

    if (!tc.isEmpty) {
      planets.foreach { p =>
        if (!tc.exists {
          case SunsetPlanetTweetContent(_, p._1, _, _, _, _) => true;
          case _ => false;
        }) {
          val azi = dayData(2 * p._4 + 0);
          val alt = dayData(2 * p._4 + 1);
          if (azi >= aziThres0 && azi <= aziThres1 && alt >= altThres0) {
            val isIncreasing = (alt >= data(day-1)(2 * p._4 + 1));
            tc = SunsetPlanetTweetContent(day, p._1, azi, alt, isIncreasing, false) :: tc;
          }
        }
      }
      if (!tc.exists {
        case SunsetMoonTweetContent(_, _, _) => true;
        case _ => false;
      }) {
        val azi = dayData(2 * Moon + 0);
        val alt = dayData(2 * Moon + 1);
        if (azi >= aziThres0 && azi <= aziThres1 && alt >= altThres0) {
          tc = SunsetMoonTweetContent(day, dayData(2 * Moon + 0), dayData(2 * Moon + 1)) :: tc;
        }
      }
    }
    if (wday == 0) {
      tc = SunsetTweetContent(day) :: tc;
    }
    if (!tc.isEmpty) {
      if (!tc.tail.isEmpty) {
        putTweet(MultiSunsetTweetContent(day, tc.reverse));
      } else {
        putTweet(tc.head);
      }
    }
  }
}

// 月相
{
  val termStrs = IndexedSeq(
    "新月",
    "月相 1/8。新月と上弦の中間です \uD83C\uDF12",
    "上弦の月 \uD83C\uDF13",
    "月相 3/8。上弦と満月の中間です \uD83C\uDF14",
    "満月 \uD83C\uDF15",
    "月相 5/8。満月と下弦の中間です \uD83C\uDF16",
    "下弦の月 \uD83C\uDF17",
    "月相 7/8。下弦と新月の中間です \uD83C\uDF18",
  );
  case class MoonPhaseTermTweetContent(rawTime: Double, term: Int) extends TweetContent {
    def time: Double = TimeLib.floor(rawTime, 24) + 1.0 / (24 * 4);
    def message: String = termStrs(term);
    def hashtags: List[String] = Nil;
  }
  case class FullMoonDistanceTweetContent(rawTime: Double, flag: Boolean) extends TweetContent {
    // true: 近い, false: 遠い
    def time: Double = TimeLib.floor(rawTime, 24) + 1.0 / (24 * 4);
    def message: String = if (flag) {
      "満月。月が地球に近く、もっとも大きい満月です";
    } else {
      "満月。月が地球から遠く、もっとも小さい満月です";
    }
    def hashtags: List[String] = Nil;
  }

  // TODO moonPhaseTerms に移行したい
  (8 until (moonPhaseTerms2.size - 8)).foreach { i =>
    val (time, term, distance) = moonPhaseTerms2(i);
    putTweet(if (term == 4) {
      if (moonPhaseTerms2(i - 8)._3 <= distance && distance > moonPhaseTerms2(i + 8)._3) {
        FullMoonDistanceTweetContent(time, false);
      } else if (moonPhaseTerms2(i - 8)._3 >= distance && distance < moonPhaseTerms2(i + 8)._3) {
        FullMoonDistanceTweetContent(time, true);
      } else {
        MoonPhaseTermTweetContent(time, term);
      }
    } else {
      MoonPhaseTermTweetContent(time, term);
    });
  }
}

val sunDistanceEventList: IndexedSeq[(Double, Int)] = Lib.findEvents2(startTime, endTime, 30, 24) { time =>
  val utc = time;
  val tdb = TimeLib.mjdutcToTdb(utc);
  val sun = jplData.calcSunFromEarth(tdb);
  VectorLib.distance(sun);
}

val solarTerms: List[(Double, Int)] = {
  def calcSunLng(time: Double): Double = {
    val utc = time;
    val tdb = TimeLib.mjdutcToTdb(utc);
    val sun = jplData.calcSunFromEarth(tdb);
    val bpnMatrix = Bpn.icrsToTrueEclipticMatrix(tdb);
    val sun2 = VectorLib.multiplyMV(bpnMatrix, sun);
    val sunLng = VectorLib.xyzToLng(sun2);
    sunLng;
  }
  var solarTerms: List[(Double, Int)] = Nil;
  var currTerm: Int = (calcSunLng(startTime) * 24 / PI2).toInt;
  var nextStartTime: Double = startTime;
  var nextRange: Int = 18;
  while (nextStartTime >= 0.0 && nextStartTime < endTime) {
    val nextTerm = if (currTerm == 23) 0 else currTerm + 1;
    val nextTime = Lib.findCrossingBoundaryTime(nextTerm.toDouble * PI2 / 24, false, true,
      nextStartTime, 24, nextRange * 24) { time =>
      calcSunLng(time);
    }
    if (nextTime >= 0.0) {
      solarTerms = (nextTime, nextTerm) :: solarTerms;
      currTerm = nextTerm;
      nextStartTime = TimeLib.floor(nextTime, 1) + 12.0;
      nextRange = 6;
    } else {
      nextStartTime = -1.0;
    }
  }
  solarTerms;
}

val planetPhase: List[(Double, String, String, Boolean, Option[Array[Double]])] = {
  var planetPhase: List[(Double, String, String, Boolean, Option[Array[Double]])] = Nil;
  def calcInnerPlanetLngEc(time: Double, targetPlanet: JplData.TargetPlanet): Double = {
    val utc = time;
    val tdb = TimeLib.mjdutcToTdb(utc);
    val sun = jplData.calcSunFromEarth(tdb);
    val planet = jplData.calcPlanetFromEarth(tdb, targetPlanet);
    val bpnMatrix = Bpn.icrsToTrueEclipticMatrix(tdb);
    val sun2 = VectorLib.multiplyMV(bpnMatrix, sun);
    val planet2 = VectorLib.multiplyMV(bpnMatrix, planet);
    val sunLng = VectorLib.xyzToLng(sun2);
    val planetLng = VectorLib.xyzToLng(planet2);
    val d1 = VectorLib.calcLngDiff(planetLng, sunLng);
    if (d1 >= PI) d1 - PI2 else d1;
  }
  def calcOuterPlanetLngEq(time: Double, targetPlanet: JplData.TargetPlanet): Double = {
    val utc = time;
    val tdb = TimeLib.mjdutcToTdb(utc);
    val sun = jplData.calcSunFromEarth(tdb);
    val planet = jplData.calcPlanetFromEarth(tdb, targetPlanet);
    val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
    val sun2 = VectorLib.multiplyMV(bpnMatrix, sun);
    val planet2 = VectorLib.multiplyMV(bpnMatrix, planet);
    val sunLng = VectorLib.xyzToLng(sun2);
    val planetLng = VectorLib.xyzToLng(planet2);
    val d1 = VectorLib.calcLngDiff(planetLng, sunLng);
    d1;
  }

  List(
    ("火星", JplData.Mars, 779.94, TimeLib.stringToModifiedJulianDay("2020-10-14T00:00:00Z")),
    ("木星", JplData.Jupiter, 398.88, TimeLib.stringToModifiedJulianDay("2020-07-14T00:00:00Z")),
    ("土星", JplData.Saturn, 378.09, TimeLib.stringToModifiedJulianDay("2020-07-21T00:00:00Z")),
  ).foreach { t =>
    val (planetName, targetPlanet, synodicPeriod, conjunctionEpoch) = t;
    val synodicPeriod5 = synodicPeriod * 0.5;

    var conjunctionEvents: List[(Double, Int)] = Nil;
    {
      val n = ((startTime - conjunctionEpoch) / synodicPeriod5).toInt;
      var time = conjunctionEpoch + n * synodicPeriod5;
      var flag = (n % 2) == 0; // true: 衝, false: 合
      while (time - synodicPeriod5 < endTime || conjunctionEvents.size < 2) {
        time = TimeLib.floor(time, 1);
        val th = if (flag) PI else 0;
        time = Lib.findCrossingBoundaryTime(th, true, true,
          time - (synodicPeriod * 0.2).toInt, 24, (synodicPeriod * 0.4).toInt * 24) { time =>
          calcOuterPlanetLngEq(time, targetPlanet);
        }
        if (time < 0.0) {
          throw new Exception();
        }
        conjunctionEvents = (time, (if (flag) 2 else 0)) :: conjunctionEvents;
        time = time + synodicPeriod5;
        flag = !flag;
      }
    }

    {
      val conjunctionEvents2 = conjunctionEvents.reverse.toIndexedSeq;
      (1 until conjunctionEvents2.size).foreach { i =>
        val range1 = conjunctionEvents2(i)._1 - conjunctionEvents2(i-1)._1;
        val time = TimeLib.floor(conjunctionEvents2(i-1)._1 + range1 * 0.2, 1);
        val range = (range1 * 0.6).toInt;
        val term = conjunctionEvents2(i-1)._2 + 1;
        val th = if (term == 1) PI + PI5 else PI5;
        val time2 = Lib.findCrossingBoundaryTime(th, true, true,
          time, 24, range * 24) { time =>
          calcOuterPlanetLngEq(time, targetPlanet);
        }
        if (time2 < 0.0) {
          throw new Exception();
        }
        conjunctionEvents = (time2, term) :: conjunctionEvents;
      }
    }

    conjunctionEvents.foreach { case (time, term) =>
      val utc = time;
      val tdb = TimeLib.mjdutcToTdb(utc);
      val xyz = jplData.calcPlanetFromEarth(tdb, targetPlanet);
      if (term == 0) {
        planetPhase = (time, planetName, "合(赤経基準)", false, None) :: planetPhase;
      } else if (term == 1) {
        planetPhase = (time, planetName, "西矩(赤経基準)", false, Some(xyz)) :: planetPhase;
      } else if (term == 2) {
        planetPhase = (time, planetName, "衝(赤経基準)", false, Some(xyz)) :: planetPhase;
      } else if (term == 3) {
        planetPhase = (time, planetName, "東矩(赤経基準)", false, Some(xyz)) :: planetPhase;
      }
    }
  }

  List(
    ("水星", JplData.Mercury),
    ("金星", JplData.Venus),
  ).foreach { t =>
    val (planetName, targetPlanet) = t;
    Lib.findEvents3(startTime, endTime, 10.0, 24) { time =>
      calcInnerPlanetLngEc(time, targetPlanet);
    }.foreach { case (time, term) =>
      if (term == 0) {
        planetPhase = (time, planetName, "外合(黄経基準)", false, None) :: planetPhase;
      } else if (term == 1) {
        planetPhase = (time, planetName, "東方最大離角(黄経基準) \uD83C\uDF13", true, None) :: planetPhase;
      } else if (term == 2) {
        planetPhase = (time, planetName, "内合(黄経基準)", false, None) :: planetPhase;
      } else if (term == 3) {
        planetPhase = (time, planetName, "西方最大離角(黄経基準) \uD83C\uDF17", true, None) :: planetPhase;
      }
    }
  }

  planetPhase;
}

// 21時・23時の惑星
val (moonCons: List[(Double, Array[Double])], planetCons: List[(Double, String, Array[Double])]) = {
  var moonCons: List[(Double, Array[Double])] = Nil;
  var planetCons: List[(Double, String, Array[Double])] = Nil;
  val altThres = 10 / PI57;
  (0 until period).foreach { d =>
    val time = startTime + d + 21.0 / 24.0;
    val wday = TimeLib.wday(time);
    val ut1 = time; // 近似的
    val tdb = TimeLib.mjdutcToTdb(time);
    val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
    {
      val moon = jplData.calcMoonFromEarth(tdb);
      val moon2 = VectorLib.multiplyMV(bpnMatrix, moon);
      val (azi, alt) = hcs.trueEquatorialXyzToAziAlt(moon2, ut1);
      if (alt >= altThres) {
        moonCons = (time, moon) :: moonCons;
      }
    }
    val (planetName, targetPlanet) = if (wday == 2) {
      ("火星", Some(JplData.Mars));
    } else if (wday == 4) {
      ("木星", Some(JplData.Jupiter));
    } else if (wday == 6) {
      ("土星", Some(JplData.Saturn));
    } else {
      ("", None);
    }
    targetPlanet match {
      case None => ;
      case Some(targetPlanet) =>
        val planet = jplData.calcPlanetFromEarth(tdb, targetPlanet);
        val planet2 = VectorLib.multiplyMV(bpnMatrix, planet);
        val (azi, alt) = hcs.trueEquatorialXyzToAziAlt(planet2, ut1);
        if (alt >= altThres) {
          planetCons = (time, planetName, planet) :: planetCons;
        } else {
          val time = startTime + d + 23.0 / 24.0;
          val ut1 = time; // 近似的
          val tdb = TimeLib.mjdutcToTdb(time);
          val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
          val planet = jplData.calcPlanetFromEarth(tdb, targetPlanet);
          val planet2 = VectorLib.multiplyMV(bpnMatrix, planet);
          val (azi, alt) = hcs.trueEquatorialXyzToAziAlt(planet2, ut1);
          if (alt >= altThres) {
            planetCons = (time, planetName, planet) :: planetCons;
          }
        }
    }
  }
  (moonCons, planetCons);
}

//==============================================================================
// ツイート文字列構築
//==============================================================================

case class LegacyTweetContent(time: Double, message: String) extends  TweetContent {
  def hashtags: List[String] = Nil;
}

def putTweet(time: Double, msg: String): Unit = {
  putTweet(LegacyTweetContent(time, msg));
}

sunDistanceEventList.foreach { case (time, flag) =>
  val s = if (flag < 0) "近日点" else "遠日点";
  putTweet(TimeLib.round(time, 24) - 1.0 / (24 * 4), "地球が%s通過".format(s));
}

{
  val termStrs = IndexedSeq(
    "春分", "清明", "穀雨", "立夏", "小満", "芒種", "夏至", "小暑", "大暑", "立秋", "処暑", "白露",
    "秋分", "寒露", "霜降", "立冬", "小雪", "大雪", "冬至", "小寒", "大寒", "立春", "雨水", "啓蟄",
  );
  solarTerms.foreach { case (time, term) =>
    if (term % 6 == 0) {
      putTweet(TimeLib.floor(time, 24) + 1.0 / (24 * 4), "%s。太陽の黄経が%d°です".format(termStrs(term), term * 15));
    } else {
      putTweet(TimeLib.floor(time, 24) + 1.0 / (24 * 4), "二十四節気の%s。太陽の黄経が%d°です".format(termStrs(term), term * 15));
    }
  }
}

planetPhase.foreach { case (time, planetName, content, timeFlag, xyzOpt) =>
  val time2 = if (timeFlag) {
    TimeLib.round(time, 24) - 1.0 / (24 * 4);
  } else {
    TimeLib.floor(time, 24) + 1.0 / (24 * 4);
  }
  xyzOpt match {
  case None =>
    putTweet(time2, "%sが%s #%s".format(planetName, content, planetName));
  case Some(xyz) =>
    val (conscomment, cons) = Constellations.icrsToConstellation(xyz);
    putTweet(time2, "%s%sが%s。%sにいます #%s".format(conscomment, planetName, content, cons, planetName));
  }
}

moonCons.foreach { case (time, xyz) =>
  val (conscomment, cons) = Constellations.icrsToConstellation(xyz);
  putTweet(TimeLib.floor(time, 24) - 1.0 / (24 * 6), "%s月は%sにいます".format(conscomment, cons));
}
planetCons.foreach { case (time, planetName, xyz) =>
  val (conscomment, cons) = Constellations.icrsToConstellation(xyz);
  putTweet(TimeLib.floor(time, 24) - 1.0 / (24 * 4), "%s%sは%sにいます #%s".format(conscomment, planetName, cons, planetName));
}

//==============================================================================

{
  case class StarTweetContent(time: Double, message: String) extends TweetContent {
    def hashtags: List[String] = List("星空");
  }
  val stars = IndexedSeq(
    // J2000
    ("00h11m", "ペガスス座の四角形の左の辺が真南です"),
    ("03h47m", "おうし座すばるは天頂付近にいます"),
    ("04h36m", "おうし座アルデバランは真南の空高くにいます。冬のダイヤモンドを構成する6個の星の1つです"),
    ("05h15m", "オリオン座リゲルは真南にいます。冬のダイヤモンドを構成する6個の星の1つです"),
    ("05h17m", "ぎょしゃ座カペラは天頂付近にいます。冬のダイヤモンドを構成する6個の星の1つです"),
    ("05h55m", "オリオン座ベテルギウスは真南にいます。冬のダイヤモンドの中心の星です"),
    ("06h45m", "おおいぬ座シリウスは真南にいます。太陽を除いてもっとも明るい恒星です。冬のダイヤモンドを構成する6個の星の1つです"),
    ("07h35m", "ふたご座の右側の明るい星カストルは天頂付近にいます。左側のポルックスよりよりやや暗いです"),
    ("07h39m", "こいぬ座プロキオンは真南にいます。冬のダイヤモンドを構成する6個の星の1つです"),
    ("07h45m", "ふたご座の左側の明るい星ポルックスは天頂付近にいます。冬のダイヤモンドを構成する6個の星の1つです"),
    ("10h08m", "しし座の前足にあるレグルスは真南にいます"),
    ("11h49m", "しし座の尾にあるデネボラは真南にいます"),
    ("13h25m", "おとめ座スピカは真南にいます"),
    ("14h16m", "うしかい座アークトゥルスは真南の空高くにいます。シリウス、カノープスに次ぐ明るさの恒星です"),
    ("16h29m", "さそり座アンタレスは真南にいます"),
    ("18h37m", "こと座ベガは天頂付近にいます。夏の大三角形の1つです"),
    ("19h51m", "わし座アルタイルは真南にいます。夏の大三角形であるベガ・アルタイル・デネブは空高くにいます"),
    ("20h41m", "はくちょう座デネブは天頂付近にいます。夏の大三角形の1つです"),
    ("23h04m", "ペガスス座の四角形の右の辺が真南です"),
  ).map { t =>
    ((t._1.substring(0, 2).toInt.toDouble + t._1.substring(3, 5).toInt.toDouble / 60) / 24 * PI2, t._2);
  }
  var index = -1;
  (0 until period).foreach { day =>
    val date = TimeLib.modifiedJulianDayToStringJSTDate(startTime + day);
    val time = startTime + day + 21.0 / 24.0;
    val sid = hcs.siderealTime(time);
    if (index < 0) {
      index = stars.indexWhere(_._1 > sid);
    } else {
      if (sid > stars(index)._1 && (index > 0 || sid < PI)) {
        putTweet(StarTweetContent(time, stars(index)._2));
        index += 1;
        if (index == stars.size) {
          index = 0;
        }
      }
    }
  }
}

def tweetConstellations(data: IndexedSeq[(String, String)], span: Int): Unit = {
  val hashtag = " #星空";
  var nextDay: Int = 11;
  (0 until period).foreach { day =>
    val date = TimeLib.modifiedJulianDayToStringJSTDate(startTime + day);
    if (day >= nextDay && !tweets.contains(date)) {
      val time = startTime + day + 21.0 / 24.0;
      val sid = hcs.siderealTime(time);
      val cons = Constellations.siderealTimeToConstellation(sid, data);
      if (cons != "") {
        putTweet(time, cons + hashtag);
      }
      nextDay = day + span;
    }
  }
}
tweetConstellations(Constellations.ecliptical, 14);
tweetConstellations(Constellations.winter, 7);
tweetConstellations(Constellations.summer, 7);
tweetConstellations(Constellations.northern, 14);

// なにもツイートのない日付をエラー出力
{
  (0 until period).foreach { day =>
    val date = TimeLib.modifiedJulianDayToStringJSTDate(startTime + day);
    if (!tweets.contains(date)) {
      val time = startTime + day + 21.0 / 24.0;
      val sid = hcs.siderealTime(time);
      System.err.println("%s #empty".format(date));
    }
  }
}

//==============================================================================
// ツイート出力

(0 until period).foreach { d =>
  val date = TimeLib.modifiedJulianDayToStringJSTDate(startTime + d);
  if (tweets.contains(date)) {
    import Ordering.Double.IeeeOrdering;
    tweets(date).sortBy(_.time).foreach { case tc =>
      val time = tc.time;
      if (time >= startTime && time < endTime) {
        val msg = tc.message + tc.hashtags.map(" #" + _).mkString;
        println("%s %s".format(TimeLib.modifiedJulianDayToStringJST(time), msg));
      }
    }
  }
}

//==============================================================================
