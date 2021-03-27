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
  object Moon    extends TargetPlanet;
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
      case JplData.Moon    => return calcMoonFromEarth(time);
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

  def angularDistance1(xyz1: Array[Double], xyz2: Array[Double]): Double = {
    val x1 = xyz1(0);
    val y1 = xyz1(1);
    val z1 = xyz1(2);
    val x2 = xyz2(0);
    val y2 = xyz2(1);
    val z2 = xyz2(2);
    (x1 * x2 + y1 * y2 + z1 * z2) / Math.sqrt((x1 * x1 + y1 * y1 + z1 * z1) * (x2 * x2 + y2 * y2 + z2 * z2));
  }

  def angularDistance(xyz1: Array[Double], xyz2: Array[Double]): Double = {
    Math.acos(angularDistance1(xyz1, xyz2));
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
    ret(6) = -sin * r(0) + cos * r(6);
    ret(7) = -sin * r(1) + cos * r(7);
    ret(8) = -sin * r(2) + cos * r(8);
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

  val icrsToJ2000Matrix: Array[Double] = {
    var r: Array[Double] = VectorLib.unitMatrix;
    r = VectorLib.rotateMatrixX(  5.1 / 1000 / 3600 / PI57, r);
    r = VectorLib.rotateMatrixY( 17.3 / 1000 / 3600 / PI57, r);
    r = VectorLib.rotateMatrixZ(-78.0 / 1000 / 3600 / PI57, r);
    r;
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
  val constellationsMap = Map[String, (String, List[String])] (
    " 0h00m, -5" -> ("うお座の西側の魚(ペガスス座の南)の南", Nil),

    " 0h20m, -5" -> ("くじら座のしっぽ付近", Nil),

    " 0h40m,  0" -> ("うお座の西側の魚のしっぽ付近", Nil),

    " 1h00m,  0" -> ("うお座の西側の魚(ペガスス座の南)のしっぽ付近", Nil),

    " 1h20m,  5" -> ("うお座の2匹の魚が接続されている付近", Nil),
    " 1h20m, 10" -> ("うお座の東側の魚(アンドロメダ座の南)のしっぽ付近", Nil),

    " 1h40m,  5" -> ("うお座の2匹の魚の付け根付近", Nil),
    " 1h40m, 10" -> ("うお座の東側の魚(アンドロメダ座の南)のしっぽ付近", Nil),

    " 2h00m,  5" -> ("くじら座の頭付近", Nil),
    " 2h00m, 10" -> ("おひつじ座の頭とくじら座の頭の間", Nil),
    " 2h00m, 15" -> ("おひつじ座の頭付近", Nil),

    " 2h20m, 10" -> ("おひつじ座とくじら座の頭の間", Nil),
    " 2h20m, 15" -> ("おひつじ座の南", Nil),

    " 2h40m, 10" -> ("おひつじ座とくじら座の頭の間", Nil),
    " 2h40m, 15" -> ("おひつじ座の南", Nil),

    " 3h00m, 15" -> ("おひつじ座のしっぽ側", Nil),

    " 3h20m, 15" -> ("おうし座の西", Nil),
    " 3h20m, 20" -> ("おうし座すばる付近", List("プレアデス星団")),

    " 3h40m, 15" -> ("おうし座すばるの南", List("プレアデス星団")),
    " 3h40m, 20" -> ("おうし座すばる付近", List("プレアデス星団")),

    " 4h00m, 15" -> ("おうし座ヒアデス星団の西", Nil),
    " 4h00m, 20" -> ("おうし座すばる付近", List("プレアデス星団")),

    " 4h20m, 20" -> ("おうし座とペルセウス座とぎょしゃ座の間", Nil),

    " 4h40m, 20" -> ("おうし座ヒアデス星団の北東", Nil),

    " 5h00m, 20" -> ("おうし座の角付近", Nil),

    " 5h20m, 20" -> ("おうし座の角とぎょしゃ座付近", Nil),

    " 5h40m, 20" -> ("おうし座の角とオリオン座の腕付近", List("オリオン座")),

    " 6h00m, 20" -> ("ふたご座の西側の子の足元とオリオン座の腕付近", List("オリオン座")),

    " 6h20m, 20" -> ("ふたご座の西側の子の足元付近", Nil),
    " 6h20m, 25" -> ("ふたご座の西側の子の腰付近", Nil),

    " 6h40m, 20" -> ("ふたご座の2人の足元付近", Nil),
    " 6h40m, 25" -> ("ふたご座の西側の子の胴体付近", Nil),

    " 7h00m, 20" -> ("ふたご座の東側の子の腰付近", Nil),
    " 7h00m, 25" -> ("ふたご座の西側の子の胴体付近", Nil),

    " 7h20m, 20" -> ("ふたご座の東側の子の胴体付近", Nil),
    " 7h20m, 25" -> ("ふたご座の東側の子の胴体付近", Nil),

    " 7h40m, 20" -> ("ふたご座ポルックスの南でふたご座の東", Nil),

    " 8h00m, 20" -> ("かに座の西", Nil),

    " 8h20m, 20" -> ("かに座", Nil),

    " 8h40m, 20" -> ("かに座", Nil),

    " 9h00m, 20" -> ("かに座の東", Nil),

    " 9h20m, 15" -> ("しし座とかに座の間", Nil),

    " 9h40m, 15" -> ("しし座の西", Nil),

    "10h00m, 15" -> ("しし座の肩付近", Nil),

    "10h20m, 10" -> ("しし座レグルスの東", Nil),
    "10h20m, 15" -> ("しし座の中央", Nil),

    "10h40m, 10" -> ("しし座の腹付近", Nil),

    "11h00m, 10" -> ("しし座の後ろ足付近", Nil),

    "11h20m,  5" -> ("しし座の後ろ足とおとめ座の間", Nil),

    "11h40m,  5" -> ("おとめ座の頭付近", Nil),

    "12h00m,  0" -> ("おとめ座の四角形の西", Nil),
    "12h00m,  5" -> ("おとめ座の頭付近", Nil),

    "12h20m,  0" -> ("おとめ座の四角形の西", Nil),

    "12h40m,  0" -> ("おとめ座の四角形", Nil),
    "12h40m, -5" -> ("おとめ座の四角形", Nil),

    "13h00m, -5" -> ("おとめ座の四角形", Nil),

    "13h20m, -5" -> ("おとめ座スピカの北", Nil),
    "13h20m,-10" -> ("おとめ座スピカの北", Nil),

    "13h40m,-10" -> ("おとめ座スピカの東", Nil),

    "14h00m,-10" -> ("おとめ座スピカの東", Nil),
    "14h00m,-15" -> ("おとめ座スピカとてんびん座の間", Nil),

    "14h20m,-15" -> ("てんびん座の西おとめ座との間", Nil),

    "14h40m,-15" -> ("てんびん座の西", Nil),

    "15h00m,-20" -> ("てんびん座", Nil),

    "15h20m,-20" -> ("てんびん座", Nil),

    "15h40m,-20" -> ("てんびん座の東さそり座との間", Nil),

    "16h00m,-25" -> ("さそり座アンタレスの北西", Nil),

    "16h20m,-25" -> ("へびつかい座でさそり座アンタレスの北", Nil),

    "16h40m,-25" -> ("へびつかい座でさそり座アンタレスの北東", Nil),

    "17h00m,-25" -> ("へびつかい座でさそり座アンタレスの東", Nil),

    "17h20m,-25" -> ("へびつかい座でさそり座アンタレスの東", Nil),

    "17h40m,-30" -> ("いて座さそり座へびつかい座の間", Nil),

    "18h00m,-30" -> ("いて座の西", Nil),

    "18h20m,-30" -> ("いて座南斗六星", Nil),

    "18h40m,-30" -> ("いて座南斗六星", Nil),

    "19h00m,-30" -> ("いて座南斗六星の先", Nil),

    "19h20m,-30" -> ("いて座南斗六星の東", Nil),

    "20h00m,-25" -> ("やぎ座の西部", Nil),

    "20h20m,-20" -> ("やぎ座の西部", Nil),

    "20h40m,-20" -> ("やぎ座の中央", Nil),
    "20h40m,-25" -> ("やぎ座の中央", Nil),

    "21h00m,-20" -> ("やぎ座の中央", Nil),
    "21h00m,-25" -> ("やぎ座の中央", Nil),

    "21h20m,-20" -> ("やぎ座の東部", Nil),

    "21h40m,-15" -> ("やぎ座とみずがめ座の間", Nil),
    "21h40m,-20" -> ("やぎ座の東部", Nil),

    "22h00m,-20" -> ("みずがめ座の南部でやぎ座の東", Nil),
    "22h00m,-15" -> ("みずがめ座の南部でやぎ座の東", Nil),

    "22h20m,-15" -> ("みずがめ座の南部でやぎ座の東", Nil),
    "22h20m,-20" -> ("みずがめ座の南部でやぎ座の東", Nil),

    "22h40m,-15" -> ("みずがめ座の南部でうお座の西側の魚の頭の南西", Nil),

    "23h00m,-15" -> ("みずがめ座の南部でうお座の西側の魚の頭の南", Nil),

    "23h20m,-10" -> ("みずがめ座とうお座の西側の魚の頭の間", Nil),

    "23h40m,-10" -> ("みずがめ座とうお座の西側の魚の頭の間", Nil),
  );

  def icrsToConstellation(lng: Double, lat: Double): (String, String, List[String]) = {
    val lng5 = (lng * PI57 / 5).toInt * 5;
    val lat5 = ((lat * PI57 + 90) / 5).toInt * 5 - 90;
    val key = "%2dh%02dm,%3d".format(lng5 / 15, lng5 % 15 * 4, lat5);
    val (cons, hashtags) = constellationsMap.getOrElse(key, ("-", Nil));
    if (cons == "") {
      ("#", "(%s)".format(key), hashtags);
    } else if (cons == "-") {
      ("##", "(%s)".format(key), hashtags);
    } else {
      ("", cons, hashtags);
    }
  }

  def icrsToConstellation(xyz: Array[Double]): (String, String, List[String]) = {
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

  val season = IndexedSeq(
    ("01h00m", "冬のダイヤモンドは東の空から昇ってきています"),
    ("03h00m", "冬のダイヤモンドは南東の空にいます"),
    ("04h40m", "冬のダイヤモンドは南の空高くにいます"),
    ("07h20m", "冬のダイヤモンドは南西の空にいます"),
    ("09h00m", "冬のダイヤモンドは西の空に沈みかけています。東の空には春の大三角形がいます"),
    ("10h00m", ""),
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

object MathLib {

  def circleAdd(a: Double, b: Double): Double = {
    val d = a + b;
    if (d > PI) {
      d - PI2;
    } else if (d > -PI) {
      d;
    } else {
      d + PI2;
    }
  }

  // 2つ目の返値 +1: 最大値
  // 2つ目の返値 -1: 最小値
  def findMaxMinListContinuous(start: Double, end: Double, step: Double, partitionCount: Int)(f: Double => Double): IndexedSeq[(Double, Int)] = {
    findMaxMinListContinuous(0, ((end - start) * partitionCount).toInt, (step * partitionCount).toInt) { x =>
      f(x / partitionCount + start);
    }.map { case (x, p) =>
      (x / partitionCount + start, p);
    }
  }

  // 2つ目の返値 +1: 最大値
  // 2つ目の返値 -1: 最小値
  def findMaxMinListContinuous(start: Int, end: Int, step: Int)(f: Double => Double): IndexedSeq[(Double, Int)] = {
    findMaxMinListDiscrete(start, end, step)(x => f(x.toDouble)).map { case (x, p) =>
      (x + findExtremum(f(x - 1), f(x), f(x + 1)), p);
    }
  }

  // 返値 0: 増加中
  // 返値 1: 最大値
  // 返値 2: 減少中
  // 返値 3: 最小値
  def getMaxMinUpDownFlagListDiscrete(start: Int, end: Int, step: Int)(f: Int => Double): IndexedSeq[Int] = {
    val maxMinList = findMaxMinListDiscrete(start, end, step)(f);
    if (maxMinList.size == 0) {
      if (f(start) > f(end)) {
        IndexedSeq.fill(end - start)(2);
      } else {
        IndexedSeq.fill(end - start)(0);
      }
    } else {
      var listIdx: Int = 0;
      var nextIdx: Int = maxMinList(0)._1;
      var flag: Int = maxMinList(0)._2;
      (start until end).map { i =>
        if (i == nextIdx) {
          val value = 2 - flag;
          listIdx += 1;
          if (listIdx == maxMinList.size) {
            nextIdx = end;
            flag = - flag;
          } else {
            nextIdx = maxMinList(listIdx)._1;
            flag = maxMinList(listIdx)._2;
          }
          value;
        } else {
          1 - flag;
        }
      }
    }
  }

  // 2つ目の返値 +1: 最大値
  // 2つ目の返値 -1: 最小値
  def findMaxMinListDiscrete(start: Int, end: Int, step: Int)(f: Int => Double): IndexedSeq[(Int, Int)] = {
    val fm = { x: Int => -f(x) }
    var result: List[(Int, Int)] = Nil;
    var s: Int = start;
    var sv: Double = f(s);
    var t: Int = s + step;
    if (t >= end) {
      return IndexedSeq.empty;
    }
    var tv: Double = f(t);
    var u: Int = t + step;
    while (u < end) {
      val uv = f(u);
      if (sv < tv && tv > uv) {
        val max = findMaxDiscrete(s, t, u, sv, tv, uv)(f);
        result = (max, +1) :: result;
      } else if (sv > tv && tv < uv) {
        val min = findMaxDiscrete(s, t, u, -sv, -tv, -uv)(fm);
        result = (min, -1) :: result;
      }
      s = t;
      sv = tv;
      t = u;
      tv = uv;
      u = u + step;
    }
    result.reverse.toIndexedSeq;
  }

  // 2つ目の返値 0: 0
  // 2つ目の返値 1: 最大値
  // 2つ目の返値 2: 0
  // 2つ目の返値 3: 最小値
  def findMaxMinCrossingListContinuous(start: Double, end: Double, step: Double, partitionCount: Int)(f: Double => Double): IndexedSeq[(Double, Int)] = {
    findMaxMinCrossingListContinuous(0, ((end - start) * partitionCount).toInt, (step * partitionCount).toInt) { x =>
      f(x / partitionCount + start);
    }.map { case (x, p) =>
      (x / partitionCount + start, p);
    }
  }

  // 2つ目の返値 0: 0
  // 2つ目の返値 1: 最大値
  // 2つ目の返値 2: 0
  // 2つ目の返値 3: 最小値
  def findMaxMinCrossingListContinuous(start: Int, end: Int, step: Int)(f: Double => Double): IndexedSeq[(Double, Int)] = {
    findMaxMinCrossingListDiscrete(start, end, step)(x => f(x.toDouble)).map { case (x, p) =>
      if (p % 2 != 0) {
        (x + findExtremum(f(x - 1), f(x), f(x + 1)), p);
      } else if (p == 0) {
        (x + findCrossing(f(x), f(x + 1)), p);
      } else {
        (x + findCrossing(-f(x), -f(x + 1)), p);
      }
    }
  }

  // 2つ目の返値 0: 0
  // 2つ目の返値 1: 最大値
  // 2つ目の返値 2: 0
  // 2つ目の返値 3: 最小値
  def findMaxMinCrossingListDiscrete(start: Int, end: Int, step: Int)(f: Int => Double): IndexedSeq[(Int, Int)] = {
    val fm = { x: Int => -f(x) }
    var result: List[(Int, Int)] = Nil;
    var s: Int = start;
    var sv: Double = f(s);
    var t: Int = s + step;
    if (t >= end) {
      return IndexedSeq.empty;
    }
    var tv: Double = f(t);
    if (sv <= 0.0 && tv > 0.0) {
      val cr = findCrossingDiscrete(s, t, sv, tv)(f);
      result = (cr, 0) :: result;
    } else if (sv >= 0.0 && tv < 0.0) {
      val cr = findCrossingDiscrete(s, t, -sv, -tv)(fm);
      result = (cr, 2) :: result;
    }
    var u: Int = t + step;
    while (u < end) {
      val uv = f(u);
      if (sv < tv && tv > uv) {
        val max = findMaxDiscrete(s, t, u, sv, tv, uv)(f);
        result = (max, 1) :: result;
      } else if (sv > tv && tv < uv) {
        val min = findMaxDiscrete(s, t, u, -sv, -tv, -uv)(x => -f(x));
        result = (min, 3) :: result;
      }
      if (tv <= 0.0 && uv > 0.0) {
        val cr = findCrossingDiscrete(t, u, tv, uv)(f);
        result = (cr, 0) :: result;
      } else if (tv >= 0.0 && uv < 0.0) {
        val cr = findCrossingDiscrete(t, u, -tv, -uv)(fm);
        result = (cr, 2) :: result;
      }
      s = t;
      sv = tv;
      t = u;
      tv = uv;
      u = u + step;
    }
    result.reverse.toIndexedSeq;
  }

  def findCyclicPhaseListContinuous(cycle: Int, start: Double, end: Double, step: Double, partitionCount: Int)(f: Double => Double): IndexedSeq[(Double, Int)] = {
    findCyclicPhaseListContinuous(cycle, 0, ((end - start) * partitionCount).toInt, (step * partitionCount).toInt) { x =>
      f(x / partitionCount + start);
    }.map { case (x, p) =>
      (x / partitionCount + start, p);
    }
  }

  def findCyclicPhaseListContinuous(cycle: Int, start: Int, end: Int, step: Int)(f: Double => Double): IndexedSeq[(Double, Int)] = {
    val cycle_r = PI2R * cycle;
    val cycle2 = cycle * 0.5;
    val fs = (0 until cycle).map { i =>
      { x: Int =>
        val y = f(x) * cycle_r - i;
        if (y >= cycle2) {
          y - cycle;
        } else {
          y;
        }
      }
    }
    val fs0 = fs(0);
    var result: List[(Double, Int)] = Nil;
    var prevX: Int = start;
    var nextI: Int = (f(prevX) * cycle_r).toInt + 1;
    if (nextI == cycle) nextI = 0;
    var fsi = fs(nextI);
    var prevValue: Double = fsi(prevX)
    while (prevX + 1 < end) {
      val x = (prevX + step) min (end - 1);
      val value = fsi(x);
      if (value > 0.0) {
        val t = findCrossingDiscrete(prevX, x, prevValue, value)(fsi);
        val t2 = t + findCrossing(fsi(t), fsi(t + 1));
        result = (t2, nextI) :: result;
        nextI += 1;
        if (nextI == cycle) nextI = 0;
        fsi = fs(nextI);
        prevValue = fsi(x);
      } else {
        prevValue = value;
      }
      prevX = x;
    }
    return result.reverse.toIndexedSeq;
  }

  @scala.annotation.tailrec
  private def findMaxDiscrete(s: Int, t: Int, u: Int, sv: Double, tv: Double, uv: Double)(f: Int => Double): Int = {
    if (s + 1 == t && t + 1 == u) {
      t;
    } else if (s + u >= 2 * t) {
      val c = (u + t) / 2;
      val cv = f(c);
      if (cv >= tv) {
        findMaxDiscrete(t, c, u, tv, cv, uv)(f);
      } else {
        findMaxDiscrete(s, t, c, sv, tv, cv)(f);
      }
    } else {
      val c = (s + t) / 2;
      val cv = f(c);
      if (cv > tv) {
        findMaxDiscrete(s, c, t, sv, cv, tv)(f);
      } else {
        findMaxDiscrete(c, t, u, cv, tv, uv)(f);
      }
    }
  }

  @scala.annotation.tailrec
  private def findCrossingDiscrete(s: Int, t: Int, sv: Double, tv: Double)(f: Int => Double): Int = {
    if (s + 1 == t) {
      s;
    } else {
      val c = (s + t) / 2;
      val cv = f(c);
      if (cv > 0.0) {
        findCrossingDiscrete(s, c, sv, cv)(f);
      } else {
        findCrossingDiscrete(c, t, cv, tv)(f);
      }
    }
  }

  private def findExtremum(sv: Double, tv: Double, uv: Double): Double = {
    val a: Double = (uv + sv) * 0.5 - tv;
    val b: Double = (uv - sv) * 0.5;
    -0.5 * b / a;
  }

  private def findCrossing(sv: Double, tv: Double): Double = {
    -sv / (tv - sv);
  }
}

object Lib2 {

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

val sunsetTimesData: IndexedSeq[(Double, Double, Array[Double])] = { // time, tdb, bpnMatrix
  val altHor = -0.90 / PI57;
  (0 until period).map { d =>
    val time = Lib2.findCrossingBoundaryTime(altHor, true, false, startTime + d + 16.0 / 24.0, 24 * 6, 4 * 6) { time =>
      val utc = time;
      val ut1 = utc; // 近似的
      val tdb = TimeLib.mjdutcToTdb(utc);
      val sun = jplData.calcSunFromEarth(tdb);
      val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
      val sun2 = VectorLib.multiplyMV(bpnMatrix, sun);
      val (azi, alt) = hcs.trueEquatorialXyzToAziAlt(sun2, ut1);
      alt;
    }
    val tdb = TimeLib.mjdutcToTdb(time);
    val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
    (time, tdb, bpnMatrix);
  }
}

def sunsetTimes(day: Int): Double = sunsetTimesData(day)._1;

def calcPlanetOnSunsetTime(day: Int, targetPlanet: JplData.TargetPlanet): (Double, Double) = { // azi, alt
  val (time, tdb, bpnMatrix) = sunsetTimesData(day);
  val ut1 = time; // 近似的
  val xyz = jplData.calcPlanetFromEarth(tdb, targetPlanet);
  val xyz2 = VectorLib.multiplyMV(bpnMatrix, xyz);
  hcs.trueEquatorialXyzToAziAlt(xyz2, ut1);
}

def calcPlanetXyzAziAlt(time: Double, targetPlanet: JplData.TargetPlanet): (Array[Double], Double, Double) = {
  val utc = time;
  val tdb = TimeLib.mjdutcToTdb(utc);
  val ut1 = time; // 近似的
  val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
  val xyz = jplData.calcPlanetFromEarth(tdb, targetPlanet);
  val xyz2 = VectorLib.multiplyMV(bpnMatrix, xyz);
  val azialt = hcs.trueEquatorialXyzToAziAlt(xyz2, ut1);
  (xyz, azialt._1, azialt._2);
}

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
def calcMoonPhase(time: Double): Double = {
  calcMoonLng(time) / PI2 * 28.0;
}
val moonPhaseTerms: IndexedSeq[(Double, Int)] = { // time, term
  MathLib.findCyclicPhaseListContinuous(8, startTime, endTime, 3.0, 24)(calcMoonLng);
}

//==============================================================================

def isNightTime0(time: Double): Boolean = {
  val day = (time - startTime).toInt;
  val s = sunsetTimes(day);
  time - time.toInt >= 6.0 / 24 && time - time.toInt < 15.0 / 24;
}

def isNightTime1(time: Double): Boolean = {
  val day = (time - startTime).toInt;
  val s = sunsetTimes(day);
  time >= s && time - time.toInt < 15.0 / 24;
}

def isNightTime2(time: Double): Boolean = {
  val day = (time - startTime).toInt;
  val s = sunsetTimes(day);
  time >= s && time - time.toInt < 16.0 / 24;
}

def getConjunctionTweetTime(time: Double, xyz: Array[Double]): Option[Double] = {
  val altThres0 = 10 / PI57;
  val sun_distanceThres = 20.0 / PI57;
  val utc = time;
  val ut1 = time; // 近似的
  val tdb = TimeLib.mjdutcToTdb(utc);
  val sun_xyz = jplData.calcSunFromEarth(tdb);
  val sun_distance = VectorLib.angularDistance(sun_xyz, xyz);
  if (sun_distance < sun_distanceThres) {
    None;
  } else {
    val time21 = time.toInt + 0.5;
    val utc21 = time21;
    val ut121 = utc21; // 近似的
    val tdb21 = TimeLib.mjdutcToTdb(utc21);
    val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb21);
    val xyz2 = VectorLib.multiplyMV(bpnMatrix, xyz);
    val (azi, alt) = hcs.trueEquatorialXyzToAziAlt(xyz2, ut121);
    if (alt >= altThres0) {
      // 21時に見える
      Some(time21);
    } else {
      val (azi, alt) = hcs.trueEquatorialXyzToAziAlt(xyz2, ut121 + 2.0 / 24);
      if (alt >= altThres0) {
        // 23時に見える
        Some(time21 + 2.0 / 24);
      } else {
        val day = (time - startTime).toInt;
        val sunsetTime = sunsetTimes(day);
        val (azi, alt) = hcs.trueEquatorialXyzToAziAlt(xyz2, sunsetTime);
        if (alt >= altThres0) {
          // 日没時に見える
          Some(sunsetTime);
        } else {
          None;
        }
      }
    }
  }
}

//==============================================================================
// ツイート管理
//==============================================================================

sealed trait TweetContent {
  def time: Double;
  def message: String;
  def hashtags: List[String];
  def starNames: List[String];

  def date: String = TimeLib.modifiedJulianDayToStringJSTDate(time);
}

case class LegacyTweetContent(time: Double, message: String) extends  TweetContent {
  def hashtags: List[String] = Nil;
  def starNames: List[String] = Nil;
}

case class DateTweets(otherTweets: List[TweetContent], sunsetTweets: List[OnSunsetTweetContent],
  nightTweets: List[TweetContent]) {
  def isEmpty: Boolean = otherTweets.isEmpty && sunsetTweets.isEmpty && nightTweets.isEmpty;
  def isEmpty2: Boolean = sunsetTweets.isEmpty && nightTweets.isEmpty;
  def added(tc: TweetContent): DateTweets = {
    tc match {
      case tc: OnSunsetTweetContent =>
        this.copy(sunsetTweets = tc :: this.sunsetTweets);
      case _ =>
        if (isNightTime0(tc.time)) {
          this.copy(nightTweets = tc :: this.nightTweets);
        } else {
          this.copy(otherTweets = tc :: this.otherTweets);
        }
    }
  }
  def tweets: List[TweetContent] = {
    import Ordering.Double.IeeeOrdering;
    ((if (sunsetTweets.isEmpty) {
      Nil;
    } else if (sunsetTweets.tail.isEmpty) {
      sunsetTweets.head :: Nil;
    } else {
      MultiSunsetTweetContent(sunsetTweets.head.day, sunsetTweets.reverse) :: Nil;
    }) ::: otherTweets ::: nightTweets).sortBy(_.time);
  }
}

var _tweets: Map[String, DateTweets] = (0 until period).map { day =>
  val date = TimeLib.modifiedJulianDayToStringJSTDate(startTime + day);
  (date, DateTweets(Nil, Nil, Nil));
}.toMap;

def getTweets(time: Double): DateTweets = {
  val date = TimeLib.modifiedJulianDayToStringJSTDate(time);
  _tweets(date);
}

def putTweet(tc: TweetContent): Unit = {
  val date = TimeLib.modifiedJulianDayToStringJSTDate(tc.time);
  if (_tweets.contains(date)) {
    val tw = _tweets(date);
    _tweets = _tweets.updated(date, tw.added(tc));
  }
}

def putTweet(time: Double, msg: String): Unit = {
  putTweet(LegacyTweetContent(time, msg));
}

//==============================================================================
// 24節気
//==============================================================================

{
  val sunPhaseTerms = MathLib.findCyclicPhaseListContinuous(24, startTime, endTime, 3.0, 24) { time =>
    val utc = time;
    val tdb = TimeLib.mjdutcToTdb(utc);
    val sun = jplData.calcSunFromEarth(tdb);
    val bpnMatrix = Bpn.icrsToTrueEclipticMatrix(tdb);
    val sun2 = VectorLib.multiplyMV(bpnMatrix, sun);
    val sunLng = VectorLib.xyzToLng(sun2);
    sunLng;
  }
  val termStrs = IndexedSeq(
    "春分", "清明", "穀雨", "立夏", "小満", "芒種", "夏至", "小暑", "大暑", "立秋", "処暑", "白露",
    "秋分", "寒露", "霜降", "立冬", "小雪", "大雪", "冬至", "小寒", "大寒", "立春", "雨水", "啓蟄",
  );
  sunPhaseTerms.foreach { case (time, term) =>
    if (term % 6 == 0) {
      putTweet(TimeLib.floor(time, 24) + 1.0 / (24 * 4), "%s。太陽の黄経が%d°です".format(termStrs(term), term * 15));
    } else {
      putTweet(TimeLib.floor(time, 24) + 1.0 / (24 * 4), "二十四節気の%s。太陽の黄経が%d°です".format(termStrs(term), term * 15));
    }
  }
}

//==============================================================================
// 日没時の西の空
//==============================================================================

sealed trait OnSunsetTweetContent extends TweetContent {
  def day: Int;
  def message: String;
  def message2: String;
  def message3: String;
  def hashtags: List[String];
  def starNames: List[String];

  def time: Double = sunsetTimes(day);
}

case class MultiSunsetTweetContent(day: Int, tc: List[OnSunsetTweetContent]) extends TweetContent {
  def time: Double = sunsetTimes(day);
  def message: String = tc.head.message2 + "。" + tc.tail.map(_.message3).mkString("。");
  def hashtags: List[String] = tc.flatMap(_.hashtags);
  def starNames: List[String] = tc.flatMap(_.starNames);
}

case class SunsetMoonTweetContent(day: Int, azi: Double, alt: Double) extends OnSunsetTweetContent {
  def azi360: Int = (azi * PI57 + 0.5).toInt;
  def alt360: Int = (alt * PI57 + 0.5).toInt;
  def moonPhase: Double = calcMoonPhase(sunsetTimes(day));
  def message: String = "新月後の細い月は、月相 %.1f/28 で、日没時に西の空高度約%d°".format(moonPhase, alt360);
  def message2: String = "新月後の細い月は、月相 %.1f/28 で、日没時に西の空高度約%d°にいます".format(moonPhase, alt360);
  def message3: String = "新月後の細い月は、月相 %.1f/28 で、約%d°にいます".format(moonPhase, alt360);
  def hashtags: List[String] = Nil;
  def starNames: List[String] = List("月");
}
case class SunsetPlanetTweetContent(day: Int, planetName: String,
  azi: Double, alt: Double, isIncreasing: Boolean, isMax: Boolean) extends OnSunsetTweetContent {
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
  def message3: String = if (isMax) {
    "%sは日没時最大高度で西の空高度約%d°です".format(planetName, alt360);
  } else if (isIncreasing) {
    "%sは日没時の高度を徐々に上げ、約%d°にいます".format(planetName, alt360);
  } else {
    "%sは日没時の高度を徐々に下げ、約%d°にいます".format(planetName, alt360);
  }
  def hashtags: List[String] = List(planetName);
  def starNames: List[String] = List(planetName);
}
case class SunsetStarTweetContent(day: Int, starName: String,
  azi: Double, alt: Double) extends OnSunsetTweetContent {
  def alt360: Int = (alt * PI57 + 0.5).toInt;
  def message: String = "%sは西の空高度約%d°にいます".format(starName, alt360);
  def message2: String = "%sは西の空高度約%d°にいます".format(starName, alt360);
  def message3: String = "%sは約%d°にいます".format(starName, alt360);
  def hashtags: List[String] = List(starName);
  def starNames: List[String] = List(starName);
}

case class SunsetTweetContent(day: Int, flag: Int) extends OnSunsetTweetContent {
  def message: String = {
    if (flag == 1) {
      "日没はこのころが最も遅く、%sごろです".format(TimeLib.modifiedJulianDayToStringJSTNaturalTime(time));
    } else if (flag == 3) {
      "日没はこのころが最も早く、%sごろです".format(TimeLib.modifiedJulianDayToStringJSTNaturalTime(time));
    } else {
      if (day >= 7) {
        val timePrev = sunsetTimes(day - 7);
        val d = Math.round((time - (timePrev + 7.0)) * (24 * 60));
        if (d > 0) {
          "日没は%sごろです。この1週間で約%d分遅くなっています".format(TimeLib.modifiedJulianDayToStringJSTNaturalTime(time), d);
        } else if (d < 0) {
          "日没は%sごろです。この1週間で約%d分早くなっています".format(TimeLib.modifiedJulianDayToStringJSTNaturalTime(time), -d);
        } else {
          "日没は%sごろ".format(TimeLib.modifiedJulianDayToStringJSTNaturalTime(time));
        }
      } else {
        "日没は%sごろ".format(TimeLib.modifiedJulianDayToStringJSTNaturalTime(time));
      }
    }
  }
  def message2: String = {
    if (flag == 1) {
      "日没はこのころが最も遅く、%sごろです".format(TimeLib.modifiedJulianDayToStringJSTNaturalTime(time));
    } else if (flag == 3) {
      "日没はこのころが最も早く、%sごろです".format(TimeLib.modifiedJulianDayToStringJSTNaturalTime(time));
    } else {
      if (day >= 7) {
        val timePrev = sunsetTimes(day - 7);
        val d = Math.round((time - (timePrev + 7.0)) * (24 * 60));
        if (d > 0) {
          "日没は%sごろで、この1週間で約%d分遅くなっています".format(TimeLib.modifiedJulianDayToStringJSTNaturalTime(time), d);
        } else if (d < 0) {
          "日没は%sごろで、この1週間で約%d分早くなっています".format(TimeLib.modifiedJulianDayToStringJSTNaturalTime(time), -d);
        } else {
          "日没は%sごろです".format(TimeLib.modifiedJulianDayToStringJSTNaturalTime(time));
        }
      } else {
        "日没は%sごろです".format(TimeLib.modifiedJulianDayToStringJSTNaturalTime(time));
      }
    }
  }
  def message3: String = message2;
  def hashtags: List[String] = Nil;
  def starNames: List[String] = Nil;
}

// 日没時最大高度
{
  val planets = IndexedSeq(("金星", JplData.Venus), ("水星", JplData.Mercury));
  planets.foreach { p =>
    MathLib.findMaxMinListDiscrete(0, period, 15) { day =>
      val (azi, alt) = calcPlanetOnSunsetTime(day, p._2);
      alt;
    }.foreach { case (day, flag) =>
      if (flag > 0) {
        val (azi, alt) = calcPlanetOnSunsetTime(day, p._2);
        putTweet(SunsetPlanetTweetContent(day, p._1, azi, alt, true, true));
      }
    }
  }
}

// 新月直後の月
{
  val aziThres0 = 200 / PI57;
  val aziThres1 = 315 / PI57;
  val altThres0 = 10 / PI57;

  (0 until moonPhaseTerms.size).foreach { i =>
    val (moonPhaseTime, term) = moonPhaseTerms(i);
    if (term == 0 && moonPhaseTime + 4 < endTime) {
      val d = (0 until 4).indexWhere { d =>
        val day = (moonPhaseTime + d - startTime).toInt;
        val (azi, alt) = calcPlanetOnSunsetTime(day, JplData.Moon);
        alt >= altThres0;
      }
      val day = (moonPhaseTime + d - startTime).toInt;
      val (azi, alt) = calcPlanetOnSunsetTime(day, JplData.Moon);
      putTweet(SunsetMoonTweetContent(day, azi, alt));
    }
  }
}

// 水曜・金曜
{
  val aziThres0 = 200 / PI57;
  val aziThres1 = 315 / PI57;
  val altThres0 = 10 / PI57;

  val planets = IndexedSeq(("金星", JplData.Venus, 5), ("水星", JplData.Mercury, 3));
  (1 until period).foreach { day =>
    val wday = TimeLib.wday(startTime + day);
    val sunsetTweets = getTweets(startTime + day).sunsetTweets;
    planets.foreach { p =>
      if (wday == p._3 && !sunsetTweets.exists(_.starNames.contains(p._1))) {
        val (azi, alt) = calcPlanetOnSunsetTime(day, p._2);
        if (azi >= aziThres0 && azi <= aziThres1 && alt >= altThres0) {
          val (_, prevAlt) = calcPlanetOnSunsetTime(day - 1, p._2);
          val isIncreasing = (alt >= prevAlt);
          putTweet(SunsetPlanetTweetContent(day, p._1, azi, alt, isIncreasing, false));
        }
      }
    }
  }
}

// 日没ツイートがある場合に他の天体のツイートも追加
{
  val aziThres0 = 200 / PI57;
  val aziThres1 = 315 / PI57;
  val altThres0 = 10 / PI57;

  val planets = IndexedSeq(("金星", JplData.Venus), ("水星", JplData.Mercury));
  val planets2 = IndexedSeq(("火星", JplData.Mars), ("木星", JplData.Jupiter), ("土星", JplData.Saturn));
  (1 until period).foreach { day =>
    val sunsetTweets = getTweets(startTime + day).sunsetTweets;
    if (sunsetTweets.nonEmpty) {
      planets.foreach { p =>
        if (!sunsetTweets.exists(_.starNames.contains(p._1))) {
          val (azi, alt) = calcPlanetOnSunsetTime(day, p._2);
          if (azi >= aziThres0 && azi <= aziThres1 && alt >= altThres0) {
            val (_, prevAlt) = calcPlanetOnSunsetTime(day - 1, p._2);
            val isIncreasing = (alt >= prevAlt);
            putTweet(SunsetPlanetTweetContent(day, p._1, azi, alt, isIncreasing, false));
          }
        }
      };
      {
        val p = ("月", JplData.Moon);
        if (!sunsetTweets.exists(_.starNames.contains(p._1))) {
          val (azi, alt) = calcPlanetOnSunsetTime(day, p._2);
          if (azi >= aziThres0 && azi <= aziThres1 && alt >= altThres0) {
            putTweet(SunsetMoonTweetContent(day, azi, alt));
          }
        }
      }
      planets2.foreach { p =>
          val (azi, alt) = calcPlanetOnSunsetTime(day, p._2);
          if (azi >= aziThres0 && azi <= aziThres1 && alt >= altThres0) {
            putTweet(SunsetStarTweetContent(day, p._1, azi, alt));
          }
      };
    }
  }
}

// 日没時間
{
  MathLib.getMaxMinUpDownFlagListDiscrete(0, period, 90) { day =>
    sunsetTimes(day) - startTime - day;
  }.zipWithIndex.foreach { case (flag, day) =>
    val wday = TimeLib.wday(startTime + day);
    if (flag == 1 || flag == 3 || wday == 0) {
      putTweet(SunsetTweetContent(day, flag));
    }
  }
}

//==============================================================================
// 月相
//==============================================================================

{
  case class MoonPhaseTermTweetContent(rawTime: Double, term: Int, distanceFlag: Int,
    cons: Option[(String, String, List[String])]) extends TweetContent {
    def time: Double = TimeLib.floor(rawTime, 24) + 1.0 / (24 * 4);
    def message: String = {
      val msg = if (term == 4 && distanceFlag < 0) {
        "満月\uD83C\uDF15。月が地球に近く、もっとも大きい満月です。月相 14/28";
      } else if (term == 4 && distanceFlag > 0) {
        "満月\uD83C\uDF15。月が地球から遠く、もっとも小さい満月です。月相 14/28";
      } else {
        termStrs(term);
      }
      cons match {
        case Some((conscomment, cons, hashtags)) => "%s%s。%sにいます".format(conscomment, msg, cons);
        case None => msg;
      }
    }
    def hashtags: List[String] = {
      cons match {
        case Some((conscomment, cons, hashtags)) => hashtags;
        case None => Nil;
      }
    }
    def starNames: List[String] = List("月");
    private[this] val termStrs = IndexedSeq(
      "新月",
      "月相 3.5/28。新月と上弦の中間です\uD83C\uDF12",
      "上弦の月\uD83C\uDF13。月相 7/28",
      "月相 10.5/28。上弦と満月の中間です\uD83C\uDF14",
      "満月\uD83C\uDF15。月相 14/28",
      "月相 17.5/28。満月と下弦の中間です\uD83C\uDF16",
      "下弦の月\uD83C\uDF17。月相 21/28",
      "月相 24.5/28。下弦と新月の中間です\uD83C\uDF18",
    );
  }

  val altThres0 = 10 / PI57;

  def calcMoonConstellation(time: Double): Option[(String, String, List[String])] = {
    if (isNightTime2(time)) {
      val ut1 = time; // 近似的
      val tdb = TimeLib.mjdutcToTdb(time);
      val xyz = jplData.calcMoonFromEarth(tdb);
      val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
      val xyz2 = VectorLib.multiplyMV(bpnMatrix, xyz);
      val (azi, alt) = hcs.trueEquatorialXyzToAziAlt(xyz2, ut1);
      if (alt >= altThres0) {
        val (conscomment, cons, hashtags) = Constellations.icrsToConstellation(xyz);
        Some((conscomment, cons, hashtags));
      } else {
        None;
      }
    } else {
      None;
    }
  }

  {
    val fullMoons = moonPhaseTerms.filter(_._2 == 4).map(_._1);
    val fullMoonsDistanceMaxMinUpDownFlags = MathLib.getMaxMinUpDownFlagListDiscrete(0, fullMoons.size, 2) { idx =>
      val utc = fullMoons(idx);
      val tdb = TimeLib.mjdutcToTdb(utc);
      val moon = jplData.calcMoonFromEarth(tdb);
      VectorLib.distance(moon);
    }
    fullMoonsDistanceMaxMinUpDownFlags.zipWithIndex.map { case (flag, idx) =>
      val time = fullMoons(idx);
      if (flag == 1) {
        MoonPhaseTermTweetContent(time, 4, +1, calcMoonConstellation(time));
      } else if (flag == 3) {
        MoonPhaseTermTweetContent(time, 4, -1, calcMoonConstellation(time));
      } else {
        MoonPhaseTermTweetContent(time, 4, 0, calcMoonConstellation(time));
      }
    }.foreach(putTweet);
  }
  moonPhaseTerms.filter(_._2 != 4).map { case (time, term) =>
    MoonPhaseTermTweetContent(time, term, 0, calcMoonConstellation(time));
  }.foreach(putTweet);
}

//==============================================================================
// 近日点・遠日点
//==============================================================================

MathLib.findMaxMinListContinuous(startTime, endTime, 30, 24) { time =>
  val utc = time;
  val tdb = TimeLib.mjdutcToTdb(utc);
  val sun = jplData.calcSunFromEarth(tdb);
  VectorLib.distance(sun);
}.foreach { case (time, flag) =>
  val s = if (flag < 0) "近日点" else "遠日点";
  putTweet(TimeLib.round(time, 24) - 1.0 / (24 * 4), "地球が%s通過".format(s));
}

//==============================================================================
// 惑星の天象
//==============================================================================

case class PlanetAstronomyTweetContent(time: Double, message: String, planetName: String) extends TweetContent {
  def hashtags: List[String] = List(planetName);
  def starNames: List[String] = List(planetName);
}

{
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

  {
    val (planetName, targetPlanet) = ("水星", JplData.Mercury);
    MathLib.findMaxMinCrossingListContinuous(startTime, endTime, 10.0, 24) { time =>
      calcInnerPlanetLngEc(time, targetPlanet);
    }.map { case (time, term) =>
      if (term == 0) {
        putTweet(PlanetAstronomyTweetContent(TimeLib.floor(time, 24) + 1.0 / (24 * 4),
          "水星が外合(黄経基準)", planetName));
      } else if (term == 1) {
        putTweet(PlanetAstronomyTweetContent(TimeLib.round(time, 24) - 1.0 / (24 * 4),
          "水星が東方最大離角(黄経基準)\uD83C\uDF13", planetName));
      } else if (term == 2) {
        putTweet(PlanetAstronomyTweetContent(TimeLib.floor(time, 24) + 1.0 / (24 * 4),
          "水星が内合(黄経基準)", planetName));
      } else {
        putTweet(PlanetAstronomyTweetContent(TimeLib.round(time, 24) - 1.0 / (24 * 4),
          "水星が西方最大離角(黄経基準)\uD83C\uDF17", planetName));
      }
    }
  }

  {
    val (planetName, targetPlanet) = ("金星", JplData.Venus);
    MathLib.findMaxMinCrossingListContinuous(startTime, endTime, 10.0, 24) { time =>
      calcInnerPlanetLngEc(time, targetPlanet);
    }.map { case (time, term) =>
      if (term == 0) {
        putTweet(PlanetAstronomyTweetContent(TimeLib.floor(time, 24) + 1.0 / (24 * 4),
          "金星が外合(黄経基準)。数か月後に夕方の西の空に現れます", planetName));
      } else if (term == 1) {
        putTweet(PlanetAstronomyTweetContent(TimeLib.round(time, 24) - 1.0 / (24 * 4),
          "金星が東方最大離角(黄経基準)\uD83C\uDF13。宵の明星として夕方に西の空にいます", planetName));
      } else if (term == 2) {
        putTweet(PlanetAstronomyTweetContent(TimeLib.floor(time, 24) + 1.0 / (24 * 4),
          "金星が内合(黄経基準)。数週間後に明け方の東の空に現れます", planetName));
      } else {
        putTweet(PlanetAstronomyTweetContent(TimeLib.round(time, 24) - 1.0 / (24 * 4),
          "金星が西方最大離角(黄経基準)\uD83C\uDF17。明けの明星として明け方に東の空にいます", planetName));
      }
    }
  }
}

{
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

  val planetPhases1: List[(Double, String, String, Boolean, Option[Array[Double]])] = List(
    ("火星", JplData.Mars),
    ("木星", JplData.Jupiter),
    ("土星", JplData.Saturn),
  ).flatMap { t =>
    val (planetName, targetPlanet) = t;
    MathLib.findCyclicPhaseListContinuous(4, startTime, endTime, 30, 24) { time =>
      PI2 - calcOuterPlanetLngEq(time, targetPlanet);
    }.map { case (time, term) =>
      val utc = time;
      val tdb = TimeLib.mjdutcToTdb(utc);
      val xyz = jplData.calcPlanetFromEarth(tdb, targetPlanet);
      if (term == 0) {
        (time, planetName, "合(赤経基準)", false, None);
      } else if (term == 1) {
        (time, planetName, "西矩(赤経基準)", false, Some(xyz));
      } else if (term == 2) {
        (time, planetName, "衝(赤経基準)", false, Some(xyz));
      } else {
        (time, planetName, "東矩(赤経基準)", false, Some(xyz));
      }
    }
  }

  planetPhases1.foreach { case (time, planetName, content, timeFlag, xyzOpt) =>
    val time2 = if (timeFlag) {
      TimeLib.round(time, 24) - 1.0 / (24 * 4);
    } else {
      TimeLib.floor(time, 24) + 1.0 / (24 * 4);
    }
    xyzOpt match {
    case None =>
      putTweet(time2, "%sが%s #%s".format(planetName, content, planetName));
    case Some(xyz) =>
      val (conscomment, cons, hashtags) = Constellations.icrsToConstellation(xyz);
      putTweet(time2, "%s%sが%s。%sにいます #%s".format(conscomment, planetName, content, cons, planetName) +
        hashtags.map(" #" + _).mkString);
    }
  }
}

//==============================================================================
// 月・惑星と恒星の会合
//==============================================================================

case class CloseStarsTweetContent(rawTime: Double, stepCountPerDay: Int, slowStarName: String, fastStarName: String,
  distance: Double, hashtags: List[String]) extends TweetContent {
  def time: Double = TimeLib.round(rawTime, stepCountPerDay);
  def message: String = "%sが%s%s".format(fastStarName, slowStarName, distanceStr);
  def distanceStr: String = {
    val distance360 = distance * PI57;
    if (distance360 < 1.0) {
      "に接近 (1°未満)";
    } else if (distance360 < 2.0) {
      "に接近 (2°未満)";
    } else if (distance360 < 3.0) {
      "に接近 (3°未満)";
    } else {
      "の近くにいます";
    }
  }
  def starNames: List[String] = List(slowStarName, fastStarName);
}

{
  val altThres0 = 10 / PI57;
  val distanceThres = 10.0 / PI57;

  def calcClosest1(slowStarName: String, fastStarName: String, hashtags: List[String],
    slowStarXyzFunc: Double => Array[Double],
    fastStarXyzFunc: Double => Array[Double]): Unit = {
    MathLib.findMaxMinListContinuous(startTime, endTime, 10.0, 24 * 6) { time =>
      val utc = time;
      val tdb = TimeLib.mjdutcToTdb(utc);
      val xyz_s = slowStarXyzFunc(tdb);
      val xyz_f = fastStarXyzFunc(tdb);
      VectorLib.angularDistance1(xyz_s, xyz_f);
    }.foreach { case (time, flag) =>
      if (flag > 0 && isNightTime2(time)) {
        val ut1 = time; // 近似的
        val tdb = TimeLib.mjdutcToTdb(time);
        val xyz_f = fastStarXyzFunc(tdb);
        val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
        val xyz_f2 = VectorLib.multiplyMV(bpnMatrix, xyz_f);
        val (azi, alt) = hcs.trueEquatorialXyzToAziAlt(xyz_f2, ut1);
        if (alt >= altThres0) {
          val xyz_s = slowStarXyzFunc(tdb);
          val distance = VectorLib.angularDistance(xyz_s, xyz_f);
          if (distance < distanceThres) {
            putTweet(CloseStarsTweetContent(time, 24 * 6, slowStarName, fastStarName, distance, hashtags));
          }
        }
      }
    }
  }
  def calcClosest2(slowStarName: String, fastStarName: String, hashtags: List[String],
    slowStarXyzFunc: Double => Array[Double],
    fastStarXyzFunc: Double => Array[Double]): Unit = {
    MathLib.findMaxMinListContinuous(startTime, endTime, 10.0, 24) { time =>
      val utc = time;
      val tdb = TimeLib.mjdutcToTdb(utc);
      val xyz_s = slowStarXyzFunc(tdb);
      val xyz_f = fastStarXyzFunc(tdb);
      VectorLib.angularDistance1(xyz_s, xyz_f);
    }.foreach { case (time, flag) =>
      if (flag > 0) {
        val ut1 = time; // 近似的
        val tdb = TimeLib.mjdutcToTdb(time);
        val xyz_f = fastStarXyzFunc(tdb);
        val xyz_s = slowStarXyzFunc(tdb);
        val distance = VectorLib.angularDistance(xyz_s, xyz_f);
        if (distance < distanceThres) {
          getConjunctionTweetTime(time, xyz_f) match {
            case Some(postTime) =>
              putTweet(CloseStarsTweetContent(postTime, 24, slowStarName, fastStarName, distance, hashtags));
            case None =>
              // nop
          }
        }
      }
    }
  }

  val stars0: IndexedSeq[(String, Array[Double], List[String])] = IndexedSeq(
    ("03h47m", "+24°06′", "おうし座すばる", List("プレアデス星団")),
    ("04h36m", "+16°31′", "おうし座アルデバラン", Nil),
    ("07h45m", "+28°02′", "ふたご座ポルックス", Nil),
    ("10h08m", "+11°58′", "しし座レグルス", Nil),
    ("13h25m", "-11°09′", "おとめ座スピカ", Nil),
    ("16h29m", "-26°26′", "さそり座アンタレス", Nil),
  ).map { case (lngStr, latStr, name, hashtags) =>
    val lng = (lngStr.substring(0, 2).toInt.toDouble + lngStr.substring(3, 5).toInt.toDouble / 60) / 24 * PI2;
    val lat = (latStr.substring(0, 3).toInt.toDouble + latStr.substring(4, 6).toInt.toDouble / 60) / 360 * PI2;
    val c = Math.cos(lat);
    val x = c * Math.cos(lng);
    val y = c * Math.sin(lng);
    val z = Math.sin(lat);
    (name, Array(x, y, z), hashtags);
  }
  val stars1 = IndexedSeq(
    ("水星", JplData.Mercury),
    ("金星", JplData.Venus),
    ("火星", JplData.Mars),
    ("木星", JplData.Jupiter),
    ("土星", JplData.Saturn),
  );
  val stars2 = IndexedSeq(
    ("月", JplData.Moon),
  );
  stars2.foreach { star2 =>
    stars0.foreach { star0 =>
      calcClosest1(star0._1, star2._1, star0._3,
        { tdb: Double =>
          star0._2;
        },
        { tdb: Double =>
          jplData.calcPlanetFromEarth(tdb, star2._2);
        });
    }
  }
  stars2.foreach { star2 =>
    stars1.foreach { star1 =>
      calcClosest1(star1._1, star2._1, List(star1._1),
        { tdb: Double =>
          jplData.calcPlanetFromEarth(tdb, star1._2);
        },
        { tdb: Double =>
          jplData.calcPlanetFromEarth(tdb, star2._2);
        });
    }
  }
  stars1.foreach { star1 =>
    stars0.foreach { star0 =>
      calcClosest2(star0._1, star1._1, star1._1 :: star0._3,
        { tdb: Double =>
          star0._2;
        },
        { tdb: Double =>
          jplData.calcPlanetFromEarth(tdb, star1._2);
        });
    }
  }
  (0 until stars1.size).foreach { i =>
    ((i + 1) until stars1.size).foreach { j =>
      val star1 = stars1(i);
      val star0 = stars1(j);
      calcClosest2(star0._1, star1._1, List(star1._1, star0._1),
        { tdb: Double =>
          jplData.calcPlanetFromEarth(tdb, star0._2);
        },
        { tdb: Double =>
          jplData.calcPlanetFromEarth(tdb, star1._2);
        });
    }
  }
}

//==============================================================================
// 21時・23時の惑星
//==============================================================================

{
  var moonCons: List[(Double, Array[Double])] = Nil;
  val altThres = 10 / PI57;
  (0 until period).foreach { d =>
    val time = startTime + d + 21.0 / 24.0 - 1.0 / (24 * 6);
    if (!getTweets(time).nightTweets.flatMap(_.starNames).contains("月")) {
      val (xyz, azi, alt) = calcPlanetXyzAziAlt(time, JplData.Moon);
      if (alt >= altThres) {
        val (conscomment, cons, hashtags) = Constellations.icrsToConstellation(xyz);
        val moonPhase = calcMoonPhase(time);
        putTweet(time, "%s月は%sにいます。月相 %.1f/28".format(conscomment, cons, moonPhase) +
          hashtags.map(" #" + _).mkString);
      } else {
        val time = startTime + d + 23.0 / 24.0;
        val (xyz, azi, alt) = calcPlanetXyzAziAlt(time, JplData.Moon);
        if (alt >= altThres) {
          val (conscomment, cons, hashtags) = Constellations.icrsToConstellation(xyz);
          val moonPhase = calcMoonPhase(time);
          putTweet(time, "%s月は%sにいます。月相 %.1f/28".format(conscomment, cons, moonPhase) +
            hashtags.map(" #" + _).mkString);
          }
      }
    }
  }
  moonCons;
}

val planetCons: List[(Double, String, Array[Double])] = {
  var planetCons: List[(Double, String, Array[Double])] = Nil;
  val altThres = 10 / PI57;
  (0 until period).foreach { d =>
    val time = startTime + d + 21.0 / 24.0;
    val wday = TimeLib.wday(time);
    val ut1 = time; // 近似的
    val tdb = TimeLib.mjdutcToTdb(time);
    val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
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
  planetCons;
}

planetCons.foreach { case (time, planetName, xyz) =>
  val (conscomment, cons, hashtags) = Constellations.icrsToConstellation(xyz);
  putTweet(TimeLib.floor(time, 24) - 1.0 / (24 * 4), "%s%sは%sにいます #%s".format(conscomment, planetName, cons, planetName) +
    hashtags.map(" #" + _).mkString);
}

//==============================================================================

val manualTweet = """
2021-11-05T12:00 天王星が衝
2021-09-15T12:00 海王星が衝
2021-07-18T12:00 冥王星が衝

2021-11-28T12:00 準惑星ケレスが衝
2021-09-09T12:00 小惑星パラスが衝
2021-06-08T12:00 小惑星ジュノーが衝
2021-03-08T12:00 小惑星ベスタが衝

2021-04-28T12:00 準惑星ハウメアが衝
2021-04-09T12:00 準惑星マケマケが衝
2021-10-22T12:00 準惑星エリスが衝
""";

manualTweet.split("\n").foreach { line =>
  val cols = line.split(" ", 2);
  if (cols.length == 2) {
    val time = TimeLib.stringToModifiedJulianDay(cols(0) + ":00+09:00");
    val msg = cols(1);
     putTweet(time, msg);
  }
}


//==============================================================================
// 星座
//==============================================================================

def tweetConstellations(data: IndexedSeq[(String, String)], span: Int, startDay: Int): Unit = {
  val altHor = -0.90 / PI57;
  val hashtag = " #星空";
  var nextDay: Int = startDay;
  (0 until period).foreach { day =>
    val date = TimeLib.modifiedJulianDayToStringJSTDate(startTime + day);
    if (day >= nextDay && getTweets(startTime + day).isEmpty2) {
      val time = startTime + day + 21.0 / 24.0;
      val sid = hcs.siderealTime(time);
      val cons = Constellations.siderealTimeToConstellation(sid, data);
      if (cons != "") {
        val msg = if (calcPlanetXyzAziAlt(time, JplData.Moon)._3 >= altHor) {
          cons;
        } else {
          cons + "。月明かりなし";
        }
        putTweet(time, msg + hashtag);
      }
      nextDay = day + span;
    }
  }
}
//tweetConstellations(Constellations.ecliptical, 14, 59);
tweetConstellations(Constellations.season, 7, 36); // PERIOD
//tweetConstellations(Constellations.northern, 14, 39);

{
  import Ordering.Double.IeeeOrdering;
  case class StarTweetContent(time: Double, message: String, hashtags2: List[String]) extends TweetContent {
    def hashtags: List[String] = "星空" :: hashtags2;
    def starNames: List[String] = Nil;
  }
  val stars = IndexedSeq(
    // J2000

    // 冬の星座
    ("03h08m", "ペルセウス座アルゴル(2等星)が南中。天頂付近です", Nil),
    ("03h24m", "ペルセウス座α(2等星)が南中。北の空高くにいます", Nil),
    ("03h47m", "おうし座すばるが南中。天頂付近です", List("プレアデス星団")),
    ("04h30m", "この時期に南の空高くにいる黄道十二星座はおうし座です。西から東に順番に、うお座、おひつじ座、おうし座、ふたご座、かに座、の順に並んでいます", Nil),
    ("04h36m", "おうし座アルデバランが南中。南の空高くにいます。冬のダイヤモンドを構成する6個の星の1つです", Nil),
    ("05h26m", "おうし座βが南中。天頂付近です。おうし座の角であり、ぎょしゃ座の五角形の1つでもあります", Nil),
    ("04h50m", "オリオン座π3(3等星)が南中。オリオン座の西にある星です", List("オリオン座")),
    ("05h15m", "オリオン座リゲルが南中。オリオン座の南西(右下)の星で、冬のダイヤモンドを構成する6個の星の1つです", List("オリオン座")),
    ("05h17m", "ぎょしゃ座カペラが南中。天頂付近です。冬のダイヤモンドを構成する6個の星の1つです", Nil),
    ("05h36m", "オリオン座εが南中。オリオン座の三つ星の中央の星です。その南(下)にはオリオン座大星雲がいます", List("オリオン座")),
    ("05h55m", "オリオン座ベテルギウスが南中。オリオン座の北東(左上)の星で、冬のダイヤモンドの中心の星です", List("オリオン座")),
    ("05h33m", "うさぎ座α(3等星)が南中。うさぎ座はオリオン座の南(下)にいる目立たない星座です", Nil),
    ("06h45m", "おおいぬ座シリウスが南中。太陽を除いてもっとも明るい恒星です。冬のダイヤモンドを構成する6個の星の1つです", Nil),
    ("06h29m", "いっかくじゅう座βが南中。いっかくじゅう座はオリオン座とこいぬ座とおおいぬ座に囲まれた領域の目立たない星座です", Nil),
    ("07h39m", "こいぬ座プロキオンが南中。冬のダイヤモンドを構成する6個の星の1つです", Nil),
    ("07h00m", "この時期に南の空高くにいる黄道十二星座はふたご座です。西から東に順番に、おひつじ座、おうし座、ふたご座、かに座、しし座の順に並んでいます", Nil),
    ("06h38m", "ふたご座γが南中。南の空高くにいます。ふたご座の東側の子の足元の星です", Nil),
    ("07h35m", "ふたご座の西側の明るい星カストルが南中。天頂付近です。東側のポルックスよりよりやや暗いです", Nil),
    ("07h45m", "ふたご座の東側の明るい星ポルックスが南中。天頂付近です。西側のカストルよりも明るく、冬のダイヤモンドを構成する6個の星の1つです", Nil),

    ("08h08m", "とも座ρ(3等星)が南中。とも座はおおいぬ座の南東(左下)にあり、日本からは南の空低くにしか見られません", Nil),

    // 春の北の空の目立たない星座
    ("08h17m", "かに座β(4等星)が南中。かに座の南西(右下)の星です", Nil),
    ("08h40m", "この時期に南の空高くにいる黄道十二星座はかに座です。西から東に順番に、おうし座、ふたご座、かに座、しし座、おとめ座の順に並んでいます", Nil),
    ("09h21m", "やまねこ座α(3等星)が南中。天頂付近です。やまねこ座は、ぎょしゃ座・ふたご座・かに座とおおぐま座の隙間にある目立たない星座で、見つけるにはヤマネコのような暗視能力が必要だとも", Nil),
    ("10h53m", "こじし座で一番明るい星「46番星」(4等星)が南中。天頂付近です。こじし座は、しし座の頭とおおぐま座の隙間にある目立たない星座です。小さな獅子(しし)の意味", Nil),
    ("12h23m", "かみのけ座の散開星団が南中。天頂付近です。りょうけん座コル・カロリ(3等星)としし座デネブの中間にあります。かみのけ座は、この散開星団を中心とした目立たない星座で、銀河北極点に位置します。かみのけ座銀河団もあります", Nil),
    ("12h56m", "りょうけん座コル・カロリ(3等星)が南中。天頂付近です。春の大三角形であるうしかい座アークトゥルス、おとめ座スピカ、しし座デネブにコル・カロリを加えると、春のダイヤモンドです", Nil),

    // うみへび座付近
    ("08h55m", "うみへび座の頭にある星ζ(3等星)が南中。うみへび座の頭はかに座の南(下)にあります。うみへび座は、全天でもっとも大きな星座で、東西に長く、尾はてんびん座付近です", Nil),
    ("09h28m", "うみへび座α(2等星)が南中。しし座レグルスの南西(右下)です。うみへび座は、全天でもっとも大きな星座で、東西に長く、頭はかに座付近、尾はてんびん座付近です", Nil),
    ("13h19m", "うみへび座γ(3等星)が南中。おとめ座スピカの南(下)です。うみへび座は、全天でもっとも大きな星座で、東西に長く、頭はかに座付近、尾はてんびん座付近です", Nil),
    ("10h24m", "ろくぶんぎ座が南中。ろくぶんぎ座は、しし座レグルスの南(下)、うみへび座αの東(左)にある目立たない星座です", Nil),
    ("11h19m", "コップ座が南中。コップ座は、しし座の後ろ足の南(下)、からす座の西(右)にある目立たない星座です。コップというよりは杯です", Nil),
    ("12h16m", "からす座が南中。からす座は、おとめ座スピカの南西(右下)にある目立たない星座です。北斗七星、うしかい座アークトゥルス、おとめ座スピカ、からす座で春の大曲線を構成します", Nil),

    // 春の大三角形
    ("10h08m", "しし座の前足にあるレグルスが南中。東(左)にはしし座の尾のデネボラがいます", Nil),
    ("10h40m", "この時期に南の空にいる黄道十二星座はしし座です。西から東に順番に、ふたご座、かに座、しし座、おとめ座の順に並んでいます", Nil),
    ("11h49m", "しし座の尾にあるデネボラが南中。おとめ座スピカ、うしかい座アークトゥルスとともに春の大三角形を構成します。西(右)にはしし座の前足のレグルスがいます", Nil),
    ("13h00m", "この時期に南の空にいる黄道十二星座はおとめ座です。西から東に順番に、ふたご座、かに座、しし座、おとめ座、てんびん座、さそり座の順に並んでいます", Nil),
    ("13h25m", "おとめ座スピカが南中。しし座デネボラ、うしかい座アークトゥルスとともに春の大三角形を構成します", Nil),
    ("14h16m", "うしかい座アークトゥルスが南中。南の空高くにいます。しし座デネボラ、おとめ座スピカとともに春の大三角形を構成します。全天でシリウス、カノープスに次ぐ明るさの恒星です", Nil),

    // へびつかい座付近
    ("15h35m", "かんむり座αが南中。天頂付近です。かんむり座は、西のうしかい座と東のヘラクレス座に挟まれ、へび座の頭の北にあります。唯一の2等星αを中心に小さなはっきりとした半円形のかんむりを描いています", Nil),
    ("17h00m", "ヘラクレス座が南中。天頂付近です。ヘラクレス座は、大きな星座ですが、明るい星はありません。ヘラクレス座の西はかんむり座、南はへびつかい座、東はこと座、北はりゅう座です", Nil),
    ("15h44m", "へび座α(3等星)が南中。へび座でもっとも明るい星で、へびの頭にあります。夏の大三角形のベガとアルタイルを軸にして、デネブと線対称の位置です", Nil),
    ("17h35m", "へびつかい座α(2等星)が南中。へびつかい座でもっとも明るい星です。へびつかい座は、さそり座の北(上)、ヘラクレス座の南にある星座です。へびつかいの持っているへびの頭はかんむり座付近、尾はわし座付近です", Nil),
    ("18h21m", "へび座η(3等星)が南中。へび座の尾にある星です", Nil),
    // TODO たて座

    // 夏の大三角形
    ("18h37m", "こと座ベガが南中。天頂付近です。夏の大三角形の1つです。ベガの東はわし座アルタイルです", Nil),
    ("19h51m", "わし座アルタイルが南中。夏の大三角形であるベガ・アルタイル・デネブは空高くにいます", Nil),
    ("20h41m", "はくちょう座デネブが南中。天頂付近です。夏の大三角形の1つです。デネブの西はこと座ベガです", Nil),
    // TODO こぎつね、いるか、こうま

    // さそり座付近
    ("15h10m", "この時期に南の空にいる黄道十二星座はてんびん座です。西から東に順番に、しし座、おとめ座、てんびん座、さそり座、いて座の順に並んでいます", Nil),
    ("16h29m", "さそり座アンタレスが南中。アンタレスは火星のように赤い星です", Nil),
    ("16h50m", "この時期に南の空にいる黄道十二星座はさそり座です。西から東に順番に、おとめ座、てんびん座、さそり座、いて座、やぎ座の順に並んでいます", Nil),
    ("18h40m", "この時期に南の空にいる黄道十二星座はいて座です。西から東に順番に、てんびん座、さそり座、いて座、やぎ座、みずがめ座の順に並んでいます。いて座は南斗六星が目印です。銀河系の中心方向であり、天の川がもっとも濃い領域です", Nil),

    // 秋の黄道星座
    ("20h21m", "やぎ座βが南中。やぎ座の西(右)端の星で、夏の大三角形のベガからアルタイルへの直線を延ばしたところにいます", Nil),
    ("21h00m", "この時期に南の空にいる黄道十二星座はやぎ座です。西から東に順番に、いて座、やぎ座、みずがめ座、うお座の順に並んでいます", Nil),
    ("21h47m", "やぎ座δが南中。やぎ座の東(左)端の星です", Nil),
    ("22h48m", "みなみのうお座フォーマルハウトが南中", Nil),
    ("22h20m", "この時期に南の空にいる黄道十二星座はみずがめ座です。西から東に順番に、やぎ座、みずがめ座、うお座、おひつじ座、おうし座の順に並んでいます", Nil),

    // ペガスス
    ("23h04m", "ペガスス座の四角形の西(右)の辺が南中。南の空高くにいます", Nil),
    // TODO とかげ
    ("00h08m", "ペガスス座の四角形の東(左)の辺が南中。南の空高くにいます。ペガスス座の四角形の左上の星はアンドロメダ座α(2等星)です", Nil),
    ("01h10m", "アンドロメダ座β(2等星)が南中。天頂付近です。アンドロメダ座は西(右)から東(左)にα、β、γと2等星が並んでいます。βの近くにアンドロメダ銀河がいます。アンドロメダ銀河は40億年後に天の川銀河と衝突します", Nil),
    ("02h04m", "アンドロメダ座γ(2等星)が南中。天頂付近です。アンドロメダ座は西(右)から東(左)にα、β、γと2等星が並んでいます", Nil),

    // ペガススの東
    // TODO さんかく
    ("23h17m", "うお座γ(4等星)が南中。うお座の東側の魚の頭にいます。ペガススの四角形の南(下)です", Nil),
    ("00h50m", "この時期に南の空にいる黄道十二星座はうお座です。西から東に順番に、やぎ座、みずがめ座、うお座、おひつじ座、おうし座、ふたご座の順に並んでいます", Nil),
    ("01h31m", "うお座η(4等星)が南中。うお座で最も明るい星で、うお座の西側の魚の中ほどにいます", Nil),
    ("02h07m", "おひつじ座α(2等星)が南中。南の空高くです", Nil),
    ("02h40m", "この時期に南の空高くにいる黄道十二星座はおひつじ座です。西から東に順番に、みずがめ座、うお座、おひつじ座、おうし座、ふたご座、かに座の順に並んでいます", Nil),
    // TODO くじら、エリダヌス

    ("00h57m", "カシオペア座の中央の星が南中。北極星の上方向にいます", Nil),
    ("05h04m", "北斗七星の端(スプーンのすくう側)の星が北極星のちょうど右方向にいます。カシオペア座は北極星の左上方向にいます", Nil),
    ("06h57m", "カシオペア座が北極星のちょうど左方向にいます。北斗七星は北極星の右方向にいます", Nil),
    ("11h04m", "北斗七星の端(スプーンのすくう側)の星が南中。北極星の上方向にいます", Nil),
    ("13h48m", "北斗七星の柄の部分の端の星が南中。北の空高く、北極星の上方向にいます", Nil),
    ("14h51m", "こぐま座βが南中。北極星から上にこぐま座が伸びています", Nil),
    ("17h04m", "北斗七星の端(スプーンのすくう側)の星が北極星のちょうど左方向にいます。カシオペア座は北極星の右下方向にいます", Nil),
    ("17h57m", "りゅう座γが南中。りゅう座の頭にある星です。りゅう座の尾は北斗七星の近くです", Nil),
    ("18h57m", "カシオペア座が北極星のちょうど右方向にいます。北斗七星は北極星の左方向にいます", Nil),
    // TODO ケフェウス、きりん

  ).map { t =>
    ((t._1.substring(0, 2).toInt.toDouble + t._1.substring(3, 5).toInt.toDouble / 60) / 24 * PI2,
    t._2, t._3);
  }.sortBy(_._1);

  val altHor = -0.90 / PI57;

  var index: Int = -1;
  (52 until period).foreach { day => // PERIOD
    val date = TimeLib.modifiedJulianDayToStringJSTDate(startTime + day);
    if (index < 0) {
      val time = startTime + day + 21.0 / 24.0; // PERIOD
      val sid = hcs.siderealTime(time);
      index = stars.indexWhere(_._1 > sid);
      if (index < 0) {
        index = 0;
      }
    } else {
      val time = startTime + day + 21.0 / 24.0;
      val sid = hcs.siderealTime(time);
      val time0 = if (MathLib.circleAdd(sid, -stars(index)._1) >= 0) {
        time;
      } else if (getTweets(time).isEmpty2) {
        val sid2 = MathLib.circleAdd(sid, PI2 / 12);
        if (MathLib.circleAdd(sid2, -stars(index)._1) >= 0) {
          time + MathLib.circleAdd(stars(index)._1, -sid) / PI2;
        } else {
          0.0;
        }
      } else {
        0.0;
      }
      if (time0 > 0.0) {
        val msg = if (calcPlanetXyzAziAlt(time0, JplData.Moon)._3 >= altHor) {
          stars(index)._2;
        } else {
          stars(index)._2 + "。月明かりなし";
        }
        putTweet(StarTweetContent(time0, msg, stars(index)._3));
        index += 1;
        if (index == stars.size) {
          index = 0;
        }
      }
    }
  }
}

//==============================================================================
// なにもツイートのない日付

{
  (0 until period).foreach { day =>
    val time = startTime + day;
    val date = TimeLib.modifiedJulianDayToStringJSTDate(time);
    if (getTweets(time).isEmpty) {
      putTweet(time, "#empty".format(date));
    }
  }
}

//==============================================================================
// ツイート出力

(0 until period).foreach { day =>
  getTweets(startTime + day).tweets.foreach { case tc =>
    val time = tc.time;
    if (time >= startTime && time < endTime) {
      val msg = tc.message + tc.hashtags.map(" #" + _).mkString;
      println("%s %s".format(TimeLib.modifiedJulianDayToStringJST(time), msg));
    }
  }
}

//==============================================================================
