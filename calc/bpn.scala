
class Bpn(nutLsDataPath: String, nutPlDataPath: String) {

  private val PI2 = Math.PI * 2.0;
  private val PI57 = 180.0 / Math.PI;
  private val PI_AS2R = Math.PI / (3600 * 180);

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

  val icrsToMeanEclipticMatrix2000: Array[Double] = {
    icrsToMeanEclipticMatrix(TimeLib.stringToModifiedJulianDay("2000-01-01T00:00:00+00:00"));
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

  private def loadNutLs(path: String): IndexedSeq[IndexedSeq[Double]] = {
    // L L' F D Om   PS PST PC EC ECT ES
    val data1 = {
      val sc = new java.util.Scanner(new java.io.BufferedInputStream(new java.io.FileInputStream(path)));
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

  private def loadNutPl(path: String): IndexedSeq[IndexedSeq[Double]] = {
    // L L' F D Om Lm Lv Le LM Lj Ls Lu Ln Pa   PS PC ES EC
    val data1 = {
      val sc = new java.util.Scanner(new java.io.BufferedInputStream(new java.io.FileInputStream(path)));
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

  private val nutLsData: IndexedSeq[IndexedSeq[Double]] = loadNutLs(nutLsDataPath);
  private val nutPlData: IndexedSeq[IndexedSeq[Double]] = loadNutPl(nutPlDataPath);

}

