
object JplData {

  private val EMRAT = 81.3005690741906200; // 地球と月の質量比
  private val EMRAT1 = 1.0 / (1.0 + EMRAT);

  private sealed trait JplPlanet {
    def dataOffset: Int;
    def subPeriodCount: Int;
    def coefficientCount: Int;
  }
  private object MERCURY_JPL extends JplPlanet {
    val dataOffset = 2;
    val subPeriodCount = 4;
    val coefficientCount = 14;
  }
  private object VENUS_JPL extends JplPlanet {
    val dataOffset = 170;
    val subPeriodCount = 2;
    val coefficientCount = 10;
  }
  private object EARTH_MOON_BARYCENTER_JPL extends JplPlanet {
    val dataOffset = 230;
    val subPeriodCount = 2;
    val coefficientCount = 13;
  }
  private object MARS_JPL extends JplPlanet {
    val dataOffset = 308;
    val subPeriodCount = 1;
    val coefficientCount = 11;
  }
  private object JUPITER_JPL extends JplPlanet {
    val dataOffset = 341;
    val subPeriodCount = 1;
    val coefficientCount = 8;
  }
  private object SATURN_JPL extends JplPlanet {
    val dataOffset = 365;
    val subPeriodCount = 1;
    val coefficientCount = 7;
  }
  private object MOON_JPL extends JplPlanet {
    val dataOffset = 440;
    val subPeriodCount = 8;
    val coefficientCount = 13;
  }
  private object SUN_JPL extends JplPlanet {
    val dataOffset = 752;
    val subPeriodCount = 2;
    val coefficientCount = 11;
  }

  sealed trait TargetPlanet;
  object Sun     extends TargetPlanet;
  object Moon    extends TargetPlanet;
  object Mercury extends TargetPlanet;
  object Venus   extends TargetPlanet;
  object Mars    extends TargetPlanet;
  object Jupiter extends TargetPlanet;
  object Saturn  extends TargetPlanet;

  private[this] val jplData: IndexedSeq[IndexedSeq[Double]] = loadJplData();

  private[this] def loadJplData(): IndexedSeq[IndexedSeq[Double]] = {
    val dataPath: String = "ascp1950.430";
    val data1 = {
      val sc = new java.util.Scanner(new java.io.BufferedInputStream(
        getClass.getClassLoader.getResourceAsStream(dataPath)));
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

  private[this] def calcChebyshev(cs: IndexedSeq[Double], x: Double): Double = {
    // x は -1.0 <= x <= +1.0
    val n = cs.size;
    if (n == 1) {
      cs(0);
    } else {
      val x2 = 2.0 * x;
      var result: Double = cs(0) + cs(1) * x;
      var i: Int = 2;
      var t0: Double = 1.0;
      var t1: Double = x;
      while (i < n) {
        val t2 = x2 * t1 - t0;
        result += cs(i) * t2;
        t0 = t1;
        t1 = t2;
        i = i + 1;
      }
      result;
    }
  }

  private[this] def calcPlanetPosition(time: Double, planet: JplData.JplPlanet): Array[Double] = {
    val (periodData, time2) = getPeriodData(time);
    val time3 = time2 * planet.subPeriodCount;
    val subPeriodIndex = time3.toInt;
    val time4 = 2.0 * (time3 - subPeriodIndex) - 1.0;
    val coefficientCount = planet.coefficientCount;
    val offset = 3 * coefficientCount * subPeriodIndex + planet.dataOffset;
    val ret = new Array[Double](3);
    ret(0) = calcChebyshev(periodData.slice(offset + 0 * coefficientCount, offset + 1 * coefficientCount), time4);
    ret(1) = calcChebyshev(periodData.slice(offset + 1 * coefficientCount, offset + 2 * coefficientCount), time4);
    ret(2) = calcChebyshev(periodData.slice(offset + 2 * coefficientCount, offset + 3 * coefficientCount), time4);
    ret;
  }

  private[this] def calcEarthPosition(time: Double): Array[Double] = {
    val em = calcPlanetPosition(time, JplData.EARTH_MOON_BARYCENTER_JPL);
    val moon = calcPlanetPosition(time, JplData.MOON_JPL);
    val emrat1 = JplData.EMRAT1;
    val ret1 = new Array[Double](3);
    ret1(0) = em(0) - moon(0) * emrat1;
    ret1(1) = em(1) - moon(1) * emrat1;
    ret1(2) = em(2) - moon(2) * emrat1;
    ret1;
  }

  def calcPlanetFromEarth(time: Double, targetPlanet: JplData.TargetPlanet): Array[Double] = {
    val target = targetPlanet match {
      case JplData.Sun     => calcPlanetPosition(time, JplData.SUN_JPL);
      case JplData.Moon    => return calcPlanetPosition(time, JplData.MOON_JPL);
      case JplData.Mercury => calcPlanetPosition(time, JplData.MERCURY_JPL);
      case JplData.Venus   => calcPlanetPosition(time, JplData.VENUS_JPL);
      case JplData.Mars    => calcPlanetPosition(time, JplData.MARS_JPL);
      case JplData.Jupiter => calcPlanetPosition(time, JplData.JUPITER_JPL);
      case JplData.Saturn  => calcPlanetPosition(time, JplData.SATURN_JPL);
    }
    val earth = calcEarthPosition(time);
    VectorLib.minus(target, earth);
  }

}

