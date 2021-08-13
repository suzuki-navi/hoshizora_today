
object Acs {

  object Icrs {

    def calcPlanetXyz(utc: Double, targetPlanet: JplData.TargetPlanet): Array[Double] = {
      val tdb = TimeLib.mjdutcToTdb(utc);
      val planet = JplData.calcPlanetFromEarth(tdb, targetPlanet);
      planet;
    }

  }

  object Ecliptic2000 {

    def calcPlanetLng(utc: Double, targetPlanet: JplData.TargetPlanet): Double = {
      val tdb = TimeLib.mjdutcToTdb(utc);
      val bpnMatrix = Bpn.icrsToMeanEclipticMatrix2000;
      val planet = JplData.calcPlanetFromEarth(tdb, targetPlanet);
      val planet2 = VectorLib.multiplyMV(bpnMatrix, planet);
      val planetLng = VectorLib.xyzToLng(planet2);
      planetLng;
    }

  }

  object EclipticTrue {

    def calcPlanetLng(utc: Double, targetPlanet: JplData.TargetPlanet): Double = {
      val tdb = TimeLib.mjdutcToTdb(utc);
      val bpnMatrix = Bpn.icrsToTrueEclipticMatrix(tdb);
      val planet = JplData.calcPlanetFromEarth(tdb, targetPlanet);
      val planet2 = VectorLib.multiplyMV(bpnMatrix, planet);
      val planetLng = VectorLib.xyzToLng(planet2);
      planetLng;
    }

  }

  object Horizontal {

    // 地平座標系 X:西 Y:南 Z:天頂

    def calcPlanetAzi(utc: Double, targetPlanet: JplData.TargetPlanet, hcs: Hcs): Double = {
      val tdb = TimeLib.mjdutcToTdb(utc);
      val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
      val planet = JplData.calcPlanetFromEarth(tdb, targetPlanet);
      val planet2 = VectorLib.multiplyMV(bpnMatrix, planet);

      var r: Array[Double] = VectorLib.unitMatrix;
      r = VectorLib.rotateMatrixZ(Const.PI5 - hcs.siderealTimeFromUt1(utc), r);
      r = VectorLib.rotateMatrixX(Const.PI5 - hcs.lat, r);
      val planet3 = VectorLib.multiplyMV(r, planet2);
      VectorLib.xyzToAzi(planet3);
    }

    def calcPlanetAlt(utc: Double, targetPlanet: JplData.TargetPlanet, hcs: Hcs): Double = {
      val tdb = TimeLib.mjdutcToTdb(utc);
      val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
      val planet = JplData.calcPlanetFromEarth(tdb, targetPlanet);
      val planet2 = VectorLib.multiplyMV(bpnMatrix, planet);

      var r: Array[Double] = VectorLib.unitMatrix;
      r = VectorLib.rotateMatrixZ(Const.PI5 - hcs.siderealTimeFromUt1(utc), r);
      r = VectorLib.rotateMatrixX(Const.PI5 - hcs.lat, r);
      val planet3 = VectorLib.multiplyMV(r, planet2);
      VectorLib.xyzToAlt(planet3);
    }

    def calcPlanetAziAlt(utc: Double, targetPlanet: JplData.TargetPlanet, hcs: Hcs): (Double, Double) = {
      val tdb = TimeLib.mjdutcToTdb(utc);
      val bpnMatrix = Bpn.icrsToTrueEquatorialMatrix(tdb);
      val planet = JplData.calcPlanetFromEarth(tdb, targetPlanet);
      val planet2 = VectorLib.multiplyMV(bpnMatrix, planet);

      var r: Array[Double] = VectorLib.unitMatrix;
      r = VectorLib.rotateMatrixZ(Const.PI5 - hcs.siderealTimeFromUt1(utc), r);
      r = VectorLib.rotateMatrixX(Const.PI5 - hcs.lat, r);
      val planet3 = VectorLib.multiplyMV(r, planet2);
      (VectorLib.xyzToAzi(planet3), VectorLib.xyzToAlt(planet3));
    }

  }

}

