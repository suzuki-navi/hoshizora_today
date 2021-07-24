
class Hcs(lng: Double, lat: Double) {

  def siderealTimeFromUtc(utc: Double): Double = {
    val ut1 = utc; // 近似的
    siderealTimeFromUt1(ut1);
  }

  def siderealTimeFromUt1(ut1: Double): Double = {
    val gmst = TimeLib.mjdut1ToGmst(ut1);
    val gast = gmst; // 近似
    val s = gast + lng;
    if (s >= Const.PI2) {
      s - Const.PI2;
    } else {
      s;
    }
  }

  def siderealTimeToUtc(sid: Double, time0: Double): Double = {
    siderealTimeToUt1(sid, time0); // 近似
  }

  def siderealTimeToUt1(sid: Double, time0: Double): Double = {
    val s = sid - lng;
    val s2 = if (s < 0) {
      s + Const.PI2;
    } else {
      s;
    }
    TimeLib.gmstToMjdut1(s2, time0);
  }

  def trueEquatorialXyzToAziAltFromUtc(xyz: Array[Double], utc: Double): (Double, Double) = {
    val ut1 = utc; // 近似的
    trueEquatorialXyzToAziAltFromUt1(xyz, ut1);
  }

  def trueEquatorialXyzToAziAltFromUt1(xyz: Array[Double], ut1: Double): (Double, Double) = {
    var r: Array[Double] = VectorLib.unitMatrix;
    r = VectorLib.rotateMatrixZ(Const.PI5 - siderealTimeFromUt1(ut1), r);
    r = VectorLib.rotateMatrixX(Const.PI5 - lat, r);
    // 地平座標系 X:西 Y:南 Z:天頂
    val hcs = VectorLib.multiplyMV(r, xyz);
    val x = hcs(0);
    val y = hcs(1);
    val z = hcs(2);
    val azi1 = Math.atan2(-x, -y);
    val azi = if (azi1 < 0) azi1 + Const.PI2 else azi1; // 北0°、東90°、南180°、西270°
    val xy = Math.sqrt(x * x + y * y);
    val alt = Math.atan2(z, xy);
    (azi, alt);
  }

}

object Hcs {

  def aziAltToNaturalString(azi: Double, alt: Double): String = {
    val altHor = -0.90 / Const.PI57;
    val alt360 = alt * Const.PI57;
    if (alt360 >= 80) {
      "天頂付近";
    } else {
      val azi360 = azi * Const.PI57;
      val s1 = if (azi360 <= 22.5 || azi360 >= 337.5) {
        "北の";
      } else if (azi360 < 67.5) {
        "北東の";
      } else if (azi360 > 292.5) {
        "北西の";
      } else if (azi360 <= 112.5) {
        "東の";
      } else if (azi360 >= 247.5) {
        "西の";
      } else if (azi360 < 157.5) {
        "南東の";
      } else if (azi360 > 202.5) {
        "南西の";
      } else {
        "南の";
      }
      val s2 = if (alt < altHor) {
        "地平線の下";
      } else if (alt360 < 10.0) {
        "空地平線近く";
      } else if (alt360 < 30.0) {
        "空低く";
      } else if (alt360 >= 60.0) {
        "空高く";
      } else {
        "空";
      }
      s1 + s2;
    }
  }

  def aziAltToNaturalStringC(azi: Double, alt: Double): (String, Double) = {
    val alt360 = alt * Const.PI57;
    if (alt360 >= 80) {
      ("天頂", 0.0);
    } else if (alt360 < 30) {
      ("", 0.0);
    } else {
      val azi360 = azi * Const.PI57;
      if (azi360 <= 45) {
        ("北の空", -alt);
      } else if (azi360 >= 315) {
        ("北の空", -alt);
      } else if (azi360 <= 135) {
        ("東の空", -alt);
      } else if (azi360 >= 225) {
        ("西の空", -alt);
      } else {
        ("南の空", -alt);
      }
    }
  }

}

