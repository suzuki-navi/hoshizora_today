
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
    if (lng < 0) lng + Const.PI2 else lng;
  }

  def xyzToLat(xyz: Array[Double]): Double = {
    val x = xyz(0);
    val y = xyz(1);
    val z = xyz(2);
    val r = Math.sqrt(x * x + y * y);
    val lat = Math.atan2(z, r);
    lat;
  }

  // 地平座標系から方位角
  // 地平座標系 X:西 Y:南 Z:天頂
  def xyzToAzi(xyz: Array[Double]): Double = {
    val x = xyz(0);
    val y = xyz(1);
    val azi = Math.atan2(-x, -y);
    if (azi < 0) azi + Const.PI2 else azi;
  }

  // 地平座標系から高度
  def xyzToAlt(xyz: Array[Double]): Double = {
    xyzToLat(xyz);
  }

  def calcLngDiff(lng1: Double, lng2: Double): Double = {
    val d = lng1 - lng2;
    if (d < 0) d + Const.PI2 else d;
  }

}

