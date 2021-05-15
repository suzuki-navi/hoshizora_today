
object Lib2 { // TODO 削除したい

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

  private val PI = Math.PI;
  private val PI2 = Math.PI * 2.0;

}

