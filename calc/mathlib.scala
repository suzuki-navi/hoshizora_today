
object MathLib {

  def circleAdd(a: Double, b: Double): Double = { // TODO 削除したい
    val d = a + b;
    if (d >= 3 * PI) {
      d - 2 * PI2;
    } else if (d >= PI) {
      d - PI2;
    } else if (d < -3 * PI) {
      d + 2 * PI2;
    } else if (d < -PI) {
      d + PI2;
    } else {
      d;
    }
  }

  // 0~2pi - 0~2pi => -pi~+pi
  def circleDiff1(a: Double, b: Double): Double = {
    val d = a - b;
    if (d >= PI) {
      d - PI2;
    } else if (d < -PI) {
      d + PI2;
    } else {
      d;
    }
  }

  // 0以上になるインデックスを探す
  def binarySearchBy[A](seq: IndexedSeq[A])(f: A => Double): Int = {
    @scala.annotation.tailrec
    def sub(start: Int, end: Int): Int = {
      if (end - start == 1) {
        start;
      } else if (end - start == 2) {
        if (f(seq(start)) >= 0.0) {
          start;
        } else {
          start + 1;
        }
      } else {
        val mid = (end + start) / 2;
        val midf = f(seq(mid));
        if (midf <= 0.0) {
          sub(mid, end);
        } else {
          sub(start, mid + 1);
        }
      }
    }
    sub(0, seq.size);
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

  private val PI = Math.PI;
  private val PI2 = Math.PI * 2.0;
  private val PI2R = 1.0 / PI2;

}

