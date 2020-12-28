require 'date'
require 'time'
require 'holiday_japan'
require 'mk_time'
require 'eph_bpn'
require 'eph_jpl'

$data_path = '../var/de430/ssd.jpl.nasa.gov/pub/eph/planets/Linux/de430/linux_p1550p2650.430'

$pi2 = 2.0 * Math::PI
$pi5 = 0.5 * Math::PI
$pi57 = 180.0 / Math::PI

def main(args)
  startTime = Time.parse("2021-01-01T00:00:00.000+09:00")
  endTime   = Time.parse("2021-01-06T00:00:00.000+09:00")

  $stdout.sync = true

  puts("# #{startTime} - #{endTime}")
  mkcal(startTime, endTime)
end

####################################################################################################

def mkcal(startTime, endTime)
  sunTerm24str = [
    "春分",
    "清明",
    "穀雨",
    "立夏",
    "小満",
    "芒種",
    "夏至",
    "小暑",
    "大暑",
    "立秋",
    "処暑",
    "白露",
    "秋分",
    "寒露",
    "霜降",
    "立冬",
    "小雪",
    "大雪",
    "冬至",
    "小寒",
    "大寒",
    "立春",
    "雨水",
    "啓蟄",
  ]

  moonTermStr = [
    "新月",
    "上弦の月",
    "満月",
    "下弦の月",
  ]

  timeJST = startTime.localtime("+09:00")

  sunTerm24 = -1
  sunDistanceFlag = nil
  sunPrevDistance = -1

  moonTerm = -1

  while timeJST <= endTime
    time = Time.at(timeJST)
    time.localtime("UTC")

    if timeJST.hour == 9
      date = Date.new(timeJST.year, timeJST.month, timeJST.mday)
      holiday = HolidayJapan.name(date)
      if holiday != nil
        mkcalPuts(time + 300, "#{holiday}")
      end
    end

    ################################
    # 太陽
    ################################

    sun = ephXyz(:sun, time)
    sunTE = xyzGCRSToTrueEcliptic(sun, time)
    sunLng = Math.atan2(sunTE[1], sunTE[0])
    sunLng += $pi2 if (sunLng < 0)

    sunDistance = calcDistance(sun)

    sunTerm24b = (sunLng * $pi57 / 15).to_i
    if sunTerm24 != sunTerm24b
      if sunTerm24 >= 0
        mkcalPuts(time - 2700, "#{sunTerm24str[sunTerm24b]}。太陽の黄経が#{sunTerm24b * 15}°です。")
      end
      sunTerm24 = sunTerm24b
    end

    if sunDistanceFlag != nil
      if sunDistanceFlag && sunPrevDistance > sunDistance
        mkcalPuts(time - 4500, "地球が遠日点通過")
        sunDistanceFlag = false
      elsif !sunDistanceFlag && sunPrevDistance < sunDistance
        mkcalPuts(time - 4500, "地球が近日点通過")
        sunDistanceFlag = true
      end
    elsif sunPrevDistance > 0
      if sunPrevDistance < sunDistance
        sunDistanceFlag = true
      else
        sunDistanceFlag = false
      end
    end
    sunPrevDistance = sunDistance

    ################################
    # 月
    ################################

    moon = ephXyz(:moon, time)
    moonTE = xyzGCRSToTrueEcliptic(moon, time)
    moonLng = Math.atan2(moonTE[1], moonTE[0])
    moonLng += $pi2 if (moonLng < 0)

    moonSunLng = moonLng - sunLng
    moonSunLng += $pi2 if (moonSunLng < 0)
    moonSunLng360 = moonSunLng * $pi57
    moonTermB = (moonSunLng360 / 90).to_i

    if moonTerm != moonTermB
      if moonTerm >= 0
        mkcalPuts(time - 2700, "#{moonTermStr[moonTermB]}")
      end
      moonTerm = moonTermB
    end

    ################################

    if timeJST.hour == 0
      puts("# #{timeJST.to_s.slice(0, 10)}")
    end
    timeJST = timeJST + 3600
  end
end

####################################################################################################

$ephXyzCache = {}

def ephXyz(target, time)
  # target: sun / moon
  # time: UTC DateTime型
  key = [target, time]
  ret = $ephXyzCache[key]
  if ret == nil
    ttjd = MkTime::Calc.new(MkTime::Calc.new(time).tt).jd
    if target == :sun
      target_id = 11
    elsif target == :moon
      target_id = 10
    else
      raise
    end
    ephJpl = EphJpl.new($data_path, target_id, 3, ttjd)
    value = ephJpl.calc
    x = value[0]
    y = value[1]
    z = value[2]
    ret = [x, y, z]
    $ephXyzCache[key] = ret
  end
  return ret
end

def xyzGCRSToTrueEquatorial(xyz, time)
  bpn = EphBpn::Ephemeris.new(time)
  bpn.apply_bias_prec_nut(xyz)
end

def xyzGCRSToTrueEcliptic(xyz, time)
  bpn = EphBpn::Ephemeris.new(time)
  xyz2 = bpn.apply_bias_prec_nut(xyz)
  rotationX(xyz2, -bpn.eps)
end

def xyzEquatorialToEcliptic(xyz)
  return rotationX(xyz, -0.4090926)
end

def calcDistance(xyz)
  x = xyz[0]
  y = xyz[1]
  z = xyz[2]
  Math.sqrt(x * x + y * y + z * z)
end

def rotationX(xyz, th)
  x = xyz[0]
  y = xyz[1]
  z = xyz[2]
  cos = Math.cos(th)
  sin = Math.sin(th)
  return [x, y * cos - z * sin, y * sin + z * cos]
end

####################################################################################################

def mkcalPuts(time, msg)
  puts("             #{mkcalTimeStr(time)} #{msg}")
end

def mkcalTimeStr(time)
  (time + 9 * 3600).strftime('%Y-%m-%d %H:%M:%S')
end

####################################################################################################

main(ARGV)

