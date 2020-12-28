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
  startTime = Time.parse("2020-12-28T00:00:00.000+09:00")
  endTime   = Time.parse("2021-02-01T00:00:00.000+09:00")

  $stdout.sync = true

  puts("# #{startTime} - #{endTime}")
  mkcal(startTime, endTime)
end

####################################################################################################

def mkcal(startTime, endTime)
  sunTerm24str = [
    "春分", "清明", "穀雨", "立夏", "小満", "芒種", "夏至", "小暑", "大暑", "立秋", "処暑", "白露",
    "秋分", "寒露", "霜降", "立冬", "小雪", "大雪", "冬至", "小寒", "大寒", "立春", "雨水", "啓蟄",
  ]

  moonTermStr = [
    "新月", "上弦の月", "満月", "下弦の月",
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
    sunLng = calcLng(sunTE)

    sunDistance = calcDistance(sun)

    sunTerm24b = (sunLng * $pi57 / 15).to_i
    if sunTerm24 != sunTerm24b
      if sunTerm24 >= 0
        msg = "#{sunTerm24str[sunTerm24b]}。太陽の黄経が#{sunTerm24b * 15}°です"
        if sunTerm24b % 6 != 0
          msg = "二十四節気の" + msg
        end
        mkcalPuts(time - 2700, msg)
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
    moonLng = calcLng(moonTE)

    moonSunLng = calcLngDiff(moonLng, sunLng)
    moonSunLng360 = moonSunLng * $pi57
    moonTermB = (moonSunLng360 / 90).to_i

    if moonTerm != moonTermB
      if moonTerm >= 0
        mkcalPuts(time - 2700, "#{moonTermStr[moonTermB]}")
      end
      moonTerm = moonTermB
    end

    moonConsFlag = false

    if moonSunLng360 >= 90 && moonSunLng360 < 225 && timeJST.hour == 21
      moonConsFlag = true
      moonCons = j2000ToConstellations(moon)
      if moonCons != nil
        mkcalPuts(time - 600, "月は#{moonCons}にいます")
      end
    end

    ################################
    # 火星
    ################################

    mars = ephXyz(:mars, time)
    marsTE = xyzGCRSToTrueEcliptic(mars, time)
    marsLng = calcLng(marsTE)

    marsSunLng = calcLngDiff(marsLng, sunLng)
    marsSunLng360 = marsSunLng * $pi57

    if marsSunLng360 >= 90 && marsSunLng360 < 225 && timeJST.wday == 2 && timeJST.hour == 21 && !moonConsFlag
      marsCons = j2000ToConstellations(mars)
      if marsCons != nil
        mkcalPuts(time - 600, "火星は#{marsCons}にいます")
      end
    end

    ################################

    if timeJST.hour == 0
      puts("# #{timeJST.to_s.slice(0, 10)}")
    end
    timeJST = timeJST + 3600
  end
end

####################################################################################################

$constellationsMap = {
  [ 25, 10] => "うお座の東側の魚(アンドロメダ座の南)のしっぽ付近",
  [ 30, 10] => "おひつじ座の頭とくじら座の頭の間",
  [ 35, 10] => "おひつじ座とくじら座の頭の間",
  [ 40, 10] => "おひつじ座とくじら座の頭の間",
  [ 45, 15] => "おひつじ座のしっぽ側",
  [ 60, 15] => "おうし座ヒアデスの西",
  [ 65, 20] => "おうし座とペルセウス座とぎょしゃ座の間",
  [ 70, 20] => "おうし座ヒアデスの北東",
  [ 75, 20] => "おうし座の角付近",
  [ 80, 20] => "おうし座の角とぎょしゃ座付近",
  [ 85, 20] => "おうし座の角とオリオン座の腕付近",
  [ 90, 20] => "ふたご座の西側の子の足元とオリオン座の腕付近",
  [ 95, 20] => "ふたご座の西側の子の足元付近",
  [100, 20] => "ふたご座の2人の足元付近",
  [100, 25] => "ふたご座の西側の子の胴体付近",
  [105, 20] => "ふたご座の東側の子の腰付近",
  [110, 20] => "ふたご座の東側の子の胴体付近",
  [115, 20] => "ふたご座ポルックスの南でふたご座の東",
  [120, 20] => "かに座の西",
  [125, 20] => "かに座",
  [130, 20] => "かに座",
  [135, 20] => "かに座の東",
  [140, 15] => "しし座とかに座の間",
  [145, 15] => "しし座の西",
  [150, 15] => "しし座の肩付近",
  [155, 15] => "しし座の中央",
  [160, 10] => "しし座の腹付近",
  [165, 10] => "しし座の後ろ足付近",
  [175,  5] => "おとめ座の頭付近",
  [190,  0] => "おとめ座の胴体付近",
  [200,-10] => "おとめ座スピカの北",
}

def j2000ToConstellations(xyz)
  lng = calcLng(xyz) * $pi57
  lat = calcLat(xyz) * $pi57
  lng5 = (lng / 5).to_i * 5
  lat5 = ((lat + 90) / 5).to_i * 5 - 90
  key = [lng5, lat5]
  cons = $constellationsMap[key]
  if cons == nil
    raise "#{lng5}(#{lng5/15}h+#{lng5%15*4}m), #{lat5}"
  end
  cons
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
    elsif target == :mars
      target_id = 4
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

def calcLng(xyz)
  lng = Math.atan2(xyz[1], xyz[0])
  lng += $pi2 if (lng < 0)
  lng
end

def calcLat(xyz)
  x = xyz[0]
  y = xyz[1]
  z = xyz[2]
  xy = Math.sqrt(x * x + y * y)
  lat = Math.atan2(z, xy)
  lat
end

def calcLngDiff(lngTarget, lngCenter)
  d = lngTarget - lngCenter
  d += $pi2 if d < 0
  d
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

