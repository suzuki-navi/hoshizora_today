require 'time'
require 'mk_time'
require 'eph_jpl'

$data_path = '../var/de430/ssd.jpl.nasa.gov/pub/eph/planets/Linux/de430/linux_p1550p2650.430'

$pi2 = 2.0 * Math::PI
$pi5 = 0.5 * Math::PI
$pi57 = 180.0 / Math::PI

$startTime = Time.parse("2020-12-31T15:00:00.000Z")
$endTime = Time.parse("2021-01-31T15:00:00.000Z")

def main
  mkcalMoon
end

def mkcalMoon
  time = $startTime
  v = -1
  while time < $endTime
    wday = time.wday
    sm = lngDistanceSunMoon(time) * $pi57
    v2 = (sm / 90).to_i
    if v != v2
      if v == 0
        mkcalPuts(time - 2700, "上弦の月")
      elsif v == 1
        mkcalPuts(time - 2700, "満月")
      elsif v == 2
        mkcalPuts(time - 2700, "下弦の月")
      elsif v == 3
        mkcalPuts(time - 2700, "新月")
      end
      v = v2
    end
    time = time + 3600
  end
end

def lngDistanceSunMoon(time)
  sun = xyzEquatorialToEcliptic(ephXyz(:sun, time))
  moon = xyzEquatorialToEcliptic(ephXyz(:moon, time))
  sunLng = Math.atan2(sun[1], sun[0])
  moonLng = Math.atan2(moon[1], moon[0])
  d = moonLng - sunLng
  d += $pi2 if (d < 0)
  d
end

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

def xyzEquatorialToEcliptic(xyz)
  return rotationX(xyz, -0.4090926)
end

def rotationX(xyz, th)
  x = xyz[0]
  y = xyz[1]
  z = xyz[2]
  cos = Math.cos(th)
  sin = Math.sin(th)
  return [x, y * cos - z * sin, y * sin + z * cos]
end

def mkcalPuts(time, msg)
  puts("#{mkcalTimeStr(time)} #{msg}")
end

def mkcalTimeStr(time)
  (time + 9 * 3600).strftime('%Y-%m-%d %H:%M:%S')
end

main

