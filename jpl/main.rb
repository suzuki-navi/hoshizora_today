require 'time'
require 'mk_time'
require 'eph_bpn'
require 'eph_jpl'
require 'mk_greenwich'

$jpl_data_path = '../var/ssd.jpl.nasa.gov/pub/eph/planets/Linux/de430/linux_p1550p2650.430'
$planets_data_path = '../var/jpl.dat.new'

$pi2 = 2.0 * Math::PI
$pi5 = 0.5 * Math::PI
$pi57 = 180.0 / Math::PI

$tokyoLng = 139.7 / $pi57
$tokyoLat = 35.7 / $pi57

def main(args)
  startTime = Time.parse(args[0] + "T00:00:00+09:00").localtime("UTC")
  endTime = Time.parse(args[1] + "T00:00:00+09:00").localtime("UTC")
  File.open($planets_data_path, "w") do |fp|
    time = startTime
    while time < endTime
      puts("time: #{time}")
      data = calcData(time)
      buf = data.pack("G#{data.length}")
      fp.write(buf)
      time = time + 600
    end
  end
end

def calcData(time)
  # time: UTC DateTime型
  ttjd = MkTime::Calc.new(MkTime::Calc.new(time).tt).jd
  bpn = EphBpn::Ephemeris.new(time)
  gw = MkGreenwich::Greenwich.new(time)

  meps = -bpn.eps
  mgast = $pi5 -gw.gast - $tokyoLng
  mlat = $pi5 - $tokyoLat

  data = []
  [:sun, :moon, :mercury, :venus, :mars, :jupiter, :saturn].each do |target|
    xyz = calcXyz(target, ttjd)
    trueEquatorialXyz = bpn.apply_bias_prec_nut(xyz)
    trueEclipticXyz = rotationX(trueEquatorialXyz, meps)
    hcsXyz = rotationX(rotationZ(trueEquatorialXyz, mgast), mlat) # 地平座標系 X:西 Y:南 Z:天頂
    distance = calcDistance(xyz)
    j2000Lng = calcLng(xyz) # J2000.0をxyzで近似
    j2000Lat = calcLat(xyz)
    trueEclipticLng = calcLng(trueEclipticXyz)
    hcsAzi = calcAzi(hcsXyz)
    hcsAlt = calcLat(hcsXyz)
    data.push(distance, j2000Lng, j2000Lat, trueEclipticLng, hcsAzi, hcsAlt)
  end
  return data
end

def calcXyz(target, ttjd)
  if target == :sun
    target_id = 11
  elsif target == :moon
    target_id = 10
  elsif target == :mercury
    target_id = 1
  elsif target == :venus
    target_id = 2
  elsif target == :mars
    target_id = 4
  elsif target == :jupiter
    target_id = 5
  elsif target == :saturn
    target_id = 6
  else
    raise
  end
  ephJpl = EphJpl.new($jpl_data_path, target_id, 3, ttjd)
  value = ephJpl.calc
  x = value[0]
  y = value[1]
  z = value[2]
  [x, y, z]
end

def rotationX(xyz, th)
  x = xyz[0]
  y = xyz[1]
  z = xyz[2]
  cos = Math.cos(th)
  sin = Math.sin(th)
  return [x, y * cos - z * sin, y * sin + z * cos]
end

def rotationZ(xyz, th)
  x = xyz[0]
  y = xyz[1]
  z = xyz[2]
  cos = Math.cos(th)
  sin = Math.sin(th)
  return [x * cos - y * sin, x * sin + y * cos, z]
end

def calcLng(xyz)
  lng = Math.atan2(xyz[1], xyz[0])
  lng += $pi2 if (lng < 0)
  lng
end

def calcAzi(xyz)
  lng = Math.atan2(-xyz[0], -xyz[1])
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

def calcDistance(xyz)
  x = xyz[0]
  y = xyz[1]
  z = xyz[2]
  Math.sqrt(x * x + y * y + z * z)
end

main(ARGV)

