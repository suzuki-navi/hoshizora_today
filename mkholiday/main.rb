require 'date'
require 'time'
require 'holiday_japan'

def main(args)
  startTime = Time.parse(args[0])
  endTime = Time.parse(args[1])
  puts "# time: #{startTime} - #{endTime}"
  time = startTime
  while time < endTime
    date = time.to_date
    holiday = HolidayJapan.name(date)
    if holiday != nil
      puts "#{date}T09:05 #{holiday}"
    end
    time = time + 86400
  end
end

main(ARGV)

