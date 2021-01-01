
set -Ceu
set -o pipefail

cd $(dirname $0)

. ../period.sh

(
    echo $start
    echo $end1
) >| ../var/holiday.period.txt.new

if [[ -e ../var/holiday.period.txt ]] && cmp ../var/holiday.period.txt ../var/holiday.period.txt.new >/dev/null; then
    exit 0
fi

bundle exec ruby ./main.rb $start $end1 | tee ../var/holiday.txt

mv ../var/holiday.period.txt.new ../var/holiday.period.txt

