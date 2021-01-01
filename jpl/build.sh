
set -Ceu
set -o pipefail

cd $(dirname $0)

. ../period.sh

if [[ ! -e ../var/ssd.jpl.nasa.gov ]]; then
    (
        cd ../var
        wget -r ftp://ssd.jpl.nasa.gov/pub/eph/planets/Linux/de430/
    )
fi

(
    echo $start
    echo $end2
) >| ../var/jpl.period.txt.new

if [[ -e ../var/jpl.period.txt ]] && cmp ../var/jpl.period.txt ../var/jpl.period.txt.new >/dev/null; then
    exit 0
fi

bundle exec ruby ./main.rb $start $end2

mv ../var/jpl.dat.new ../var/jpl.dat
mv ../var/jpl.period.txt.new ../var/jpl.period.txt

