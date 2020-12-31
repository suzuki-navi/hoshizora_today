
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

bundle exec ruby ./main.rb $start $end

