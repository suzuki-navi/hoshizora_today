
set -Ceu
set -o pipefail

cd $(dirname $0)

if [ ! -e src/main/resources/ascp1950.430 ]; then
    wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de430/ascp1950.430 -O src/main/resources/ascp1950.430
fi

exit

echo "Starting calculating..."

time sbt run

echo "Finished calculating..."

cat ../data.txt | perl -nle '/^[-:T0-9]+ #/ and print' >&2

