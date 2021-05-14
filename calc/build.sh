
set -Ceu
set -o pipefail

cd $(dirname $0)

. ./period.sh

echo "Starting calculating..."

time scala -deprecation ./main.scala $start0  $start1 $end1 $end0

echo "Finished calculating..."

cat ../data.txt | perl -nle '/^[-:T0-9]+ #/ and print' >&2

