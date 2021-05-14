
set -Ceu
set -o pipefail

cd $(dirname $0)

echo "Starting calculating..."

time scala -deprecation ./main.scala

echo "Finished calculating..."

cat ../data.txt | perl -nle '/^[-:T0-9]+ #/ and print' >&2

