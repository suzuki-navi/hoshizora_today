
set -Ceu
set -o pipefail

cd $(dirname $0)

. ./period.sh

echo "Starting calculating..."

scala ./main.scala $start0 $end0 >| ../var/statuses.txt

echo "Finished calculating..."

cat diff.txt | grep '^20' | LC_ALL=C sort >| ../var/diff.txt

cat ../var/statuses.txt | node patch.js $start1 $end1 >| ../var/data.txt

mv ../var/data.txt ../data.txt

