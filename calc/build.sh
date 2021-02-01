
set -Ceu
set -o pipefail

cd $(dirname $0)

. ./period.sh

echo "started calculating..."

scala ./main.scala $start $end >| ../var/statuses.txt

echo "finished calculating..."

cat diff.txt | grep '^20' | LC_ALL=C sort >| ../var/diff.txt

cat ../var/statuses.txt | node patch.js >| ../var/data.txt

mv ../var/data.txt ../data.txt

