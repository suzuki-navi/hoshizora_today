
set -Ceu
set -o pipefail

cd $(dirname $0)

. ../period.sh

echo "start mkstatuses..."

scala ./main.scala $start $end | tee ../var/statuses.txt

