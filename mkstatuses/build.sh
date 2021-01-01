
set -Ceu
set -o pipefail

cd $(dirname $0)

. ../period.sh

(
    echo $start
    echo $end1
) >| ../var/statuses.period.txt.new

#if [[ -e ../var/statuses.period.txt ]] && cmp ../var/statuses.period.txt ../var/statuses.period.txt.new >/dev/null; then
#    exit 0
#fi

echo "start mkstatuses..."

scala ./main.scala $start $end1 >| ../var/statuses.txt

mv ../var/statuses.period.txt.new ../var/statuses.period.txt

