
set -Ceu
set -o pipefail

cd $(dirname $0)

. ../period.sh

bundle exec ruby ./main.rb $start $end | tee ../var/holiday.txt

