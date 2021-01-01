
set -Ceu
set -o pipefail

cd $(dirname $0)

cat ../var/data.1.txt | node patch.js >| ../var/data.2.txt

