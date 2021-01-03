
set -Ceu
set -o pipefail

cd $(dirname $0)

cat diff.txt | grep '^20' | LC_ALL=C sort >| ../var/diff.txt

cat ../var/data.1.txt | node patch.js >| ../var/data.2.txt

