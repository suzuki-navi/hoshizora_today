
set -Ceu

cd $(dirname $0)
cd mkcal

bundle exec ruby ./mkcal.rb | tee var/data.txt

cat var/data.txt | grep '^\s*2' | sed -E 's/^\s*(.+)$/\1/' | LC_ALL=C sort >| ../data.txt

