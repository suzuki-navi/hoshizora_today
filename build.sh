
set -Ceu
set -o pipefail

cd $(dirname $0)

if [[ ! -e var/planets.dat || var/planets.dat -ot jpl/main.rb || var/planets.dat -ot period.sh ]]; then
    bash jpl/build.sh
fi

bash mkholiday/build.sh
bash mkstatuses/build.sh
cat var/holiday.txt var/statuses.txt | grep '^20' | LC_ALL=C sort >| data.txt

#aws_profile=$(cat tweet/config.yml | yq -r '.profile')
#s3_path="s3://$(cat tweet/config.yml | yq -r '.data_s3_bucket')/$(cat tweet/config.yml | yq -r '.data_s3_key')"
#aws --profile $aws_profile s3 cp data.txt $s3_path

