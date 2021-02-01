
set -Ceu
set -o pipefail

cd $(dirname $0)

bash calc/build.sh

#aws_profile=$(cat tweet/config.yml | yq -r '.profile')
#s3_path="s3://$(cat tweet/config.yml | yq -r '.data_s3_bucket')/$(cat tweet/config.yml | yq -r '.data_s3_key')"
#aws --profile $aws_profile s3 cp data.txt $s3_path

