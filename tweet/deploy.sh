
set -Ceu
set -o pipefail

cd $(dirname $0)

serverless deploy --verbose

