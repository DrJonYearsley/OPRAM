#  bash script to delete specific files from aws s3 bucket
#

s3Dir="s3://ucd-pest-risk/irish-hectad-map/data" # the s3 bucket path
txt="base*"


aws s3 ls "$s3Dir"  --recursive  --exclude "$txt"  
#aws s3 ls --recursive --exclude "*" --include "average*"  "$s3Dir"
