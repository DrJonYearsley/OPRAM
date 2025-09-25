#  bash script to delete specific files from aws s3 bucket
#

s3Dir="s3://ucd-pest-risk/irish-hectad-map/data" # the s3 bucket path



#aws s3 ls "$s3Dir"  --recursive
#aws s3 ls --recursive --exclude "*" --include "average*"  "$s3Dir"

# Always do a dryrun before running a delete command!!
#aws s3 rm "$s3Dir"  --dryrun --recursive --exclude "*" --include "*average*"


#aws s3 rm "$s3Dir"  --recursive --exclude "*" --include "*average*"
