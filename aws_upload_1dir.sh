#  bash script to copy files from computer to aws s3 bucket
# This copies a single
#
#path="/media/jon/Seagate_5TB/OPRAM_results/data" # the path containing the data
path="/Users/jon/OPRAM/results/halyomorpha_halys/" # the path containing the data
s3Dir="s3://ucd-pest-risk/irish-hectad-map/data/halyomorpha_halys/" # the s3 bucket path


    
# Always use --dryrun before commiting to run a command
#aws s3 cp  --dryrun --recursive --exclude "*" --include "*.csv" --exclude "average*.csv" "$path" "$s3Dir"
aws s3 cp  --recursive --exclude "*" --include "*.csv" --exclude "average*.csv" "$path" "$s3Dir"


# List contents
#aws s3 ls "$s3Dir"  
