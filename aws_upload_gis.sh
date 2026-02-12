#  bash script to copy files from computer to aws s3 bucket
# This copies multiple directories
#
path="/media/jon/Seagate_5TB/OPRAM_sonic/gisdata" # the path containing the data
#path="/Users/jon/OPRAM/gisdata" # the path containing the data
s3Dir="s3://ucd-pest-risk/irish-hectad-map/gisdata" # the s3 bucket path

for entry in "$path"/*; do
  name=`echo $entry | sed 's/.*\///'`  # getting the name of the file or directory
    
  if [[ -d  "$entry" ]]; then  # if it is a directory
      # Always use --dryrun before commiting to run a command
#      aws s3 cp  --dryrun --recursive --exclude "*" --include "*.tif"  "$entry" "$s3Dir/$name/"
      aws s3 cp  --recursive --exclude "*" --include "*.tif"  "$entry" "$s3Dir/$name/"

  fi
done


# List contents
#aws s3 ls "$s3Dir"  
