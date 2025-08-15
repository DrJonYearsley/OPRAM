#  bash script to copy files from computer to aws s3 bucket
#
path="/media/jon/Seagate_5TB/OPRAM_results/data" # the path containing the data
#path="/Users/jon/Desktop/gisdata" # the path containing the data
s3Dir="s3://ucd-pest-risk/irish-hectad-map/data" # the s3 bucket path

for entry in "$path"/*; do
  name=`echo $entry | sed 's/.*\///'`  # getting the name of the file or directory
    
  if [[ "$name" != "average"* ]] && [[ "$name" != "base0"* ]];then    # continue if name doesn't contain this text
    if [[ -d  "$entry" ]]; then  # if it is a directory
        aws s3 cp  --recursive --exclude "*" --include "*.csv" "$entry" "$s3Dir/$name/"
    fi
  fi
done
