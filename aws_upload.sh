#  bash script to copy files from computer to was s3 bucket
#
#path="/media/jon/Seagate_5TB/OPRAM_results/data" # the path containing the data
path="/Users/jon/Desktop/gisdata" # the path containing the data
s3Dir="s3://ucd-pest-risk/irish-hectad-map/gisdata" # the s3 bucket path

for entry in "$path"/*; do
    name=`echo $entry | sed 's/.*\///'`  # getting the name of the file or directory
    if [[ -d  "$entry" ]]; then  # if it is a directory
        aws s3 cp  --recursive --exclude "*" --include "*.tif" "$entry" "$s3Dir/$name/"
    fi

done
