#  bash script to copy files from computer to was s3 bucket
#
path="/media/jon/Seagate_5TB/OPRAM_results/data" # the path of the directory where the files and directories that need to be copied are located
s3Dir="s3://ucd-pest-risk/irish-hectad-map/data" # the s3 bucket path

for entry in "$path"/base0*; do
    name=`echo $entry | sed 's/.*\///'`  # getting the name of the file or directory
#    if [[ -d  $entry ]]; then  # if it is a directory
#        aws s3 cp  --recursive --exclude "*" --include "*.csv" "$name" "$s3Dir/$name/"
#    fi

  echo $name
done
