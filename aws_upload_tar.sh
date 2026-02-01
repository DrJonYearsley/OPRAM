#  bash script to extract tar file into a directory, copy files from extracted file to aws s3 bucket
#  and then remove the extracted files
#
# This copies multiple directories
#
#path="/media/jon/Seagate_5TB/OPRAM_results/data" # the path containing the data
path="/Users/jon/OPRAM/results_tar" # the path containing the data
s3Dir="s3://ucd-pest-risk/irish-hectad-map/data/" # the s3 bucket path

#ls $path

for entry in "$path"/*.tar.gz; do

    # Split tar file into directory and filename
    dir_path=`echo $entry | sed 's/.tar.gz//'`  # getting the name of the file or directory
    name=`echo $dir_path | sed 's/.*\///'`

    echo ""
    echo "======= Extracting ${dir_path} ================"
    echo ""
    
    # Extract tar file
    tar -xf $entry -C $path

    
    if [[ -d  "$dir_path" ]]; then  # if it is a directory
       # Always use --dryrun before commiting to run a command
       aws s3 cp  --dryrun --recursive --exclude "*" --include "*.csv" --exclude "average*.csv" "$dir_path" "$s3Dir/$name/"
    #      aws s3 cp  --recursive --exclude "*" --include "*.csv" --exclude "average*.csv" "$entry" "$s3Dir/$name/"

        # Remove extracted directory
      rm -r $dir_path
    fi


done


# List contents
#aws s3 ls "$s3Dir"  
