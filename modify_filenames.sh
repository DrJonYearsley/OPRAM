#  bash script to modify filenames
#
path="/media/jon/Seagate_5TB/OPRAM_results/data" # the path of the directory where the files and directories that need to be copied are located

for d in "$path"/*; do
  echo $d
  
  # Look at files within the directory
  for entry in "$d"/*; do
    newentry=`echo "$entry" | sed -r 's/thesh/thresh/g'`

    if [ "$entry" != "$newentry" ]; then
    #  echo $newentry
      mv "$entry" "$newentry"
    fi
  done
done
