#  bash script to modify filenames
#
path="/media/jon/Seagate_5TB/OPRAM_results/gisdata" # the path of the directory where the files and directories that need to be copied are located

for d in "$path"/*; do
  echo $d
  
  if [ -d "$d" ]; then 
    # Update dir name
    newd=`echo "$d" | sed -r 's/thesh/thresh/g'`

    if [ "$d" != "$newd" ]; then
      #echo $newd
      mv "$d" "$newd"
    fi
  
    # Look at files within the directory
    for entry in "$newd"/*; do
      newentry=`echo "$entry" | sed -r 's/thesh/thresh/g'`

      if [ "$entry" != "$newentry" ]; then
        #echo $newentry
        mv "$entry" "$newentry"
      fi
    done
  fi
done
