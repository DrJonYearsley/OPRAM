# Script to download Met Eirean's gridded data
#
# Downloads Tmax and Tmin daily data
#
# Jon Yearsley
# Aug 2026
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++

# Array of years to download
years=("2025")

output_path="/Volumes/Macintosh HD/Users/jon/Google Drive/My Drive/Projects/DAFM_OPRAM/Data/Irish Climate Data"

# URL of Met Eireann server containing the data
data_url="https://clidata.met.ie/cli/grids_daily/latest/"
months=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")


# Import max daily temp
echo "Importing daily max temperature"
for y in "${years[@]}"; do
  for m in "${months[@]}"; do
      echo "${data_url}/max/IRL_DLY_TX_${y}${m}_grid_IDW_09Z.csv.gz"
      curl -fsSL "${data_url}/max/IRL_DLY_TX_${y}${m}_grid_IDW_09Z.csv.gz" | gunzip > "${output_path}/maxtemp_grids/IRL_DLY_TX_${y}${m}_grid_IDW_09Z.csv"
  done
done


# Import min daily temp
echo "Importing daily min temperature"
for y in "${years[@]}"; do
  for m in "${months[@]}"; do
      echo "${data_url}/min/IRL_DLY_TN_${y}${m}_grid_IDW_09Z.csv.gz"
      curl -fsSL "${data_url}/min/IRL_DLY_TN_${y}${m}_grid_IDW_09Z.csv.gz" | gunzip > "${output_path}/mintemp_grids/IRL_DLY_TN_${y}${m}_grid_IDW_09Z.csv"
  done
done


