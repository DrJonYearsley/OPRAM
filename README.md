# OPRAM

The primary goal of this project is to provide a practical visual tool that communicates scientific evidence to stakeholders involved in Ireland's plant biosecurity policy and pest/pathogen-risk management. This supports the objectives of Ireland's Plant Health and Biosecurity Strategy 2020-2025 (DAFM 2019): "To minimise the threat posed to plants by the potential introduction and establishment of plant pests and diseases". 

More details at https://DrJonYearsley.github.io/pesttool.html




==============================

## File Handling Scripts


| Filename | Description |
| -----------| ----------------|
| `run_opram.sh` | Main script to run the OPRAM model and produce CSV and geotiff output (requires Julia and R to be installed) |
| `run_main_test_UK_NI.sh` | Run the OPRAM model for Northern Ireland on either UK Met Office or Met Eireann data |
| `compress_jld2_files.sh`  | Compress OPRAM model results (i.e. JLD2 files) |
| `modify_filenames.sh` | Bulk-modify filenames |

