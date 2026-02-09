# OPRAM

The primary goal of this project is to provide a practical visual tool that communicates scientific evidence to stakeholders involved in Ireland's plant biosecurity policy and pest/pathogen-risk management. This supports the objectives of Ireland's Plant Health and Biosecurity Strategy 2020-2025 (DAFM 2019): "To minimise the threat posed to plants by the potential introduction and establishment of plant pests and diseases". 

More details at https://DrJonYearsley.github.io/pesttool.html




==============================

## File Handling Scripts


| Filename | Description |
| -----------| ----------------|
| `run_opram.sh` | Main script to run the OPRAM model (requires Julia and R to be installed) |
| `compress_jld2_files.sh`  | Compress OPRAM model results (i.e. JLD2 files) |
| `aws_upload_tar.sh`  | Extract CSV files from compressed tar file and upload them to OPRAM Amazon Web Services bucket (extracted files are deleted after upload) |
| `aws_upload.sh`  | Upload multiple directories of CSV files to OPRAM Amazon Web Services bucket |
| `aws_upload_gis.sh`  | Upload geoTIFF files to OPRAM Amazon Web Services bucket |
| `aws_delete.sh` |  Delete specific files from OPRAM Amazon Web Services bucket |
| `aws_upload_1dir.sh` | Upload a single directory of CSV files to OPRAM Amazon Web Services bucket |
| `modify_filenames.sh` | Bulk-modify filenames |

