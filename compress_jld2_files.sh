#!/bin/bash

# Shell script to package and compress jld2 files
# from the OPRAM degree day models

cd /media/jon/Seagate_5TB/OPRAM_results


for f in base*_thresh*
do
        echo $f
        # Produce a compressed tar file wil maximum compression
        tar -cvf ${f}_jld2.tar.gz -z  ${f}/${f}*.jld2

        rm -r ${f}/${f}*.jld2
        rm -r ${f}
done

