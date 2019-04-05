#!/bin/bash

for entry in `less bam_files`
do
	echo ${entry}
	nohup python cbam_deepest_coverage.py -d /dataset/MBIE_genomics4production/scratch/rudi_gbsathon/sheep/PstI-MspI_v2/bwa_Oarv3.1/2_lanes/ -i $entry -c chromosomes.txt -o ${entry}_coverage_2lanes.txt -l 2_lane &
done
