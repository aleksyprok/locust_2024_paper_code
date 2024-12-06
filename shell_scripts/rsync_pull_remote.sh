#!/bin/bash

# REMOTE_DIRECTORY="ir-prok1@login-gpu.hpc.cam.ac.uk:/rds/project/rds-aSo1XX0UOlw/ir-prok1/GitHub/locust_2024_paper_code/input_data"
REMOTE_DIRECTORY="aprokopy@login.leonardo.cineca.it:/leonardo_scratch/large/userexternal/aprokopy/locust.STEP/OutputFiles/lrdn0514"
OUTPUT_DATA_DIR="output_data"
RSYNC_OPTIONS="-avHhe ssh --progress"

echo "Syncing ${REMOTE_DIRECTORY} to ${OUTPUT_DATA_DIR}..."
rsync $RSYNC_OPTIONS "$REMOTE_DIRECTORY" "$OUTPUT_DATA_DIR/."