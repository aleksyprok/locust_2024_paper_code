#!/bin/bash

CSD3_DIRECTORY="full_3d_field_spr_045_14_and_16_processed"
CSD3_PATH="ir-prok1@login-gpu.hpc.cam.ac.uk:/rds/project/rds-aSo1XX0UOlw/ir-prok1/locust.STEP/OutputFiles/"$CSD3_DIRECTORY
OUTPUT_DATA_DIR="output_data"
RSYNC_OPTIONS="-avHhe ssh --progress"

echo "Syncing ${CSD3_PATH} to ${OUTPUT_DATA_DIR}..."
rsync $RSYNC_OPTIONS "$CSD3_PATH" "$OUTPUT_DATA_DIR/."