#!/bin/bash

INPUT_DATA_DIR="input_data"
CSD3_DIRECTORY="ir-prok1@login-gpu.hpc.cam.ac.uk:/rds/project/rds-aSo1XX0UOlw/ir-prok1/GitHub/locust_2024_paper_code/input_data"
RSYNC_OPTIONS="-avHhe ssh --progress"

echo "Syncing ${INPUT_DATA_DIR} to ${CSD3_DIRECTORY}..."
rsync $RSYNC_OPTIONS "$INPUT_DATA_DIR/" "$CSD3_DIRECTORY/"