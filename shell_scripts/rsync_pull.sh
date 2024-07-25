#!/bin/bash

# Define source directories
GOOGLE_DRIVE_INPUT_DATA_DIR="/Users/alexander.prokopyszyn/Library/CloudStorage/GoogleDrive-githubburner002@gmail.com/My Drive/locust_2024_paper_code/input_data"
GOOGLE_DRIVE_PLOTS_DIR="/Users/alexander.prokopyszyn/Library/CloudStorage/GoogleDrive-githubburner002@gmail.com/My Drive/locust_2024_paper_code/plots"
ONEDRIVE_OUTPUT_DATA_DIR="/Users/alexander.prokopyszyn/Library/CloudStorage/OneDrive-ScienceandTechnologyFacilitiesCouncil/GitHub/locust_2024_paper_code_output_data/output_data"

# Define destination directories
INPUT_DATA_DIR="input_data"
PLOTS_DIR="plots"
OUTPUT_DATA_DIR="output_data"

# Rsync options
RSYNC_OPTIONS="-avh --progress"

# Sync Google Drive input_data directory to local
echo "Syncing ${GOOGLE_DRIVE_INPUT_DATA_DIR} to ${INPUT_DATA_DIR}..."
rsync $RSYNC_OPTIONS "$GOOGLE_DRIVE_INPUT_DATA_DIR/" $INPUT_DATA_DIR/

# Sync Google Drive plots directory to local
echo "Syncing ${GOOGLE_DRIVE_PLOTS_DIR} to ${PLOTS_DIR}..."
rsync $RSYNC_OPTIONS "$GOOGLE_DRIVE_PLOTS_DIR/" $PLOTS_DIR/

# Sync OneDrive output_data directory to local
echo "Syncing ${ONEDRIVE_OUTPUT_DATA_DIR} to ${OUTPUT_DATA_DIR}..."
rsync $RSYNC_OPTIONS "$ONEDRIVE_OUTPUT_DATA_DIR/" $OUTPUT_DATA_DIR/

echo "Syncing completed."
