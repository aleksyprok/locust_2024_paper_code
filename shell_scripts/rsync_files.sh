#!/bin/bash

# Define source directories
INPUT_DATA_DIR="input_data"
PLOTS_DIR="plots"
OUTPUT_DATA_DIR="output_data"

# Define destination directories
GOOGLE_DRIVE_INPUT_DATA_DIR="/Users/alexander.prokopyszyn/Library/CloudStorage/GoogleDrive-githubburner002@gmail.com/My Drive/locust_2024_paper_code/input_data"
GOOGLE_DRIVE_PLOTS_DIR="/Users/alexander.prokopyszyn/Library/CloudStorage/GoogleDrive-githubburner002@gmail.com/My Drive/locust_2024_paper_code/plots"
ONEDRIVE_OUTPUT_DATA_DIR="/Users/alexander.prokopyszyn/Library/CloudStorage/OneDrive-ScienceandTechnologyFacilitiesCouncil/GitHub/locust_2024_paper_code_output_data/output_data"

# Rsync options
RSYNC_OPTIONS="-avh --progress"

# Sync input_data directory to Google Drive
echo "Syncing ${INPUT_DATA_DIR} to ${GOOGLE_DRIVE_INPUT_DATA_DIR}..."
rsync $RSYNC_OPTIONS $INPUT_DATA_DIR/ "$GOOGLE_DRIVE_INPUT_DATA_DIR/"

# Sync plots directory to Google Drive
echo "Syncing ${PLOTS_DIR} to ${GOOGLE_DRIVE_PLOTS_DIR}..."
rsync $RSYNC_OPTIONS $PLOTS_DIR/ "$GOOGLE_DRIVE_PLOTS_DIR/"

# Sync output_data directory to OneDrive
echo "Syncing ${OUTPUT_DATA_DIR} to ${ONEDRIVE_OUTPUT_DATA_DIR}..."
rsync $RSYNC_OPTIONS $OUTPUT_DATA_DIR/ "$ONEDRIVE_OUTPUT_DATA_DIR/"

echo "Syncing completed."
