#!/bin/bash

SRC_DIR="/rds/project/rds-aSo1XX0UOlw/ir-prok1/locust.STEP/OutputFiles/scan_poincare_toggle_i3dr"
DEST_DIR="/rds/project/rds-aSo1XX0UOlw/ir-prok1/locust.STEP/OutputFiles/scan_poincare_toggle_i3dr_processed"

# Create the destination directory
mkdir -p "$DEST_DIR"

# Function to process LOG*.out files
process_log_file() {
    local file_path="$1"
    echo "Processing LOG file: $file_path"
    sed -i '/:LOCUST      : OpenMP MAIN loop/,/:LOCUST      : End of OpenMP MAIN loop/d' "$file_path"
    echo "Finished processing LOG file: $file_path"
}

# Copy relevant files and process LOG*.out files
find "$SRC_DIR" -type f \( -name 'LOG*.out' -o -name 'Poincare_map_*_formatted.dat' \) | while read -r file; do
    # Construct the destination file path
    dest_file="$DEST_DIR/${file#$SRC_DIR/}"
    dest_dir="$(dirname "$dest_file")"

    # Create destination directory and copy the file
    mkdir -p "$dest_dir"
    cp "$file" "$dest_file"

    # Process LOG*.out files
    if [[ $file == *"LOG"*".out" ]]; then
        process_log_file "$dest_file"
    fi
done

echo "Script execution completed."
