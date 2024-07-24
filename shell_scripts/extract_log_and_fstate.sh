#!/bin/bash

SRC_DIR="../output_data_raw"
DEST_DIR="../output_data"

# Create the destination directory
mkdir -p "$DEST_DIR"

# Function to process LOG*.out files
process_log_file() {
    local file_path="$1"
    echo "Processing LOG file: $file_path"
    sed -i '/:LOCUST      : OpenMP MAIN loop/,/:locust_stat : GPU processing time/d' "$file_path"
    echo "Finished processing LOG file: $file_path"
}

# Copy relevant files and process LOG*.out files
find "$SRC_DIR" -type f \( -name 'FINAL_STATE*.dat' -o -name 'LOG*.out' \) | while read -r file; do
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
