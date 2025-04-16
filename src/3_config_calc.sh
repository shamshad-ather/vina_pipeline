#!/bin/bash
# src/03_generate_configs.sh

# Ensure the script exits if any command fails
set -e

# --- Determine location of this script and the python script ---
# Get the directory where this script itself resides (src)
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) 
# Define the python script name
PYTHON_SCRIPT_NAME="calc_box_center.py"
# Construct the full path to the python script within the src directory
PYTHON_SCRIPT_PATH="${SCRIPT_DIR}/${PYTHON_SCRIPT_NAME}" 

# --- Define input/output directories relative to the parent folder ---
# Directory containing prepared receptor PDBQT files (expected in parent folder)
RECEPTOR_PDBQT_DIR="receptors_pdbqt" 
# Output directory for config files (expected in parent folder)
CONFIG_DIR="configs"
# Buffer size for the box
BUFFER_SIZE=10.0 

# --- Create output directory in the parent folder ---
# Note: Working directory is the parent folder when called by master.sh
mkdir -p "$CONFIG_DIR"

echo "--- Starting Config File Generation ---"

# --- Check if the Python script exists at its expected location ---
if [ ! -f "$PYTHON_SCRIPT_PATH" ]; then
    echo "Error: Python script '$PYTHON_SCRIPT_NAME' not found at '$PYTHON_SCRIPT_PATH'."
    exit 1
fi

# --- Check if input directory exists in the parent folder ---
if [ ! -d "$RECEPTOR_PDBQT_DIR" ] || [ -z "$(find "$RECEPTOR_PDBQT_DIR" -maxdepth 1 -name '*.pdbqt' -print -quit)" ]; then
    echo "Error: Input directory '$RECEPTOR_PDBQT_DIR' not found or contains no .pdbqt files (relative to execution directory)."
    exit 1
fi

# --- Loop through receptors and call Python script ---
# Loop through all prepared receptor PDBQT files in the parent folder's receptor dir
for receptor_pdbqt in "$RECEPTOR_PDBQT_DIR"/*.pdbqt; do
  echo "Generating config for: $(basename "$receptor_pdbqt")"
  # Execute python using the full path to the script in src/
  # Pass arguments (receptor path, buffer size) which are relative to the parent folder
  python "$PYTHON_SCRIPT_PATH" "$receptor_pdbqt" "$BUFFER_SIZE" 
  if [ $? -ne 0 ]; then
      echo "  Error generating config for $(basename "$receptor_pdbqt"). Stopping."
      exit 1 # Stop if config generation fails
  fi
done

echo "--- Config File Generation Finished ---"