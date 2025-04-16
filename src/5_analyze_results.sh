#!/bin/bash
# src/05_extract_analyze_results.sh

# Ensure the script exits if any command fails
# set -e # Commenting out set -e to allow processing remaining files even if one fails

# --- Determine location of this script and the python script ---
# Get the directory where this script itself resides (src)
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Define the python script name
PYTHON_ANALYSIS_SCRIPT_NAME="calculate_docking_metrics.py"
# Construct the full path to the python script within the src directory
PYTHON_ANALYSIS_SCRIPT_PATH="${SCRIPT_DIR}/${PYTHON_ANALYSIS_SCRIPT_NAME}" 

# --- Define input/output directories/files relative to the parent folder ---
# Directory containing docking results (subdirectories per pair)
RESULTS_DIR="docking_results"
# Output CSV file (in parent folder)
OUTPUT_CSV="docking_summary_metrics.csv"
# Directories to search for ligand files (relative to parent folder)
LIGAND_SOURCE_DIRS="ligands ligands_pdbqt" 

echo "--- Starting Results Extraction and Analysis ---"

# --- Check prerequisites ---
# Check if results directory exists in parent folder
if [ ! -d "$RESULTS_DIR" ]; then
    echo "Error: Docking results directory '$RESULTS_DIR' not found (relative to execution dir)."
    exit 1
fi
# Check if the Python script exists at its expected location in src/
if [ ! -f "$PYTHON_ANALYSIS_SCRIPT_PATH" ]; then
    echo "Error: Python analysis script '$PYTHON_ANALYSIS_SCRIPT_NAME' not found at '$PYTHON_ANALYSIS_SCRIPT_PATH'."
    exit 1
fi

# --- Prepare Output CSV ---
# Write header row to CSV in parent folder (overwrite existing file)
echo "Receptor,Ligand,Affinity_kcal_mol,pKi,MW,LogP,TPSA,NHA,NRB,HBD,HBA,SASA_A2,QED,LE,LLE,SILE_N,SILE_SASA" > "$OUTPUT_CSV"

# --- Process Logs ---
# Find all Vina log files within the results directory (relative to parent)
# Use find ... -print0 | while ... read -d $'\0' for safer handling of filenames
find "$RESULTS_DIR" -name "vina_output.log" -print0 | while IFS= read -r -d $'\0' logfile; do
    echo "Processing log: $logfile"
    
    # Run the Python script using its full path from src/
    # Pass arguments (log file path, ligand source dirs) which are relative to the parent folder
    # Append the output (a single CSV row) directly to the output CSV file in parent folder
    python "$PYTHON_ANALYSIS_SCRIPT_PATH" "$logfile" --ligand_dirs $LIGAND_SOURCE_DIRS >> "$OUTPUT_CSV"
    
    # Optional: Check python script exit status
    if [ $? -ne 0 ]; then
        echo "  Warning: Python script failed for $logfile. Check stderr output. Continuing..."
        # Write a placeholder or error marker to the CSV if needed
        # echo "$(basename $logfile),Error,,Processing Failed,,,,,,,,,,,,,," >> "$OUTPUT_CSV" 
    fi
done

echo "--- Results Extraction and Analysis Finished ---"
echo "Summary saved to: $OUTPUT_CSV"

# Optional: Display first few lines of the result file
echo "--- First 10 lines of $OUTPUT_CSV ---"
head -n 10 "$OUTPUT_CSV"
echo "--------------------------------------"