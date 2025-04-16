#!/bin/bash

# Ensure the script exits if any command fails
set -e

# Directories
RECEPTOR_DIR="receptors_pdbqt"
LIGAND_DIR="ligands_pdbqt"
CONFIG_DIR="configs"
RESULTS_DIR="docking_results"
VINA_EXECUTABLE="vina" # Make sure 'vina' is in your PATH or provide full path

# Vina parameters
EXHAUSTIVENESS=8
NUM_MODES=9 # Often useful to save more than just the best pose

# Create base output directory
mkdir -p "$RESULTS_DIR"

echo "--- Starting AutoDock Vina Docking ---"

# Check directories and Vina executable
if [ ! -d "$RECEPTOR_DIR" ]; then echo "Error: Receptor directory '$RECEPTOR_DIR' not found."; exit 1; fi
if [ ! -d "$LIGAND_DIR" ]; then echo "Error: Ligand directory '$LIGAND_DIR' not found."; exit 1; fi
if [ ! -d "$CONFIG_DIR" ]; then echo "Error: Config directory '$CONFIG_DIR' not found."; exit 1; fi
if ! command -v $VINA_EXECUTABLE &> /dev/null; then echo "Error: Vina executable '$VINA_EXECUTABLE' not found in PATH."; exit 1; fi

# Loop through all receptors
for receptor in "$RECEPTOR_DIR"/*.pdbqt; do
  receptor_name=$(basename "${receptor%.pdbqt}")
  config_file="${CONFIG_DIR}/${receptor_name}.txt"

  # Check if config file exists
  if [ ! -f "$config_file" ]; then
      echo "Warning: Config file '$config_file' not found for receptor '$receptor_name'. Skipping this receptor."
      continue
  fi

  echo "--- Processing Receptor: $receptor_name ---"

  # Loop through all ligands
  for ligand in "$LIGAND_DIR"/*.pdbqt; do
    ligand_name=$(basename "${ligand%.pdbqt}")
    
    # Create a specific output directory for this pair
    # Using underscore separation for clarity
    output_dir="${RESULTS_DIR}/${receptor_name}_${ligand_name}" 
    mkdir -p "$output_dir"

    # Define output file paths
    output_pdbqt="${output_dir}/docked_poses.pdbqt"
    output_log="${output_dir}/vina_output.log"

    echo "  Docking $ligand_name to $receptor_name..."

    # Run AutoDock Vina
    $VINA_EXECUTABLE \
      --receptor "$receptor" \
      --ligand "$ligand" \
      --config "$config_file" \
      --exhaustiveness "$EXHAUSTIVENESS" \
      --num_modes "$NUM_MODES" --out "$output_pdbqt" \
      2>&1 | tee "$output_log" # Use tee to capture stdout/stderr to file and terminal

    # Check Vina exit status
    if [ $? -ne 0 ]; then
        echo "  Error: Vina failed for $ligand_name docking to $receptor_name. Check log: $output_log"
        # Decide whether to continue or exit (current: continue)
        # exit 1 # Uncomment to stop the whole script on a single Vina failure
    else
        echo "  Docking successful. Results saved to $output_dir"
    fi
    
  done # End ligand loop
done # End receptor loop

echo "--- AutoDock Vina Docking Finished ---"