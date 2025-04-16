#!/bin/bash
# src/02_prepare_ligands.sh
# Prepares ligands, checking file extensions to determine input format.

# Ensure the script exits if any command fails
set -e

# --- Define input/output directories relative to the parent folder ---
LIGAND_INPUT_DIR="ligands"
LIGAND_PDBQT_DIR="ligands_pdbqt"

# --- Create output directory in the parent folder ---
mkdir -p "$LIGAND_PDBQT_DIR"

echo "--- Starting Ligand Preparation (with format check) ---"

# --- Check if input directory exists in the parent folder ---
if [ ! -d "$LIGAND_INPUT_DIR" ]; then
    echo "Error: Input ligand directory '$LIGAND_INPUT_DIR' not found (relative to execution directory)."
    exit 1
fi

# Check if input directory is empty
if [ -z "$(ls -A "$LIGAND_INPUT_DIR")" ]; then
    echo "Warning: Input ligand directory '$LIGAND_INPUT_DIR' is empty. No ligands to prepare."
    echo "--- Ligand Preparation Finished (No files processed) ---"
    exit 0
fi

# --- Loop through all files in the input directory ---
shopt -s extglob # Enable extended globbing for case-insensitive matching

for ligand_file in "$LIGAND_INPUT_DIR"/*; do
    
    # Skip if not a regular file (e.g., directory, symlink)
    if [ ! -f "$ligand_file" ]; then
        # Optional: Warn about non-files found
        # echo "  Skipping non-file item: $(basename "$ligand_file")"
        continue
    fi

    # Get the base name and extension
    base_name=$(basename "$ligand_file")
    extension="${base_name##*.}"
    # Get base name without extension for output file naming
    base_name_no_ext="${base_name%.*}" 
    
    # Convert extension to lowercase for case-insensitive check
    extension_lower=$(echo "$extension" | tr '[:upper:]' '[:lower:]')

    # Define output PDBQT file path (in parent folder's output dir)
    ligand_pdbqt="${LIGAND_PDBQT_DIR}/${base_name_no_ext}.pdbqt"
    
    # Variable to hold the input format flag for obabel
    input_format_flag=""

    echo "Processing file: $base_name"

    # --- Determine input format based on extension ---
    case "$extension_lower" in
        sdf)
            input_format_flag="-isdf"
            ;;
        mol2)
            input_format_flag="-imol2"
            ;;
        mol)
            # Treat .mol like .sdf
            input_format_flag="-imol" 
            ;;
        pdb)
            input_format_flag="-ipdb"
            ;;
        # Add more formats if needed (e.g., smi for SMILES)
        # smi | smiles)
        #    input_format_flag="-ismi"
        #    # SMILES often requires generating 3D coordinates
        #    gen3d_flag="--gen3d -best" 
        #    ;;
        *)
            echo "  Warning: Skipping '$base_name'. Unrecognized or unsupported extension '.$extension_lower'."
            continue # Skip to the next file in the loop
            ;;
    esac
    
    # Optional: Reset gen3d_flag if not set by case (e.g. for SMILES)
    : ${gen3d_flag:=} # Sets to empty string if undefined

    # --- Run Open Babel conversion ---
    echo "  Detected format: .$extension_lower. Converting to PDBQT..."
    obabel \
      "$input_format_flag" "$ligand_file" \
      -opdbqt \
      -O "$ligand_pdbqt" \
      -h -p 7.4 --partialcharge gasteiger \
      $gen3d_flag \
      # Optional: Add -m if your source files might contain MULTIPLE molecules 
      # and you want to process them all individually. Note: This will change 
      # output naming (e.g., base_1.pdbqt, base_2.pdbqt).
      # -m 

    # Check obabel exit status
    if [ $? -ne 0 ]; then
        echo "  Error during PDBQT conversion for '$base_name'. Check Open Babel output. Skipping."
        rm -f "$ligand_pdbqt" # Clean up potentially incomplete output file
        continue # Skip to the next ligand
    else
        echo "  Successfully prepared: $ligand_pdbqt"
    fi
    
    # Clear gen3d_flag for next iteration if it was set
    gen3d_flag="" 

done

shopt -u extglob # Disable extended globbing

echo "--- Ligand Preparation Finished ---"