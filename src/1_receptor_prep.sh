#!/bin/bash

# Ensure the script exits if any command fails
set -e

# Input directory for original receptor PDB files
RECEPTOR_PDB_DIR="receptors"
# Output directory for prepared receptor PDBQT files
RECEPTOR_PDBQT_DIR="receptors_pdbqt" 

# Create output directory if it doesn't exist
mkdir -p "$RECEPTOR_PDBQT_DIR"

echo "--- Starting Receptor Preparation ---"

# Check if input directory exists and has PDB files
if [ ! -d "$RECEPTOR_PDB_DIR" ] || [ -z "$(find "$RECEPTOR_PDB_DIR" -maxdepth 1 -name '*.pdb' -print -quit)" ]; then
    echo "Error: Input directory '$RECEPTOR_PDB_DIR' not found or contains no .pdb files."
    exit 1
fi

# Loop through all PDB files in the input directory
for receptor_pdb in "$RECEPTOR_PDB_DIR"/*.pdb; do
    # Get the base name without extension
    base_name=$(basename "${receptor_pdb%.pdb}")
    
    # Define intermediate and final file names
    intermediate_pdb="${RECEPTOR_PDBQT_DIR}/${base_name}_intermediate.pdb"
    receptor_pdbqt="${RECEPTOR_PDBQT_DIR}/${base_name}.pdbqt"
    
    echo "Processing receptor: $base_name"
    
    # Step 1: Add hydrogens at pH 7.4 and remove heteroatoms (water, ligands etc.) -> Intermediate PDB
    # Use -d to delete input hydrogens before adding new ones
    echo "  Adding hydrogens and cleaning..."
    obabel "$receptor_pdb" -O "$intermediate_pdb" -xr -d -p 7.4
    if [ $? -ne 0 ]; then
        echo "  Error during hydrogen addition/cleaning for $base_name. Skipping."
        rm -f "$intermediate_pdb" # Clean up intermediate file on error
        continue # Skip to the next receptor
    fi
    
    # Step 2: Convert intermediate PDB to PDBQT, assigning partial charges
    echo "  Converting to PDBQT and assigning charges..."
    obabel "$intermediate_pdb" -O "$receptor_pdbqt" -xr --partialcharge eem 
    if [ $? -ne 0 ]; then
        echo "  Error during PDBQT conversion for $base_name. Skipping."
        rm -f "$intermediate_pdb" "$receptor_pdbqt" # Clean up files on error
        continue # Skip to the next receptor
    fi

    # Remove the intermediate files
    rm "$intermediate_pdb"
    
    echo "  Successfully prepared: $receptor_pdbqt"
done

echo "--- Receptor Preparation Finished ---"