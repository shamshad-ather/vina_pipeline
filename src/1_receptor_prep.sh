#!/bin/bash

# Ensure the script exits if any command fails
set -e

# --- CONFIGURATION ---
# Input/Output Directories
RECEPTOR_PDB_DIR="receptors"
RECEPTOR_PDBQT_DIR="receptors_pdbqt"

# COFACTOR WHITELIST
# Add the 3-letter PDB codes of cofactors you want to KEEP.
# Separate them with the pipe symbol '|' inside the quotes.
ALLOWED_COFACTORS="HEM|NAD|FAD|NAP|ADP|ATP|MG|ZN|MN|CA|FE"

# Create output directory if it doesn't exist
mkdir -p "$RECEPTOR_PDBQT_DIR"

echo "--- Starting Receptor Preparation ---"
echo "--- Keeping Protein and Cofactors: ${ALLOWED_COFACTORS//|/, } ---"

# Check if input directory exists and has PDB files
if [ ! -d "$RECEPTOR_PDB_DIR" ] || [ -z "$(find "$RECEPTOR_PDB_DIR" -maxdepth 1 -name '*.pdb' -print -quit)" ]; then
    echo "Error: Input directory '$RECEPTOR_PDB_DIR' not found or contains no .pdb files."
    exit 1
fi

# Loop through all PDB files in the input directory
for receptor_pdb in "$RECEPTOR_PDB_DIR"/*.pdb; do
    # Get the base name without extension
    base_name=$(basename "${receptor_pdb%.pdb}")
    
    # Define file names
    clean_pdb="${RECEPTOR_PDBQT_DIR}/${base_name}_clean.pdb"
    receptor_pdbqt="${RECEPTOR_PDBQT_DIR}/${base_name}.pdbqt"
    
    echo "Processing receptor: $base_name"
    
    # Step 1: Filter the PDB File (The 'Cleaning' Step)
    # Logic: Keep lines starting with 'ATOM' (Protein/DNA) OR 'TER' (Chain breaks)
    #        OR 'HETATM' lines that match our ALLOWED_COFACTORS list.
    #        Everything else (Water 'HOH', Ligands 'LIG', unknown ions) is dropped.
    echo "  Filtering out water and unwanted ligands..."
    grep -E "^ATOM|^TER|^HETATM.*($ALLOWED_COFACTORS)" "$receptor_pdb" > "$clean_pdb"

    # Step 2: Process with Open Babel
    # -d: Delete existing hydrogens (to ensure clean slate)
    # -p 7.4: Add hydrogens appropriate for pH 7.4
    # -xr: Output as 'rigid' molecule (prevents torsion tree generation on the receptor)
    # --partialcharge eem: Calculate charges
    echo "  Adding hydrogens, charges, and converting to PDBQT..."
    obabel "$clean_pdb" -O "$receptor_pdbqt" -xr -p 7.4 --partialcharge gasteiger    
    # Error checking for Obabel
    if [ $? -ne 0 ]; then
        echo "  Error during Open Babel processing for $base_name. Skipping."
        rm -f "$clean_pdb" "$receptor_pdbqt"
        continue
    fi

    # Delete the intermediate cleaned PDB
    rm "$clean_pdb"
    
    echo "  Successfully prepared: $receptor_pdbqt"
done

echo "--- Receptor Preparation Finished ---"
