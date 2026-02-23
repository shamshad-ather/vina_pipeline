import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import os
import concurrent.futures
from multiprocessing import cpu_count

# ================= CONFIGURATION =================
INPUT_FILE = 'input_smiles.xlsx'  # Your Excel file
SHEET_NAME = 0                    # Sheet name or index (0 is first sheet)
CID_COL    = 'CID'                # Column name for the ID (used for filename)
SMILES_COL = 'SMILES'             # Column name for the Structure
OUTPUT_DIR = 'ligands'     # Folder to save the .sdf files
NUM_CONFS  = 50                   # Number of conformers to generate
# =================================================

def generate_best_conformer(data):
    """
    Worker function to process a single molecule.
    Args:
        data (tuple): (cid, smiles)
    Returns:
        str: Status message
    """
    cid, smiles = data
    
    # 1. Sanitize Input
    if not smiles or pd.isna(smiles):
        return f"[{cid}] Skipped: Empty SMILES"
    
    # Clean CID for filename usage
    safe_cid = str(cid).strip().replace('/', '_').replace('\\', '_')
    output_path = os.path.join(OUTPUT_DIR, f"{safe_cid}.sdf")

    # 2. Create Molecule
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return f"[{cid}] Failed: Invalid SMILES"

    try:
        mol = Chem.AddHs(mol) # Add Hydrogens (Crucial for 3D)
        
        # 3. Embed Conformers (ETKDGv3 is the state-of-the-art in RDKit)
        params = AllChem.ETKDGv3()
        params.useSmallRingTorsions = True
        
        # Embed 50 random conformers
        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=NUM_CONFS, params=params)
        
        if not conf_ids:
            return f"[{cid}] Failed: Embedding Error (Structure might be too complex)"

        # 4. Optimize All Conformers (MMFF94s)
        # MMFF94s is generally best for drug-like organic molecules
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
        
        if props is None:
            # Fallback to UFF if MMFF fails
            force_field_type = "UFF"
            res = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=200)
        else:
            force_field_type = "MMFF94s"
            res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=200, mmffVariant='MMFF94s')


        # 5. Find the Lowest Energy Conformer
        # res is a list of tuples: (is_converged, energy)
        lowest_energy = float('inf')
        best_conf_id = -1

        for i, (not_converged, energy) in enumerate(res):
            if energy < lowest_energy:
                lowest_energy = energy
                best_conf_id = conf_ids[i]

        if best_conf_id == -1:
            return f"[{cid}] Failed: Optimization could not converge"

        # 6. Save ONLY the best conformer
        # Create a new molecule object with only the best conformer
        best_mol = Chem.Mol(mol)
        best_mol.RemoveAllConformers()
        best_mol.AddConformer(mol.GetConformer(best_conf_id))
        
        # Set the title of the molecule to the CID
        best_mol.SetProp("_Name", str(cid))
        best_mol.SetProp("Energy_Min", f"{lowest_energy:.4f}")
        best_mol.SetProp("Force_Field", force_field_type)

        writer = Chem.SDWriter(output_path)
        writer.write(best_mol)
        writer.close()

        return None # Success (returning None means no error)

    except Exception as e:
        return f"[{cid}] Error: {str(e)}"

def main():
    # Setup Output Directory
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        print(f"Created directory: {OUTPUT_DIR}")

    # Load Excel
    print(f"Reading {INPUT_FILE}...")
    try:
        df = pd.read_excel(INPUT_FILE, sheet_name=SHEET_NAME)
    except Exception as e:
        print(f"CRITICAL ERROR: Could not read Excel file.\n{e}")
        return

    # Validate Columns
    if CID_COL not in df.columns or SMILES_COL not in df.columns:
        print(f"CRITICAL ERROR: Columns '{CID_COL}' or '{SMILES_COL}' not found.")
        print(f"Found columns: {list(df.columns)}")
        return

    # Prepare Data List
    tasks = []
    for _, row in df.iterrows():
        tasks.append((row[CID_COL], row[SMILES_COL]))

    total_mols = len(tasks)
    workers = cpu_count() # Use all CPU cores
    
    print(f"Starting processing of {total_mols} molecules.")
    print(f"Engine: RDKit ETKDGv3 + MMFF94s")
    print(f"Parallel Workers: {workers}")
    print("-" * 60)

    # Execute in Parallel
    success_count = 0
    errors = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        # Submit all tasks
        future_to_cid = {executor.submit(generate_best_conformer, task): task[0] for task in tasks}
        
        # Monitor progress
        completed = 0
        for future in concurrent.futures.as_completed(future_to_cid):
            completed += 1
            result = future.result()
            
            # Simple progress bar
            percent = (completed / total_mols) * 100
            print(f"\rProgress: [{completed}/{total_mols}] {percent:.1f}%", end="")

            if result is None:
                success_count += 1
            else:
                errors.append(result)

    print("\n" + "-" * 60)
    print("PROCESSING COMPLETE")
    print(f"Successfully Generated: {success_count}")
    print(f"Failed: {len(errors)}")
    
    if errors:
        print("\n--- Error Log ---")
        for err in errors[:10]: # Print first 10 errors
            print(err)
        if len(errors) > 10:
            print(f"...and {len(errors)-10} more.")

if __name__ == "__main__":
    main()