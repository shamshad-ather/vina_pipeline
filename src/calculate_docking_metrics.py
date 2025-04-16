import sys
import os
import re
import math
import argparse
import warnings
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, rdMolDescriptors, Crippen

# Constants for pKi calculation
R = 1.987204259e-3  # Ideal gas constant in kcal/(mol·K)
T = 298.15          # Standard temperature in Kelvin (approx. 25 C)
RT = R * T

def calculate_pKi(delta_G_kcal_mol):
    """Converts Vina Affinity (ΔG in kcal/mol) to pKi."""
    if delta_G_kcal_mol is None:
        return None
    try:
        # Ki = exp(ΔG / RT) - ΔG must be in kcal/mol here
        delta_G_kcal_mol = float(delta_G_kcal_mol)
        ki = math.exp(delta_G_kcal_mol / RT)
        # pKi = -log10(Ki)
        # Avoid math domain error if Ki is zero or negative (shouldn't happen with exp)
        if ki <= 0: 
             warnings.warn(f"Calculated Ki is non-positive ({ki}) from Affinity {delta_G_kcal_mol}. Cannot calculate pKi.")
             return None
        pki = -math.log10(ki)
        return pki
    except (ValueError, TypeError, OverflowError) as e:
        warnings.warn(f"Could not calculate pKi for affinity {delta_G_kcal_mol}: {e}")
        return None

def calculate_ligand_metrics(mol):
    """Calculates various metrics for an RDKit molecule object."""
    if mol is None:
        return {k: None for k in ["MW", "LogP", "TPSA", "NHA", "NRB", "HBD", "HBA", "SASA", "QED"]}
        
    try:
        Chem.SanitizeMol(mol) # Ensure sanitization
        metrics = {}
        metrics["MW"] = Descriptors.MolWt(mol)
        metrics["LogP"] = Crippen.MolLogP(mol)
        metrics["TPSA"] = Descriptors.TPSA(mol)
        metrics["NHA"] = mol.GetNumHeavyAtoms()
        metrics["NRB"] = Descriptors.NumRotatableBonds(mol)
        metrics["HBD"] = Descriptors.NumHDonors(mol)
        metrics["HBA"] = Descriptors.NumHAcceptors(mol)
        # SASA calculation requires 3D coordinates. Use LabuteASA as a common choice.
        # Ensure conformer exists - might need explicit generation if input was 2D/bad.
        if mol.GetNumConformers() > 0:
             metrics["SASA"] = rdMolDescriptors.CalcLabuteASA(mol)
        else:
             warnings.warn("Molecule has no conformer, cannot calculate SASA.")
             metrics["SASA"] = None # Or could try generating one: AllChem.EmbedMolecule(mol); AllChem.UFFOptimizeMolecule(mol)
        metrics["QED"] = QED.qed(mol)
        return metrics
    except Exception as e:
        warnings.warn(f"Error calculating RDKit metrics: {e}")
        return {k: None for k in ["MW", "LogP", "TPSA", "NHA", "NRB", "HBD", "HBA", "SASA", "QED"]}


def parse_vina_log(log_file):
    """Parses Vina log file to extract receptor, ligand, and best affinity.
       Uses log content first, falls back to directory structure for names.
    """
    receptor = None
    ligand = None
    affinity = None
    receptor_parsed_from_log = False
    ligand_parsed_from_log = False
    
    # --- Attempt 1: Parse Info from Log Content ---
    # Corrected regex for receptor line
    receptor_pattern = re.compile(r"Rigid receptor:\s*(.*\.pdbqt)") 
    # Ligand regex should be correct
    ligand_pattern = re.compile(r"Ligand:\s*(.*\.pdbqt)")
    # Affinity regex for mode 1
    affinity_pattern = re.compile(r"^\s+1\s+([-+]?\d*\.?\d+)") 

    try:
        with open(log_file, 'r') as f:
            for line in f:
                # Find receptor path from log
                if not receptor_parsed_from_log: # Only parse once
                    receptor_match = receptor_pattern.search(line)
                    if receptor_match:
                        # Extract base name without path/extension
                        receptor = os.path.basename(receptor_match.group(1).strip()).replace(".pdbqt", "")
                        receptor_parsed_from_log = True
                    
                # Find ligand path from log
                if not ligand_parsed_from_log: # Only parse once
                    ligand_match = ligand_pattern.search(line)
                    if ligand_match:
                         # Extract base name without path/extension
                        ligand = os.path.basename(ligand_match.group(1).strip()).replace(".pdbqt", "")
                        ligand_parsed_from_log = True

                # Find the first mode's affinity
                if affinity is None: # Only parse once
                     affinity_match = affinity_pattern.match(line)
                     if affinity_match:
                         affinity = float(affinity_match.group(1))
                         # Optimization: if we have affinity and both names from log, can stop early
                         if receptor_parsed_from_log and ligand_parsed_from_log:
                             break 
            
            # Check if affinity was found even if loop finished
            if affinity is None:
                 # Might need to re-read or search again if file is large and break was hit early?
                 # For typical Vina logs, affinity table is near the end, so this is unlikely needed.
                 pass # Keep affinity as None if not found by the end

    except FileNotFoundError:
        warnings.warn(f"Log file not found: {log_file}")
        return None, None, None
    except Exception as e:
        warnings.warn(f"Error reading log file {log_file}: {e}")
        # Proceed to fallback for names, affinity might be None

    # --- Attempt 2: Fallback for Names using Directory Structure ---
    # If parsing names from log failed, try using the directory name
    if not receptor_parsed_from_log or not ligand_parsed_from_log:
        try:
            log_dir = os.path.dirname(log_file)
            dir_basename = os.path.basename(log_dir)
            parts = dir_basename.split('_')
            if len(parts) >= 2:
                # Only overwrite if not already parsed from log
                if not receptor_parsed_from_log:
                    receptor = parts[0] 
                if not ligand_parsed_from_log:
                     # Handle ligand names potentially containing underscores
                    ligand = "_".join(parts[1:]) 
            else:
                 # Only warn if names are *still* missing after fallback attempt
                 if receptor is None or ligand is None:
                      warnings.warn(f"Could not parse receptor/ligand name from log content OR directory structure: {log_dir}")
        except Exception as e:
            # Only warn if names are *still* missing
            if receptor is None or ligand is None:
                warnings.warn(f"Error parsing directory path {log_file} for names: {e}")

    # Final checks and warnings
    if receptor is None or ligand is None:
         warnings.warn(f"Failed to determine Receptor or Ligand name for {log_file}")
         
    if affinity is None:
         warnings.warn(f"Could not parse affinity score (mode 1) from log: {log_file}")

    return receptor, ligand, affinity

def calculate_efficiency_metrics(affinity, pki, ligand_metrics):
    """Calculates LE, LLE, SILE_N, SILE_SASA."""
    metrics = { "LE": None, "LLE": None, "SILE_N": None, "SILE_SASA": None }
    
    nha = ligand_metrics.get("NHA")
    logp = ligand_metrics.get("LogP")
    sasa = ligand_metrics.get("SASA")

    # Ensure affinity is positive for LE, SILE_N calculations
    neg_affinity = -affinity if affinity is not None else None

    # LE = -Affinity / NHA
    if neg_affinity is not None and nha is not None and nha > 0:
        metrics["LE"] = neg_affinity / nha
        
    # LLE = pKi - LogP
    if pki is not None and logp is not None:
        metrics["LLE"] = pki - logp

    # SILE_N = -Affinity / (NHA ^ 0.3)
    nha_exp = 0.3 # Exponent for SILE based on Heavy Atoms
    if neg_affinity is not None and nha is not None and nha > 0:
        try:
             metrics["SILE_N"] = neg_affinity / (nha ** nha_exp)
        except ZeroDivisionError: # Should be caught by nha > 0, but defensive
             warnings.warn(f"Heavy atom count is zero for SILE_N calculation.")
        except ValueError: # e.g. negative base with non-integer exponent
             warnings.warn(f"ValueError during SILE_N calculation (NHA={nha}).")


    # SILE_SASA = pKi / SASA 
    if pki is not None and sasa is not None and sasa > 1e-6: # Avoid division by zero/tiny SASA
        metrics["SILE_SASA"] = pki / sasa
    elif sasa is not None and sasa <= 1e-6:
         warnings.warn(f"SASA is near zero ({sasa}), cannot calculate SILE_SASA.")


    return metrics


def find_ligand_file(ligand_name, search_dirs=["ligands", "ligands_pdbqt"]):
    """Tries to find the original ligand file (SDF preferred)."""
    # Prefer SDF for potentially better initial structure info
    for directory in search_dirs:
        potential_path_sdf = os.path.join(directory, f"{ligand_name}.sdf")
        if os.path.exists(potential_path_sdf):
            return potential_path_sdf, "sdf"
            
    # Fallback to PDBQT
    for directory in search_dirs:
         potential_path_pdbqt = os.path.join(directory, f"{ligand_name}.pdbqt")
         if os.path.exists(potential_path_pdbqt):
            return potential_path_pdbqt, "pdbqt"

    warnings.warn(f"Could not find source file for ligand '{ligand_name}' in {search_dirs}")
    return None, None

def load_molecule(file_path, file_type):
    """Loads molecule using RDKit based on file type."""
    mol = None
    try:
        if file_type == "sdf":
            # Use supplier to handle potential multiple molecules, but take the first
            suppl = Chem.SDMolSupplier(file_path, removeHs=False, sanitize=False)
            if suppl and len(suppl) > 0:
                 mol = suppl[0]
                 if len(suppl) > 1:
                     warnings.warn(f"SDF file {file_path} contains multiple molecules. Using only the first.")
        elif file_type == "pdbqt":
            # RDKit's PDBQT reader is less common, might need custom or rely on basic PDB reader
            # Try reading as PDB, might work for simple cases if ATOM/HETATM lines are standard
             mol = Chem.MolFromPDBFile(file_path, removeHs=False, sanitize=False)
             if mol is None: # Fallback if MolFromPDBFile fails for PDBQT specifics
                  warnings.warn(f"Could not directly read PDBQT {file_path} with RDKit MolFromPDBFile. Ligand metrics might be missing.")
                  # Potential future enhancement: Use OpenBabel API via Pybel if installed

        if mol:
             # Basic sanitization after loading
             try:
                 Chem.SanitizeMol(mol)
             except Exception as sanitize_error:
                  warnings.warn(f"RDKit sanitization failed for {file_path}: {sanitize_error}")
                  return None # Return None if sanitization fails

    except Exception as e:
        warnings.warn(f"Error loading molecule from {file_path}: {e}")
        return None # Return None if loading fails
        
    return mol


# --- Main Execution ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate docking metrics from a Vina log file.")
    parser.add_argument("vina_log_file", help="Path to the AutoDock Vina output log file.")
    parser.add_argument("--ligand_dirs", nargs='+', default=["ligands", "ligands_pdbqt"], 
                        help="Directories to search for original ligand files (SDF preferred, then PDBQT).")
    
    args = parser.parse_args()

    # 1. Parse Vina Log
    receptor_name, ligand_name, affinity = parse_vina_log(args.vina_log_file)

    # Initialize results dictionary with Nones
    results = {
        "Receptor": receptor_name, "Ligand": ligand_name, "Affinity": affinity,
        "pKi": None, "MW": None, "LogP": None, "TPSA": None, "NHA": None, "NRB": None,
        "HBD": None, "HBA": None, "SASA": None, "QED": None, "LE": None, "LLE": None,
        "SILE_N": None, "SILE_SASA": None
    }

    if receptor_name and ligand_name: # Proceed only if basic info was parsed
        # 2. Calculate pKi
        results["pKi"] = calculate_pKi(affinity)

        # 3. Find and Load Ligand File
        ligand_file_path, ligand_file_type = find_ligand_file(ligand_name, args.ligand_dirs)
        
        mol = None
        if ligand_file_path:
            mol = load_molecule(ligand_file_path, ligand_file_type)
            if mol is None:
                 warnings.warn(f"Failed to load molecule object for {ligand_name} from {ligand_file_path}")

        # 4. Calculate Ligand Metrics (intrinsic)
        if mol:
            ligand_metrics = calculate_ligand_metrics(mol)
            results.update(ligand_metrics) # Add MW, LogP etc. to results
        else:
            # Ensure ligand metrics stay None if molecule loading failed
             ligand_metrics = calculate_ligand_metrics(None) 
             results.update(ligand_metrics)

        # 5. Calculate Efficiency Metrics (dependent on affinity/pKi and ligand props)
        if affinity is not None: # Only calculate if affinity exists
            efficiency_metrics = calculate_efficiency_metrics(affinity, results["pKi"], ligand_metrics)
            results.update(efficiency_metrics)

    # 6. Format Output as CSV row (handle None values gracefully)
    header = [
        "Receptor", "Ligand", "Affinity_kcal_mol", "pKi", "MW", "LogP", "TPSA", "NHA", "NRB",
        "HBD", "HBA", "SASA_A2", "QED", "LE", "LLE", "SILE_N", "SILE_SASA"
    ]
    
    # Use this line if you want the python script to print the header (only useful if running it once)
    # print(",".join(header)) 

    # Format values, replacing None with empty strings or specific markers like 'NA'
    formatted_values = []
    for key in header:
        # Map internal keys to header keys if they differ (e.g., Affinity -> Affinity_kcal_mol)
        internal_key = key.replace('_kcal_mol','').replace('_A2','') # Basic mapping
        value = results.get(internal_key) 
        
        if value is None:
            formatted_values.append("NA") # Use NA for missing values
        elif isinstance(value, float):
             # Format floats nicely
             if key in ["Affinity_kcal_mol", "LE", "SILE_N"]:
                 formatted_values.append(f"{value:.3f}")
             elif key in ["pKi", "LogP", "TPSA", "SASA_A2", "QED", "LLE", "SILE_SASA"]:
                 formatted_values.append(f"{value:.2f}")
             elif key == "MW":
                 formatted_values.append(f"{value:.2f}")
             else:
                 formatted_values.append(f"{value:.3f}") # Default float format
        else:
            formatted_values.append(str(value)) # Convert other types (int, str) to string

    print(",".join(formatted_values))