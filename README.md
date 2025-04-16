# Automated Molecular Docking Pipeline with Vina and RDKit

## Description

This project provides a streamlined pipeline for performing automated molecular docking using AutoDock Vina. It takes receptor structures (PDB) and ligand files (SDF, MOL2, MOL, PDB, SMILES) as input, prepares them, runs docking simulations for all receptor-ligand pairs, extracts binding affinities, calculates various physicochemical properties and ligand efficiency metrics using RDKit, and generates summary plots for analysis.

The pipeline is designed to be modular and run via a master script.

## Prerequisites

*   **Conda:** Anaconda or Miniconda package manager installed. ([https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html))
*   **Core Libraries:** `openbabel`, `python`, `numpy`, `rdkit`, `pandas`, `matplotlib`, `seaborn`. These will be installed via Conda.
*   **AutoDock Vina:** Can be installed via Conda (recommended) or separately. The `vina` command must be accessible from the command line when the Conda environment is activated. ([http://vina.scripps.edu/](http://vina.scripps.edu/) or [https://github.com/ccsb-scripps/AutoDock-Vina](https://github.com/ccsb-scripps/AutoDock-Vina))

## Setup

1.  **Clone or Download:** Get the project files onto your system. If using Git:
    ```bash
    git clone <repository_url>
    cd vina_pipeline
    ```
    If you downloaded a ZIP, extract it.

2.  **Create Conda Environment:** Create a dedicated environment and install necessary packages. This command includes Vina:
    ```bash
    conda create --name docking_env -c conda-forge openbabel rdkit vina python=3.12.9 numpy pandas matplotlib seaborn -y
    ```
    

3.  **Activate Environment:** Before running the pipeline, always activate the environment:
    ```bash
    conda activate docking_env
    ```

4.  **Verify Installations (Optional but recommended):**
    ```bash
    conda list openbabel # Check if listed
    conda list rdkit     # Check if listed
    vina --version       # Check Vina version
    python --version     # Check Python version
    ```

## Directory Structure

The project expects the following structure. You need to create the `receptors` and `ligands` directories and populate them with your input files. The other directories will be created automatically by the pipeline.

```
vina_pipeline/
├── ligands/ # INPUT: Place your original ligand files here (SDF, MOL2, etc.)
├── receptors/ # INPUT: Place your original receptor PDB files here
├── src/ # Contains all executable pipeline scripts
│ ├── 1_receptor_prep.sh
│ ├── 2_ligand_prep.sh
│ ├── 3_config_calc.sh
│ ├── 4_run_vina.sh
│ ├── 5_analyze_results.sh
│ ├── 6_plot_results.py
│ ├── calc_box_center.py
│ └── calculate_docking_metrics.py
└── master.sh # Main script to run the entire pipeline
```
--- Output directories will be created here ---
```
├── receptors_pdbqt/ # OUTPUT: Prepared receptors in PDBQT format
├── ligands_pdbqt/ # OUTPUT: Prepared ligands in PDBQT format
├── configs/ # OUTPUT: Vina configuration files (.txt)
├── docking_results/ # OUTPUT: Raw docking results (poses, logs) per pair
├── plots/ # OUTPUT: Generated analysis plots (.png)
└── docking_summary_metrics.csv # OUTPUT: Final summary table
```
## Input Data Preparation

1.  **Receptors:** Place your receptor protein structures in PDB format (`.pdb`) inside the `receptors/` directory.
2.  **Ligands:** Place your ligand molecule files inside the `ligands/` directory.
    *   Supported formats: SDF (`.sdf`), MOL2 (`.mol2`), MOL (`.mol`), PDB (`.pdb`), SMILES (`.smi`, `.smiles`).
    *   The script uses the filename (without extension) as the unique identifier for each ligand in the final results table. Ensure filenames are meaningful and unique.
    *   If using SMILES, ensure they are valid and represent the intended molecules. 3D coordinates will be generated automatically.

## Running the Pipeline

1.  **Activate Conda Environment:**
    ```bash
    conda activate docking_env
    ```

2.  **Navigate to Project Root:** Make sure you are in the main project directory (`vina_pipeline/`).

3.  **Make Scripts Executable (if needed):** This usually only needs to be done once.
    ```bash
    chmod +x master.sh src/*.sh
    ```

4.  **Run the Master Script:**
    ```bash
    ./master.sh
    ```

The `master.sh` script will execute each step of the pipeline sequentially. Progress messages will be printed to the terminal. Check for any warnings or errors.

## Output Description

Upon successful completion, the following outputs will be generated in the main project directory:

*   **`receptors_pdbqt/`**: Contains receptor files prepared for Vina (PDBQT format, with hydrogens and charges).
*   **`ligands_pdbqt/`**: Contains ligand files prepared for Vina (PDBQT format, with hydrogens and charges).
*   **`configs/`**: Contains `.txt` configuration files for Vina, one per receptor, defining the docking search box (calculated automatically to encompass the entire receptor).
*   **`docking_results/`**: Contains subdirectories for each `receptor_ligand` pair. Each subdirectory includes:
    *   `docked_poses.pdbqt`: Predicted binding poses from Vina.
    *   `vina_output.log`: The raw text output from the Vina run for that pair, including the affinity score table.
*   **`docking_summary_metrics.csv`**: A comma-separated values file containing the aggregated results and calculated metrics for all successful docking runs. Columns include:
    *   `Receptor`, `Ligand`
    *   `Affinity_kcal_mol`: Best Vina score
    *   `pKi`: Calculated pKi from Affinity
    *   `MW`: Molecular Weight
    *   `LogP`: Calculated LogP
    *   `TPSA`: Topological Polar Surface Area
    *   `NHA`: Number of Heavy Atoms
    *   `NRB`: Number of Rotatable Bonds
    *   `HBD`: Number of Hydrogen Bond Donors
    *   `HBA`: Number of Hydrogen Bond Acceptors
    *   `SASA_A2`: Solvent Accessible Surface Area (Å²)
    *   `QED`: Quantitative Estimate of Drug-likeness
    *   `LE`: Ligand Efficiency (Affinity / NHA)
    *   `LLE`: Lipophilic Ligand Efficiency (pKi - LogP)
    *   `SILE_N`: Size-Independent LE (Affinity / NHA^0.3)
    *   `SILE_SASA`: Surface-based SILE (pKi / SASA)
*   **`plots/`**: Contains various plots in `.png` format visualizing the results, such as:
    *   Affinity distributions (overall and per receptor)
    *   Top ligand rankings per receptor
    *   Scatter plots (Affinity vs MW, pKi vs LogP, LE vs Affinity, SILE_SASA vs pKi)
    *   Correlation heatmap of metrics
    *   Pair plot of selected metrics

## Script Descriptions

*   **`master.sh`**: Orchestrates the execution of all pipeline steps in order.
*   **`src/1_receptor_prep.sh`**: Prepares receptor PDB files using Open Babel.
*   **`src/2_ligand_prep.sh`**: Prepares ligand files (detecting format) using Open Babel.
*   **`src/3_config_calc.sh`**: Calls `calc_box_center.py` for each receptor to generate Vina config files.
*   **`src/calc_box_center.py`**: Calculates Vina docking box parameters covering the entire receptor.
*   **`src/4_run_vina.sh`**: Runs AutoDock Vina for all receptor-ligand pairs.
*   **`src/5_analyze_results.sh`**: Calls `calculate_docking_metrics.py` for each docking log to parse results and calculate metrics, compiling them into the summary CSV.
*   **`src/calculate_docking_metrics.py`**: Parses individual Vina logs, uses RDKit to calculate properties and efficiency metrics.
*   **`src/6_plot_results.py`**: Reads the summary CSV and generates various analysis plots using Matplotlib and Seaborn.

## Customization

*   **Vina Parameters:** Modify `EXHAUSTIVENESS` and `NUM_MODES` in `src/4_run_vina.sh`.
*   **Docking Box:** Change the `BUFFER_SIZE` in `src/3_config_calc.sh`. For **focused docking**, you would need to significantly modify `src/calc_box_center.py` or generate config files manually based on known binding sites.
*   **Ligand Formats:** Add support for more input formats in the `case` statement within `src/2_ligand_prep.sh`.
*   **Multi-Molecule Ligand Files:** If your input files (e.g., SDF) contain multiple molecules per file and you want to dock *all* of them, uncomment the `-m` flag in the `obabel` command within `src/2_ligand_prep.sh`. Note this will change output naming in `ligands_pdbqt/`.
*   **Plotting:** Adjust `TOP_N_LIGANDS`, plot styles, selected metrics for pair plots, etc., in `src/6_plot_results.py`.
*   **Error Handling:** The pipeline uses `set -e` in most bash scripts to stop on error. `5_analyze_results.sh` currently has this commented out to allow processing remaining logs even if one fails; you can change this behaviour.

## Citation

If you use this pipeline or its components in your research, please cite the relevant software:

*   **AutoDock Vina**
*   **Open Babel**
*   **RDKit**
*   **NumPy**
*   **Pandas**
*   **Matplotlib**
*   **Seaborn**
