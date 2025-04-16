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
    cd MyDockingProject # Or your chosen project directory name
    ```
    If you downloaded a ZIP, extract it.

2.  **Create Conda Environment:** Create a dedicated environment and install necessary packages. This command includes Vina:
    ```bash
    conda create --name docking_env -c conda-forge openbabel rdkit vina python=3.9 numpy pandas matplotlib seaborn -y
    ```
    *(You can adjust the `python=3.9` version if needed, but ensure compatibility with the libraries. Using a specific version like 3.9, 3.10, or 3.11 is often more stable than the latest).*

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


MyDockingProject/
├── ligands/ # INPUT: Place your original ligand files here (SDF, MOL2, etc.)
├── receptors/ # INPUT: Place your original receptor PDB files here
├── src/ # Contains all executable pipeline scripts
│ ├── 01_prepare_receptors.sh
│ ├── 02_prepare_ligands.sh
│ ├── 03_generate_configs.sh
│ ├── 04_run_docking.sh
│ ├── 05_extract_analyze_results.sh
│ ├── 06_plot_results.py
│ ├── calc_box_center.py
│ └── calculate_docking_metrics.py
├── master.sh # Main script to run the entire pipeline

--- Output directories will be created here ---

├── receptors_pdbqt/ # OUTPUT: Prepared receptors in PDBQT format
├── ligands_pdbqt/ # OUTPUT: Prepared ligands in PDBQT format
├── configs/ # OUTPUT: Vina configuration files (.txt)
├── docking_results/ # OUTPUT: Raw docking results (poses, logs) per pair
├── plots/ # OUTPUT: Generated analysis plots (.png)
└── docking_summary_metrics.csv # OUTPUT: Final summary table

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

2.  **Navigate to Project Root:** Make sure you are in the main project directory (`MyDockingProject/`).

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
*   **`src/01_prepare_receptors.sh`**: Prepares receptor PDB files using Open Babel.
*   **`src/02_prepare_ligands.sh`**: Prepares ligand files (detecting format) using Open Babel.
*   **`src/03_generate_configs.sh`**: Calls `calc_box_center.py` for each receptor to generate Vina config files.
*   **`src/calc_box_center.py`**: Calculates Vina docking box parameters covering the entire receptor.
*   **`src/04_run_docking.sh`**: Runs AutoDock Vina for all receptor-ligand pairs.
*   **`src/05_extract_analyze_results.sh`**: Calls `calculate_docking_metrics.py` for each docking log to parse results and calculate metrics, compiling them into the summary CSV.
*   **`src/calculate_docking_metrics.py`**: Parses individual Vina logs, uses RDKit to calculate properties and efficiency metrics.
*   **`src/06_plot_results.py`**: Reads the summary CSV and generates various analysis plots using Matplotlib and Seaborn.

## Customization

*   **Vina Parameters:** Modify `EXHAUSTIVENESS` and `NUM_MODES` in `src/04_run_docking.sh`.
*   **Docking Box:** Change the `BUFFER_SIZE` in `src/03_generate_configs.sh`. For **focused docking**, you would need to significantly modify `src/calc_box_center.py` or generate config files manually based on known binding sites.
*   **Ligand Formats:** Add support for more input formats in the `case` statement within `src/02_prepare_ligands.sh`.
*   **Multi-Molecule Ligand Files:** If your input files (e.g., SDF) contain multiple molecules per file and you want to dock *all* of them, uncomment the `-m` flag in the `obabel` command within `src/02_prepare_ligands.sh`. Note this will change output naming in `ligands_pdbqt/`.
*   **Plotting:** Adjust `TOP_N_LIGANDS`, plot styles, selected metrics for pair plots, etc., in `src/06_plot_results.py`.
*   **Error Handling:** The pipeline uses `set -e` in most bash scripts to stop on error. `05_extract_analyze_results.sh` currently has this commented out to allow processing remaining logs even if one fails; you can change this behaviour.

## Citation

If you use this pipeline or its components in your research, please cite the relevant software:

*   **AutoDock Vina:**
    *   J. Eberhardt, D. Santos-Martins, A. F. Tillack, and S. Forli, AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings, J. Chem. Inf. Model. (2021) DOI: 10.1021/acs.jcim.1c00203
    *   O. Trott, A. J. Olson, AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization and multithreading, J. Comp. Chem. (2010) DOI: 10.1002/jcc.21334
*   **Open Babel:**
    *   N M O'Boyle et al., "Open Babel: An open chemical toolbox." J Cheminform (2011) 3:33. DOI: 10.1186/1758-2946-3-33
*   **RDKit:**
    *   RDKit: Open-source cheminformatics; http://www.rdkit.org
*   **NumPy:**
    *   Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2
*   **Pandas:**
    *   Wes McKinney. Data Structures for Statistical Computing in Python, Proceedings of the 9th Python in Science Conference, 51-56 (2010)
    *   pandas development team. pandas-dev/pandas: Pandas (DOI: 10.5281/zenodo.3509134)
*   **Matplotlib:**
    *   J. D. Hunter, "Matplotlib: A 2D Graphics Environment", Computing in Science & Engineering, vol. 9, no. 3, pp. 90-95, (2007) DOI: 10.1109/MCSE.2007.55
*   **Seaborn:**
    *   Michael L. Waskom, "Seaborn: statistical data visualization", Journal of Open Source Software, 6(60), 3021, (2021) DOI: 10.21105/joss.03021

## License

(Optional: Consider adding a license file, e.g., MIT, Apache 2.0, if you plan to share this project.)
IGNORE_WHEN_COPYING_START
content_copy
download
Use code with caution.
IGNORE_WHEN_COPYING_END
