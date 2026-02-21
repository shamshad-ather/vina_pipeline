# Automated Molecular Docking Pipeline with Vina and RDKit

## Description

This project provides a streamlined pipeline for performing automated molecular docking using AutoDock Vina. It takes receptor structures (PDB) and ligand files as input, prepares them, runs docking simulations for all receptor-ligand pairs, extracts binding affinities, calculates various physicochemical properties and ligand efficiency metrics using RDKit, and generates summary plots for analysis.

The pipeline is designed to be modular and run via a master script, and supports **two input workflows** for maximum flexibility.

### Quick Start: Choose Your Workflow

| **Workflow**                     | **Input Type**                 | **When to Use**                                                                              |
| -------------------------------- | ------------------------------ | -------------------------------------------------------------------------------------------- |
| **Option 1: SMILES**             | Excel file with SMILES strings | You have chemical structures as SMILES notation or want to define compounds programmatically |
| **Option 2: Pre-prepared Files** | SDF, MOL2, MOL, or PDB files   | You already have 3D structure files ready for docking                                        |

**See [Input Data Preparation](#input-data-preparation) below for detailed setup for each workflow.**

## Prerequisites

- **Conda:** Anaconda or Miniconda package manager installed. ([https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html))
- **Core Libraries:** `openbabel`, `python`, `numpy`, `rdkit`, `pandas`, `matplotlib`, `seaborn`. These will be installed via Conda.
- **AutoDock Vina:** Can be installed via Conda (recommended) or separately. The `vina` command must be accessible from the command line when the Conda environment is activated. ([http://vina.scripps.edu/](http://vina.scripps.edu/) or [https://github.com/ccsb-scripps/AutoDock-Vina](https://github.com/ccsb-scripps/AutoDock-Vina))

## Setup

1.  **Clone or Download:** Get the project files onto your system. If using Git:

    ```bash
    git clone https://github.com/shamshad-ather/vina_pipeline.git
    cd vina_pipeline
    ```

    If you downloaded a ZIP, extract it.

2.  **Create Conda Environment:** Create a dedicated environment and install necessary packages. This command includes Vina:

    ```bash
    conda create --name docking -c conda-forge openbabel rdkit vina python=3.12.9 numpy pandas matplotlib seaborn -y
    ```

3.  **Activate Environment:** Before running the pipeline, always activate the environment:

    ```bash
    conda activate docking
    ```

4.  **Verify Installations (Optional but recommended):**
    ```bash
    conda list openbabel # Check if listed
    conda list rdkit     # Check if listed
    vina --version       # Check Vina version
    python --version     # Check Python version
    ```
5.  **Install openpyxl**
    ```bash
    pip install openpyxl
    ```

## Directory Structure

The project expects the following structure. You may need to create some directories and populate them based on your chosen input workflow. Other directories will be created automatically by the pipeline.

```
vina_pipeline/
├── input_smiles.xlsx # OPTIONAL INPUT: Excel file for SMILES workflow (if using Workflow Option 1)
├── ligands/ # INPUT: Ligand files directory
│   │       # Choice A (Workflow 1): Will be auto-populated by 1_5_ligand_gen.py from SMILES
│   │       # Choice B (Workflow 2): Place your pre-prepared ligand files here (SDF, MOL2, MOL, PDB)
│   ├── Ellagic acid.sdf (example)
│   ├── Kanzonol W.sdf (example)
│   └── ...
├── receptors/ # INPUT: Place your PDB receptor files here
│   ├── CDK6.pdb (example)
│   ├── Protein_A.pdb (example)
│   └── ...
└── src/ # Contains all executable pipeline scripts
    ├── 1_5_ligand_gen.py
    ├── 1_receptor_prep.sh
    ├── 2_ligand_prep.sh
    ├── 3_config_calc.sh
    ├── 4_run_vina.sh
    ├── 5_analyze_results.sh
    ├── 6_plot_results.py
    ├── calc_box_center.py
    ├── calculate_docking_metrics.py
    └── master.sh # Main script to run the entire pipeline
```

--- Output directories (created automatically) ---

```
├── receptors_pdbqt/ # OUTPUT: Prepared receptors in PDBQT format
├── ligands_pdbqt/ # OUTPUT: Prepared ligands in PDBQT format
├── configs/ # OUTPUT: Vina configuration files (.txt)
├── docking_results/ # OUTPUT: Raw docking results (poses, logs) per pair
├── plots/ # OUTPUT: Generated analysis plots (.png)
└── docking_summary_metrics.csv # OUTPUT: Final summary table
```

## Input Data Preparation

The pipeline supports **two input workflows** for ligands:

### Workflow Option 1: SMILES Input (via Excel File)

If you want to provide your ligands as SMILES strings:

1. **Create Excel File**: Prepare an Excel file named `input_smiles.xlsx` in the project root directory with the following structure:
   - **Column 1: `CID`** – Compound identifier (used for filenames and results table). Example: `Ellagic_acid`, `Compound_001`
   - **Column 2: `SMILES`** – Valid SMILES string representation of the molecule. Example: `CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O`

   Example Excel format:

   ```
   CID              | SMILES
   Ellagic acid     | C1=CC(=C(C(=C1)O)O)C(=O)O
   Kanzonol W       | CC(C)=CCC(C)C(=CC(=CC(=CC(=CC(=CC(=CC(=C(...))O)O)O)O)O)O)O
   ```

2. **Run Pipeline**: Execute `master.sh`. Step 1.5 (Ligand Structure Generation) will:
   - Read `input_smiles.xlsx`
   - Generate 3D structures using RDKit ETKDGv3 with 50 conformers
   - Optimize with MMFF94s force field
   - Save best conformer as `.sdf` files to the `ligands/` directory
   - Convert to PDBQT format for docking

**Remove the Excel file after generation to avoid processing the outdated file in future runs.**

### Workflow Option 2: Pre-prepared Ligand Files

If you already have ligand files in standard formats:

1.  **Place Ligand Files**: Put your ligand molecule files directly into the `ligands/` directory.
    - Supported formats: SDF (`.sdf`), MOL2 (`.mol2`), MOL (`.mol`), PDB (`.pdb`)
    - The filename (without extension) is used as the unique identifier in the results table. Ensure filenames are meaningful and unique.
    - Ensure files are already in their correct molecular form (correct protonation, stereochemistry, etc.)

2.  **Run Pipeline**: Execute `master.sh`. Steps 2+ will process the files directly without the ligand generation step.

**Note:** Step 1.5 (Ligand Structure Generation) will only read `input_smiles.xlsx` if it exists. You can safely skip this step by ensuring no `input_smiles.xlsx` file is present.

### Common Input Setup

1.  **Receptors:** Place your receptor protein structures in PDB format (`.pdb`) inside the `receptors/` directory.
    - The filename (without extension) identifies the receptor in results.
    - Ensure PDB files contain valid protein structures.

## Running the Pipeline

### Workflow Selection & Setup

**Choose ONE of the following workflows:**

#### Workflow Option 1: SMILES Input

1. Create `input_smiles.xlsx` with columns `CID` and `SMILES` (see [Input Data Preparation](#input-data-preparation) for example)
2. Place receptor PDB files in `receptors/` directory
3. Ensure `receptors/` and `ligands/` directories exist (create if needed)
4. Do NOT place any files in `ligands/` for this workflow (they will be generated)

#### Workflow Option 2: Pre-prepared Files

1. Place ligand files (SDF, MOL2, MOL, PDB) in `ligands/` directory
2. Place receptor PDB files in `receptors/` directory
3. Ensure NO `input_smiles.xlsx` file exists
4. Ensure both directories exist and contain your files

### Execute Pipeline

1.  **Activate Conda Environment:**

    ```bash
    conda activate docking
    ```

2.  **Navigate to Project Root:** Make sure you are in the main project directory.

3.  **Make Scripts Executable (if needed):** This usually only needs to be done once.

    ```bash
    chmod +x master.sh src/*.sh
    ```

4.  **Run the Master Script:**
    ```bash
    ./master.sh
    ```

The `master.sh` script will execute each step of the pipeline sequentially. Progress messages will be printed to the terminal. Check for any warnings or errors.

**Important:** After completing a run with Workflow Option 1 (SMILES), delete or move the `input_smiles.xlsx` file before the next run to avoid unintended ligand regeneration.

## Output Description

Upon successful completion, the following outputs will be generated in the main project directory:

- **`receptors_pdbqt/`**: Contains receptor files prepared for Vina (PDBQT format, with hydrogens and charges).
- **`ligands_pdbqt/`**: Contains ligand files prepared for Vina (PDBQT format, with hydrogens and charges).
- **`configs/`**: Contains `.txt` configuration files for Vina, one per receptor, defining the docking search box (calculated automatically to encompass the entire receptor).
- **`docking_results/`**: Contains subdirectories for each `receptor_ligand` pair. Each subdirectory includes:
  - `docked_poses.pdbqt`: Predicted binding poses from Vina.
  - `vina_output.log`: The raw text output from the Vina run for that pair, including the affinity score table.
- **`docking_summary_metrics.csv`**: A comma-separated values file containing the aggregated results and calculated metrics for all successful docking runs. Columns include:
  - `Receptor`, `Ligand`
  - `Affinity_kcal_mol`: Best Vina score
  - `pKi`: Calculated pKi from Affinity
  - `MW`: Molecular Weight
  - `LogP`: Calculated LogP
  - `TPSA`: Topological Polar Surface Area
  - `NHA`: Number of Heavy Atoms
  - `NRB`: Number of Rotatable Bonds
  - `HBD`: Number of Hydrogen Bond Donors
  - `HBA`: Number of Hydrogen Bond Acceptors
  - `SASA_A2`: Solvent Accessible Surface Area (Å²)
  - `QED`: Quantitative Estimate of Drug-likeness
  - `LE`: Ligand Efficiency (Affinity / NHA)
  - `LLE`: Lipophilic Ligand Efficiency (pKi - LogP)
  - `SILE_N`: Size-Independent LE (Affinity / NHA^0.3)
  - `SILE_SASA`: Surface-based SILE (pKi / SASA)
- **`plots/`**: Contains various plots in `.svg` format visualizing the results, such as:
  - Affinity distributions (overall and per receptor)
  - Top ligand rankings per receptor
  - Scatter plots (Affinity vs MW, pKi vs LogP, LE vs Affinity, SILE_SASA vs pKi)
  - Correlation heatmap of metrics
  - Pair plot of selected metrics

## Script Descriptions

- **`master.sh`**: Orchestrates the execution of all pipeline steps in order:
  - Step 1: Receptor preparation
  - Step 1.5: Ligand structure generation (only runs if `input_smiles.xlsx` exists; specific to Workflow Option 1)
  - Step 2: Ligand preparation
  - Steps 3-6: Config generation, docking, analysis, and plotting

- **`src/1_5_ligand_gen.py`**: **[Workflow Option 1 only]** Processes SMILES from Excel file:
  - Reads `input_smiles.xlsx` with columns `CID` (compound ID) and `SMILES` (SMILES string)
  - Generates 3D conformers using RDKit ETKDGv3 (state-of-the-art 3D embedding)
  - Optimizes with MMFF94s or UFF force field (50 conformers → selects lowest energy)
  - Outputs `.sdf` files to `ligands/` directory for subsequent processing
  - Uses parallel processing (all CPU cores) for speed

- **`src/1_receptor_prep.sh`**: Prepares receptor PDB files using Open Babel.
  - Converts to PDBQT format with hydrogens and Gasteiger charges

- **`src/2_ligand_prep.sh`**: Prepares ligand files (applies to both workflow options):
  - **Workflow Option 1**: Processes `.sdf` files generated by `1_5_ligand_gen.py`
  - **Workflow Option 2**: Processes user-provided ligand files (SDF, MOL2, MOL, PDB)
  - Converts to PDBQT format with hydrogens and Gasteiger charges

- **`src/3_config_calc.sh`**: Calls `calc_box_center.py` for each receptor to generate Vina config files.

- **`src/calc_box_center.py`**: Calculates Vina docking box parameters covering the entire receptor.

- **`src/4_run_vina.sh`**: Runs AutoDock Vina for all receptor-ligand pairs in parallel.

- **`src/5_analyze_results.sh`**: Calls `calculate_docking_metrics.py` for each docking log to parse results and calculate metrics, compiling them into the summary CSV.

- **`src/calculate_docking_metrics.py`**: Parses individual Vina logs, uses RDKit to calculate properties and efficiency metrics.

- **`src/6_plot_results.py`**: Reads the summary CSV and generates various analysis plots using Matplotlib and Seaborn.

## Customization

### Workflow Selection

- **To use Workflow Option 1 (SMILES Input):** Create an `input_smiles.xlsx` file in the project root with columns `CID` and `SMILES`. Run `master.sh`. Ensure the file is deleted or renamed before running the pipeline again.
- **To use Workflow Option 2 (Pre-prepared Files):** Place ligand files directly in the `ligands/` directory. Ensure NO `input_smiles.xlsx` file exists. Run `master.sh`.

### Advanced Customization

- **SMILES Generation Parameters** (in `src/1_5_ligand_gen.py`):
  - `NUM_CONFS`: Number of conformers to generate (default: 50). Higher values = longer generation time but potentially better conformer sampling.
  - `SHEET_NAME`: Excel sheet index or name (default: 0 for first sheet).
  - `CID_COL` and `SMILES_COL`: Column names in Excel (default: 'CID' and 'SMILES').

- **Vina Parameters:** Modify `EXHAUSTIVENESS` and `NUM_MODES` in `src/4_run_vina.sh`.

- **Docking Box:** Change the `BUFFER_SIZE` in `src/3_config_calc.sh`. For **focused docking**, you would need to significantly modify `src/calc_box_center.py` or generate config files manually based on known binding sites.

- **Ligand Preparation:** Add support for more input formats in the `case` statement within `src/2_ligand_prep.sh`.

- **Multi-Molecule Ligand Files:** If your input files (e.g., SDF) contain multiple molecules per file and you want to dock _all_ of them, uncomment the `-m` flag in the `obabel` command within `src/2_ligand_prep.sh`. Note this will change output naming in `ligands_pdbqt/`.

- **Plotting:** Adjust `TOP_N_LIGANDS`, plot styles, selected metrics for pair plots, etc., in `src/6_plot_results.py`.

- **Error Handling:** The pipeline uses `set -e` in most bash scripts to stop on error. `5_analyze_results.sh` currently has this commented out to allow processing remaining logs even if one fails; you can change this behaviour.

## Citation

If you use this pipeline or its components in your research, please cite the relevant software:

- **AutoDock Vina**
- **Open Babel**
- **RDKit**
- **NumPy**
- **Pandas**
- **Matplotlib**
- **Seaborn**
