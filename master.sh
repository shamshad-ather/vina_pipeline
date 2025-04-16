#!/bin/bash

# Master script to run the entire docking pipeline
# This script should be executed from the project's root directory
# (e.g., MyDockingProject/)

# Exit immediately if a command exits with a non-zero status.
set -e 

# Define the source directory where scripts reside
SRC_DIR="src"

echo "========================================="
echo " STEP 1: Preparing Receptors             "
echo "========================================="
bash "${SRC_DIR}/1_receptor_prep.sh"
echo " Receptor preparation finished."
echo ""

echo "========================================="
echo " STEP 2: Preparing Ligands               "
echo "========================================="
bash "${SRC_DIR}/2_ligand_prep.sh"
echo " Ligand preparation finished."
echo ""

echo "========================================="
echo " STEP 3: Generating Vina Configs       "
echo "========================================="
bash "${SRC_DIR}/3_config_calc.sh"
echo " Config generation finished."
echo ""

echo "========================================="
echo " STEP 4: Running AutoDock Vina Docking   "
echo "========================================="
bash "${SRC_DIR}/4_run_vina.sh"
echo " Docking finished."
echo ""

echo "========================================="
echo " STEP 5: Extracting & Analyzing Results"
echo "========================================="
bash "${SRC_DIR}/5_analyze_results.sh"
echo " Results analysis finished."
echo ""

echo "========================================="
echo " STEP 6: Generating Plots              "
echo "========================================="
# Plotting script is Python, execute with python interpreter
python "${SRC_DIR}/6_plot_results.py"
echo " Plot generation finished."
echo ""

echo "========================================="
echo " PIPELINE COMPLETED SUCCESSFULLY        "
echo "========================================="

exit 0
