import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np


# --- Configuration ---
CSV_FILE = "docking_summary_metrics.csv"
PLOT_DIR = "plots" 
TOP_N_LIGANDS = 10 


# Create plot directory
os.makedirs(PLOT_DIR, exist_ok=True)

# --- Load Data ---
try:
    df = pd.read_csv(CSV_FILE)
    print(f"Successfully loaded {CSV_FILE}")
    
    # 1. Clean Numeric Columns
    numeric_cols = ['Affinity_kcal_mol', 'pKi', 'MW', 'LogP', 'TPSA', 'NHA', 'NRB', 
                    'HBD', 'HBA', 'SASA_A2', 'QED', 'LE', 'LLE', 'SILE_N', 'SILE_SASA']
    
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # 2. Drop invalid rows
    df.dropna(subset=['Affinity_kcal_mol', 'Receptor', 'Ligand'], inplace=True)
    

except FileNotFoundError:
    print(f"Error: CSV file not found at {CSV_FILE}")
    exit(1)
except Exception as e:
    print(f"Error loading CSV: {e}")
    exit(1)

# --- Plotting Functions ---

# Set global style
sns.set_theme(style="whitegrid")
plt.rcParams['svg.fonttype'] = 'none'

# 1. Distribution of Affinity Scores
plt.figure(figsize=(8, 5))
sns.histplot(df['Affinity_kcal_mol'], kde=True, bins=20, color='teal')
plt.title('Distribution of Docking Affinity Scores')
plt.xlabel('Affinity (kcal/mol)')
plt.savefig(os.path.join(PLOT_DIR, '1_affinity_distribution.svg'))
plt.close()

# 2. Box Plot per Receptor
sorted_receptors = sorted(df['Receptor'].unique())
if len(sorted_receptors) > 0:
    plt.figure(figsize=(max(2, len(sorted_receptors)*0.75), 6))
    sns.boxplot(x='Receptor', y='Affinity_kcal_mol', data=df, order=sorted_receptors, hue='Receptor', palette="icefire", legend=False)
    sns.stripplot(x='Receptor', y='Affinity_kcal_mol', data=df, order=sorted_receptors, 
                  color='black', alpha=0.3, size=3, jitter=True) # Add points for detail
    plt.title('Affinity Scores per Receptor')
    plt.xticks(rotation=0, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, '2_affinity_boxplot.svg'))
    plt.close()

# 3. Pareto Plot (Affinity vs MW) - The "Drug Discovery" Plot
# We want to highlight the bottom-left corner (High Affinity [more negative], Low MW)
plt.figure(figsize=(10, 7))
sns.scatterplot(data=df, x='MW', y='Affinity_kcal_mol', hue='Receptor', style='Receptor', s=100, alpha=0.8)
plt.title('Pareto Efficiency: Affinity vs. Molecular Weight')
plt.xlabel('Molecular Weight (Da)')
plt.ylabel('Affinity (kcal/mol)')
plt.axhline(y=-7.0, color='r', linestyle='--', alpha=0.5, label='Hit Cutoff (-7 kcal/mol)')
plt.axvline(x=500, color='r', linestyle='--', alpha=0.5, label='Lipinski Cutoff (500 Da)')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, '3_pareto_affinity_vs_mw.svg'))
plt.close()

# 4. Top N Ligands per Receptor (Bar Chart)
for receptor in sorted_receptors:
    # Get top binders (most negative affinity)
    top_df = df[df['Receptor'] == receptor].sort_values('Affinity_kcal_mol', ascending=True).head(TOP_N_LIGANDS)
    
    if top_df.empty: continue
        
    plt.figure(figsize=(8, max(4, len(top_df)*0.4)))
    sns.barplot(x='Affinity_kcal_mol', y='Ligand', data=top_df, hue='Ligand', palette='icefire', legend=False)
    plt.title(f'Top {len(top_df)} Ligands: {receptor}')
    plt.xlabel('Affinity (kcal/mol)')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, f'4_top_ligands_{receptor}.svg'))
    plt.close()

# 5. Correlation Heatmap
numeric_df = df.select_dtypes(include=[np.number])
if not numeric_df.empty:
    plt.figure(figsize=(10, 8))
    corr = numeric_df.corr()
    sns.heatmap(corr, annot=True, cmap='icefire', fmt=".2f", vmin=-1, vmax=1)
    plt.title('Correlation of Metrics')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, '5_correlation_heatmap.svg'))
    plt.close()

# --- NEW: Generate Text Summary Report ---
summary_file = os.path.join(PLOT_DIR, 'top_candidates_summary.txt')
with open(summary_file, 'w') as f:
    f.write("=== DOCKING CANDIDATE SUMMARY ===\n")
    f.write(f"Total pairs analyzed: {len(df)}\n\n")
    
    for receptor in sorted_receptors:
        rec_df = df[df['Receptor'] == receptor].sort_values('Affinity_kcal_mol', ascending=True)
        best_ligand = rec_df.iloc[0]
        
        f.write(f"--- Receptor: {receptor} ---\n")
        f.write(f"Top Candidate: {best_ligand['Ligand']}\n")
        f.write(f"  Affinity: {best_ligand['Affinity_kcal_mol']} kcal/mol\n")
        if 'pKi' in best_ligand:
            f.write(f"  pKi: {best_ligand['pKi']:.2f}\n")
        if 'LE' in best_ligand:
            f.write(f"  Ligand Efficiency: {best_ligand['LE']:.3f}\n")
        
        f.write("\nTop 5 List:\n")
        for i, row in rec_df.head(5).iterrows():
            f.write(f"  {i+1}. {row['Ligand']} ({row['Affinity_kcal_mol']} kcal/mol)\n")
        f.write("\n")

print(f"\nAnalysis Complete!")
print(f"Plots saved in: {PLOT_DIR}")
print(f"Summary report: {summary_file}")
