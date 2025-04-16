import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np # For handling potential inf/-inf during correlation

# --- Configuration ---
CSV_FILE = "docking_summary_metrics.csv"
PLOT_DIR = "plots" 
TOP_N_LIGANDS = 10 # How many top ligands to show in ranking plots

# Create plot directory if it doesn't exist
os.makedirs(PLOT_DIR, exist_ok=True)

# --- Load Data ---
try:
    df = pd.read_csv(CSV_FILE)
    print(f"Successfully loaded {CSV_FILE}")
    print(f"Data shape: {df.shape}")
    print("Columns:", df.columns.tolist())
    
    # Convert relevant columns to numeric, coercing errors to NaN
    numeric_cols = ['Affinity_kcal_mol', 'pKi', 'MW', 'LogP', 'TPSA', 'NHA', 'NRB', 
                    'HBD', 'HBA', 'SASA_A2', 'QED', 'LE', 'LLE', 'SILE_N', 'SILE_SASA']
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        else:
             print(f"Warning: Expected numeric column '{col}' not found in CSV.")

    # Drop rows where key metrics like Affinity are missing (optional, but good for plotting)
    df.dropna(subset=['Affinity_kcal_mol', 'Receptor', 'Ligand'], inplace=True)
    print(f"Data shape after dropping rows with missing Affinity/Receptor/Ligand: {df.shape}")
    
except FileNotFoundError:
    print(f"Error: CSV file not found at {CSV_FILE}")
    exit(1)
except Exception as e:
    print(f"Error loading or processing CSV: {e}")
    exit(1)

# --- Plotting Functions ---

# 1. Distribution of Affinity Scores (Overall)
plt.figure(figsize=(8, 5))
sns.histplot(df['Affinity_kcal_mol'], kde=True, bins=20)
plt.title('Distribution of Docking Affinity Scores (All Pairs)')
plt.xlabel('Affinity (kcal/mol)')
plt.ylabel('Frequency')
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, 'affinity_distribution_overall.png'))
plt.close()
print("Generated: affinity_distribution_overall.png")

# 2. Affinity Scores per Receptor (Box Plot)
# Sort receptors alphabetically for consistent plot order
sorted_receptors = sorted(df['Receptor'].unique())
if len(sorted_receptors) > 1: # Only useful if more than one receptor
    plt.figure(figsize=(min(15, len(sorted_receptors)*2), 6)) # Adjust width based on number of receptors
    sns.boxplot(x='Receptor', y='Affinity_kcal_mol', data=df, order=sorted_receptors)
    plt.title('Affinity Scores Distribution per Receptor')
    plt.xlabel('Receptor')
    plt.ylabel('Affinity (kcal/mol)')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, 'affinity_boxplot_per_receptor.png'))
    plt.close()
    print("Generated: affinity_boxplot_per_receptor.png")
else:
    print("Skipping boxplot per receptor (only one receptor found).")

# 3. Top N Ligands per Receptor (Ranked Bar Chart)
for receptor in sorted_receptors:
    receptor_df = df[df['Receptor'] == receptor].sort_values('Affinity_kcal_mol', ascending=True)
    top_ligands_df = receptor_df.head(TOP_N_LIGANDS)
    
    if top_ligands_df.empty:
        print(f"No data to plot for receptor {receptor}")
        continue
        
    plt.figure(figsize=(10, max(5, len(top_ligands_df)*0.5))) # Adjust height based on N
    sns.barplot(x='Affinity_kcal_mol', y='Ligand', data=top_ligands_df, palette='viridis')
    plt.title(f'Top {len(top_ligands_df)} Ligands for Receptor {receptor} by Affinity')
    plt.xlabel('Affinity (kcal/mol)')
    plt.ylabel('Ligand')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, f'top_ligands_{receptor}.png'))
    plt.close()
    print(f"Generated: top_ligands_{receptor}.png")

# 4. Scatter Plot: Affinity vs. Molecular Weight (colored by Receptor)
plt.figure(figsize=(10, 7))
sns.scatterplot(data=df, x='MW', y='Affinity_kcal_mol', hue='Receptor', alpha=0.7)
plt.title('Affinity vs. Molecular Weight')
plt.xlabel('Molecular Weight (Da)')
plt.ylabel('Affinity (kcal/mol)')
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, 'affinity_vs_mw.png'))
plt.close()
print("Generated: affinity_vs_mw.png")

# 5. Scatter Plot: pKi vs. LogP (colored by Receptor) - Lipinski-like view
plt.figure(figsize=(10, 7))
# Drop rows where pKi or LogP is NaN before plotting
plot_df = df.dropna(subset=['pKi', 'LogP']) 
if not plot_df.empty:
    sns.scatterplot(data=plot_df, x='LogP', y='pKi', hue='Receptor', alpha=0.7)
    plt.title('pKi vs. LogP')
    plt.xlabel('LogP')
    plt.ylabel('pKi')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, 'pki_vs_logp.png'))
    plt.close()
    print("Generated: pki_vs_logp.png")
else:
    print("Skipping pKi vs LogP plot (missing data).")

# 6. Scatter Plot: Ligand Efficiency (LE) vs Affinity
plt.figure(figsize=(10, 7))
plot_df = df.dropna(subset=['Affinity_kcal_mol', 'LE'])
if not plot_df.empty:
    sns.scatterplot(data=plot_df, x='Affinity_kcal_mol', y='LE', hue='Receptor', alpha=0.7)
    plt.title('Ligand Efficiency (LE) vs. Affinity')
    plt.xlabel('Affinity (kcal/mol)')
    plt.ylabel('LE (Affinity / NHA)')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, 'le_vs_affinity.png'))
    plt.close()
    print("Generated: le_vs_affinity.png")
else:
    print("Skipping LE vs Affinity plot (missing data).")
    
# 7. Scatter Plot: Surface-based Efficiency (SILE_SASA) vs pKi
plt.figure(figsize=(10, 7))
plot_df = df.dropna(subset=['pKi', 'SILE_SASA'])
if not plot_df.empty:
    sns.scatterplot(data=plot_df, x='pKi', y='SILE_SASA', hue='Receptor', alpha=0.7)
    plt.title('Surface Ligand Efficiency (SILE_SASA) vs. pKi')
    plt.xlabel('pKi')
    plt.ylabel('SILE_SASA (pKi / SASA)')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, 'sile_sasa_vs_pki.png'))
    plt.close()
    print("Generated: sile_sasa_vs_pki.png")
else:
    print("Skipping SILE_SASA vs pKi plot (missing data).")

# 8. Correlation Heatmap of Numerical Metrics
# Select only numeric columns for correlation calculation
numeric_df = df[numeric_cols].copy()
# Replace potential infinite values with NaN before calculating correlation
numeric_df.replace([np.inf, -np.inf], np.nan, inplace=True)
# Calculate correlation matrix (dropping columns/rows with all NaNs if any)
corr_matrix = numeric_df.corr()

# Plot heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', fmt=".2f", linewidths=.5, linecolor='black')
plt.title('Correlation Matrix of Calculated Metrics')
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, 'metrics_correlation_heatmap.png'))
plt.close()
print("Generated: metrics_correlation_heatmap.png")


# 9. Pair Plot for selected important metrics (can be slow for many points/metrics)
selected_metrics = ['Affinity_kcal_mol', 'pKi', 'MW', 'LogP', 'LE', 'SILE_SASA', 'QED']
# Check which selected metrics are actually present and numeric
available_metrics = [m for m in selected_metrics if m in df.columns and pd.api.types.is_numeric_dtype(df[m])]
if len(available_metrics) > 1:
    print("Generating Pair Plot (this might take a moment)...")
    # Use a sample if the dataframe is very large to speed it up
    sample_df = df if len(df) < 1000 else df.sample(1000) 
    pair_plot_df = sample_df[available_metrics + ['Receptor']].dropna() # Drop NaNs for plotting
    if not pair_plot_df.empty:
        sns.pairplot(pair_plot_df, hue='Receptor', diag_kind='kde', plot_kws={'alpha':0.6})
        plt.suptitle('Pair Plot of Selected Metrics', y=1.02) # Adjust title position
        plt.tight_layout()
        plt.savefig(os.path.join(PLOT_DIR, 'metrics_pairplot.png'))
        plt.close()
        print("Generated: metrics_pairplot.png")
    else:
        print("Skipping Pair Plot (no data after dropping NaNs).")
else:
    print("Skipping Pair Plot (not enough numeric metrics selected or available).")

print(f"\nAll plots saved to the '{PLOT_DIR}' directory.")