"""
Creating and Sorting Matrix for Cell-Cell Communication Analysis
---------------------------------------------------------------

This script processes raw single-cell expression matrices and prepares 
a clean `final_matrix.tsv` for downstream tools like CellCall.

Main Steps:
1. Merge raw expression matrices from multiple input files.
2. Rename cells based on their assigned class (from Excel metadata).
3. Filter out unwanted cell types (T_ and LD_).
4. Standardize SSC (spermatogonial stem cell) column naming.
5. Save processed matrices (`all_merged_cells.tsv` and `final_matrix.tsv`).
6. Report summary statistics (missing cells, removed columns, etc.).

Inputs:
- cell_classes.xlsx : Mapping of individual cell IDs to class names.
- matrix_folder     : Folder containing raw .txt gene expression matrices.

Outputs:
- all_merged_cells.tsv : Merged matrix of all raw cells.
- final_matrix.tsv     : Cleaned and standardized matrix for analysis.

Author: <Your Name>
Date: <Insert Date>
"""

import os
import pandas as pd
from collections import defaultdict


# ========================
# --- Configuration ------
# ========================
class_file = '/mnt/d_manual/computional biology/machine learning for biology/data/cell_classes.xlsx'
matrix_folder = '/mnt/d_manual/computional biology/machine learning for biology/data/GSE106487_RAW_ testicales'

merged_output = 'all_merged_cells.tsv'
filtered_output = 'final_matrix.tsv'


# =======================================================
# Step 1: Load cell class information (Excel annotation)
# =======================================================
class_df = pd.read_excel(class_file)
assert 'class' in class_df.columns and 'cell' in class_df.columns, \
    "Excel file must contain 'class' and 'cell' columns."

# Map original cell names → class names
cell_to_class = dict(zip(class_df['cell'], class_df['class']))


# ==========================================================
# Step 2: Merge all raw expression matrices by gene (row index)
# ==========================================================
merged_df = None

for filename in os.listdir(matrix_folder):
    if not filename.endswith('.txt'):
        continue

    path = os.path.join(matrix_folder, filename)
    print(f"📥 Merging {path}")

    try:
        df = pd.read_csv(path, sep='\t', index_col=0)
    except Exception as e:
        print(f"⚠️ Could not read {filename}: {e}")
        continue

    # Merge across gene names
    if merged_df is None:
        merged_df = df
    else:
        merged_df = merged_df.join(df, how='outer')

# Save merged matrix
merged_df.to_csv(merged_output, sep='\t')
print(f"\n✅ Merged matrix saved to: {merged_output}")


# ========================================================
# Step 3: Rename cells based on class & ensure uniqueness
# ========================================================
class_counter = defaultdict(int)
renamed_df = pd.DataFrame(index=merged_df.index)
matched_cells = set()

for cell_name_orig, class_name in cell_to_class.items():
    # Handle cases where pandas prefixes 'X' to numeric names
    candidates = [cell_name_orig, f"X{cell_name_orig}"]

    matched = None
    for candidate in candidates:
        if candidate in merged_df.columns:
            matched = candidate
            break

    if matched:
        class_counter[class_name] += 1
        new_col_name = f"{class_name}_{class_counter[class_name]}"
        renamed_df[new_col_name] = merged_df[matched]
        matched_cells.add(cell_name_orig)

# Save renamed and filtered matrix
renamed_df.to_csv(filtered_output, sep='\t')
print(f"✅ Renamed & filtered matrix saved to: {filtered_output}")


# ========================================================
# Step 4: Report missing cells
# ========================================================
expected_cells = set(class_df['cell'])
missing_cells = expected_cells - matched_cells

print(f"\n📊 Total expected cells: {len(expected_cells)}")
print(f"✅ Found and renamed: {len(matched_cells)}")
print(f"❌ Missing: {len(missing_cells)}")

if missing_cells:
    print("🔍 Missing cell names:")
    for cell in sorted(missing_cells):
        print(f"  - {cell}")


# ===========================================================
# Extra Processing: Sorting, filtering, SSC standardization
# ===========================================================
class_file = '/mnt/d/computional biology/machine learning for biology/data/cell_classes.xlsx'
input_matrix = 'input_matrix.tsv'   # Input matrix file (already merged)
output_matrix = 'final_matrix.tsv'  # Output cleaned file

# --- Reload class info ---
class_df = pd.read_excel(class_file)

# --- Load matrix ---
try:
    df = pd.read_csv(input_matrix, sep='\t', index_col=0)
    print(f"✅ Loaded input matrix: {df.shape}")
except Exception as e:
    print(f"Failed to load {input_matrix}: {e}")
    exit()

# --- Map original → standardized names ---
cell_name_mapping = {}
class_counters = defaultdict(int)

for _, row in class_df.iterrows():
    cell_class = row['class']
    cell_name = row['cell']
    class_counters[cell_class] += 1
    new_name = f"{cell_class}_{class_counters[cell_class]}"
    cell_name_mapping[cell_name] = new_name

# --- Apply renaming ---
columns_to_keep = []
new_column_names = []

for col in df.columns:
    if col in cell_name_mapping:
        columns_to_keep.append(col)
        new_column_names.append(cell_name_mapping[col])
    else:
        columns_to_keep.append(col)
        new_column_names.append(col)

renamed_df = df[columns_to_keep].copy()
renamed_df.columns = new_column_names
print(f"📊 After initial renaming: {renamed_df.shape}")

# --- Remove unwanted T_ and LD_ cells ---
columns_to_keep_final = []
for col in renamed_df.columns:
    if col.startswith('T_') or col.startswith('LD_'):
        print(f"🗑️ Deleting column: {col}")
        continue
    columns_to_keep_final.append(col)

filtered_df = renamed_df[columns_to_keep_final].copy()
print(f"📊 After filtering T_ and LD_: {filtered_df.shape}")

# --- Standardize SSC columns ---
ssc_columns = [col for col in filtered_df.columns if col.startswith('SSC')]
non_ssc_columns = [col for col in filtered_df.columns if not col.startswith('SSC')]

print(f"🔍 Found {len(ssc_columns)} SSC columns to combine")

if ssc_columns:
    ssc_rename_map = {}
    for i, col in enumerate(ssc_columns, 1):
        ssc_rename_map[col] = f"SSC_{i}"
        print(f"🔄 SSC renaming: {col} -> SSC_{i}")
    final_df = filtered_df.rename(columns=ssc_rename_map)
else:
    final_df = filtered_df

print(f"📊 Final matrix shape: {final_df.shape}")

# --- Save processed matrix ---
try:
    final_df.to_csv(output_matrix, sep='\t')
    print(f"✅ Saved processed matrix to: {output_matrix}")
except Exception as e:
    print(f"Failed to save {output_matrix}: {e}")
    exit()

# --- Final report ---
original_cols = len(df.columns)
final_cols = len(final_df.columns)
deleted_cols = original_cols - final_cols

print(f"\n📋 Processing Summary:")
print(f"📊 Original columns: {original_cols}")
print(f"🗑️ Deleted columns: {deleted_cols}")
print(f"✅ Final columns: {final_cols}")
print(f"🧬 SSC columns unified: {len(ssc_columns)}")
