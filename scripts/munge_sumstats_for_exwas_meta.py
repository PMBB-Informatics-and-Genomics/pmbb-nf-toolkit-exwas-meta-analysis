import pandas as pd
import argparse as ap
import pandas as pd
from scipy.stats import norm
import numpy as np
import argparse as ap
import sys

TEST_ROWS = None

def make_arg_parser():
    parser = ap.ArgumentParser(description=".")

    # Required Args
    parser.add_argument('-p', '--pheno', required=True, help='Phenotype')
    parser.add_argument('-c', '--cohort', required=True, help='Cohort')
    parser.add_argument('-s', '--sumstats', required=True, help='Path to sumstats file')
    parser.add_argument('--colnames', required=True, help='File with column name mappings')

    # Non-Requird Args
    parser.add_argument('-o', '--outDir', default='./', help='Path to output directory. Default: current working directory')
    parser.add_argument('--singles', required=False, action='store_true', help='Use flag if munging single-variant tests')
    
    return parser

args = make_arg_parser().parse_args()
print(args)

# Unpack args into variables
cohort = args.cohort
pheno = args.pheno
ss_file = args.sumstats
colnames_file = args.colnames
is_singles = args.singles
out_dir = args.outDir

# Make expected output file names
output_name = f'{out_dir}/{cohort}.{pheno}.munged_meta_{"singles" if is_singles else "regions"}.txt.gz'
output_name = output_name.replace('//', '/')
print(output_name)

# Parse column names dict
colnames_rows = open(colnames_file).read().splitlines()
col_map = dict(zip([r.split('=')[1] for r in colnames_rows],
                   [r.split('=')[0] for r in colnames_rows]))

print(col_map)

# Read and Rename
df = pd.read_table(ss_file)
df = df.rename(columns=col_map)
df = df.rename(columns={'p_value': 'p_single' if is_singles else 'p_burden',
                        'beta': 'beta_single' if is_singles else 'beta_burden',
                        'se': 'se_single' if is_singles else 'se_burden'})
print(df)
print(df.columns)

# Do column checks for various conditions
if 'null' in df.columns:
    # Drop any columns that were mapped to null from the params
    df = df.drop(columns=[c for c in df.columns if 'null' in c])

# Check required info columns
info_cols = ['chr', 'pos', 'effect_allele', 'other_allele'] if is_singles else ['region', 'annot_group', 'max_maf']
if not np.all([c in df.columns for c in info_cols]):
    sys.exit(f'Cannot properly munge: missing info columns {[c for c in info_cols if c not in df.columns]}')

# Check effect columns (if provided)
check_effects = [f'{c}_{"single" if is_singles else "burden"}' for c in ['p', 'beta', 'se']]
if np.any([c in df.columns for c in check_effects]):
    # If any columns are present, all must be present
    if not np.all([c in df.columns for c in check_effects]):
        sys.exit(f'Cannot properly munge: missing {"singles" if is_singles else "burden"} effect columns {[c for c in check_effects if c not in df.columns]}')

# Fill missing values for sample sizes if not available
sample_cols = ['n', 'n_case', 'n_control']
for c in sample_cols:
    if c not in df.columns:
        df[c] = np.nan

# Write to output
df.to_csv(output_name, index=False, sep='\t', na_rep='NA')
