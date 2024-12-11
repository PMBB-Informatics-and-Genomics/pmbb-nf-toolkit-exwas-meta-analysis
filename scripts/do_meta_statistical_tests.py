import argparse as ap
import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.stats import chi2
import sys

def make_arg_parser():
    parser = ap.ArgumentParser()

    # All Required Args:
    parser.add_argument('--sumstats', '-s', nargs='+', help='List of summary stats files to be meta-analyzed')
    parser.add_argument('--singles', required=False, action='store_true', help='Use flag if doing single-variant tests')
    parser.add_argument('--analysis', '-a', required=True, help='The name of the meta-analysis group')
    parser.add_argument('--pheno', '-p', required=True, help='The name of the phenotype')
    parser.add_argument('--pthresh', required=True, help='The p-value threshold for filtering top hits', type=float)
    parser.add_argument('--geneFile', required=False, help='File with at least four columns: ')
    parser.add_argument('-o', '--outDir', default='./', help='Path to output directory. Default: current working directory')

    return parser

def perform_two_tailed_stouffers_meta(all_sumstats, group_info_cols):
    print('\nDoing Stouffer Test')
    all_sumstats = all_sumstats.copy()

    # Detect which p-value columns were provided
    all_p_cols = ['p_single', 'p_burden', 'skat_p_col', 'skato_p_col']
    use_p_cols = [c for c in all_p_cols if c in all_sumstats.columns]
    input_p_cols = [f'{c}_input' for c in use_p_cols]

    # Divide by two because the association test p-values are two-sided
    all_sumstats[input_p_cols] = all_sumstats[use_p_cols].values / 2

    # Set up an empty DataFrame to hold results
    meta_results = pd.DataFrame(dtype=float, index=pd.MultiIndex.from_frame(all_sumstats[group_info_cols].drop_duplicates()))

    # Iterate over p-value columns to be tested
    for p_col in input_p_cols:
        # Transform into two dimensions
        # Rows = Unique combinations of info columns
        # Columns = Cohorts
        p_vals_2D = all_sumstats.pivot(index=group_info_cols, columns='cohort', values=p_col)
        n_eff_2D = all_sumstats.pivot(index=group_info_cols, columns='cohort', values='n_eff')

        # Z Score = SQRT(N effective) * (P to Z normal distribution function)
        weighted_Z_2D = (n_eff_2D ** 0.5) * -norm.ppf(p_vals_2D)

        # Get the overall Z Score
        # Sum weighted Z and N_effective across each row
        # Then divide
        weighted_Z_sum = weighted_Z_2D.sum(axis=1)
        n_eff_sum = n_eff_2D.sum(axis=1)
        overall_Z = weighted_Z_sum / (n_eff_sum ** 0.5)

        # Get the P-value from the overall Z
        meta_p_col = p_col.replace('_input', '_stouffer_meta')
        meta_results.loc[overall_Z.index, meta_p_col] = 2 * norm.cdf(-overall_Z.astype(float))

        # Save the effective sample size and the number of studies included
        n_eff_col = p_col.replace('_input', '_stouffer_N_eff')
        meta_results[n_eff_col] = n_eff_sum.values
        n_studies_col = p_col.replace('_input', '_stouffer_N_studies')
        meta_results[n_studies_col] = p_vals_2D.count(axis=1).values

    return meta_results

def perform_inv_var_weighted_meta(all_sumstats, group_info_cols, is_singles):
    print('\nDoing Inverse Variance Weighted Test')
    all_sumstats = all_sumstats.copy()
    meta_results = pd.DataFrame(dtype=float, index=pd.MultiIndex.from_frame(all_sumstats[group_info_cols].drop_duplicates()))

    # Figure out which columns to use
    p_col = 'p_single' if is_singles else 'p_burden'
    beta_col = 'beta_single' if is_singles else 'beta_burden'
    se_col = 'se_single' if is_singles else 'se_burden'

    # Transform effect sizes, standard errors, and effective sample size into 2D
    # Rows = Unique combinations of info columns
    # Columns = Cohorts
    beta_2D = all_sumstats.pivot(index=group_info_cols, columns='cohort', values=beta_col)
    se_2D = all_sumstats.pivot(index=group_info_cols, columns='cohort', values=se_col)
    n_eff_2D = all_sumstats.pivot(index=group_info_cols, columns='cohort', values='n_eff').mask(np.isnan(beta_2D))
    
    # Weight the Effects by the Inverse of the Variance
    inv_var_2D = 1 / (se_2D ** 2)
    weighted_effects_2D = beta_2D * inv_var_2D

    # Take weighted average of the effects
    meta_results[beta_col + '_inv_var_meta'] = weighted_effects_2D.sum(axis=1) / inv_var_2D.sum(axis=1)
    meta_results[se_col + '_inv_var_meta'] = np.sqrt(1 / inv_var_2D.sum(axis=1))

    # Get Z score and P-Value
    Z = meta_results[beta_col + '_inv_var_meta'] / meta_results[se_col + '_inv_var_meta']
    meta_results.loc[Z.index, p_col + '_inv_var_meta'] = 2 * norm.cdf(-Z.abs())

    # Save number of samples and number of studies
    n_eff_col = 'N_eff_inv_var_meta'
    meta_results[n_eff_col] = n_eff_2D.sum(axis=1).values
    n_studies_col = 'N_studies_inv_var_meta'
    meta_results[n_studies_col] = beta_2D.count(axis=1).values

    return meta_results

def perform_chi_sq_heterogeneity_meta(all_sumstats, group_info_cols, is_singles):
    print('\nDoing Inverse Variance Weighted Test')

    # Copy the full summart stats, set up an empty one for results
    all_sumstats = all_sumstats.copy()
    meta_results = pd.DataFrame(dtype=float, index=pd.MultiIndex.from_frame(all_sumstats[group_info_cols].drop_duplicates()))

    # Figure out which columns to use
    p_col = 'p_single' if is_singles else 'p_burden'
    beta_col = 'beta_single' if is_singles else 'beta_burden'
    se_col = 'se_single' if is_singles else 'se_burden'

    # Transform the data into 2D DataFrames for computations
    # Rows = info columns
    # Columns = cohorts
    beta_2D = all_sumstats.pivot(index=group_info_cols, columns='cohort', values=beta_col)
    beta_2D = beta_2D[beta_2D.count(axis=1) > 1]
    se_2D = all_sumstats.pivot(index=group_info_cols, columns='cohort', values=se_col).loc[beta_2D.index]
    n_eff_2D = all_sumstats.pivot(index=group_info_cols, columns='cohort', values='n_eff').mask(np.isnan(beta_2D)).loc[beta_2D.index]
    
    # Weight the effects by the SQRT of the effective N
    weights_2D = n_eff_2D ** 0.5
    weighted_betas_2D = weights_2D * beta_2D

    # Take the weighted average of the effects
    sum_betas = weighted_betas_2D.sum(axis=1)
    sum_weights = weights_2D.sum(axis=1)
    meta_beta = sum_betas / sum_weights
    meta_beta_2D = np.broadcast_to(meta_beta, tuple(reversed(beta_2D.shape))).T

    # Compute the CHI2 test statistic
    deviation = weights_2D * ((beta_2D - meta_beta_2D) ** 2)
    sum_deviation = deviation.sum(axis=1)
    meta_results[p_col + '_chi2_stat'] = sum_deviation.reindex(meta_results.index)

    # Compute the CHI2 P-Value
    n_studies = beta_2D.sum(axis=1)
    chi2_pval = pd.Series(data=chi2.sf(sum_deviation.astype(float), n_studies-1), index=sum_deviation.index)
    meta_results[p_col + '_heterogeneity'] = chi2_pval.reindex(meta_results.index)

    return meta_results

# Parse script arguments
args = make_arg_parser().parse_args()
sumstats_list = args.sumstats
is_singles = args.singles
analysis = args.analysis
pheno = args.pheno
p_thresh = args.pthresh
out_dir = args.outDir

# Set up target output names from the beginning
output_name = f'{out_dir}{analysis}.{pheno}.meta_{"singles" if is_singles else "regions"}.txt.gz'
output_filtered_name = f'{out_dir}{analysis}.{pheno}.meta_{"singles" if is_singles else "regions"}.filtered.txt'
output_N_name = f'{out_dir}{analysis}.{pheno}.meta_{"singles" if is_singles else "regions"}.Ns.txt'
output_effects_name = f'{out_dir}{analysis}.{pheno}.meta_{"singles" if is_singles else "regions"}.effects.txt.gz'

# Info columns define unique tests to meta-analyze
info_cols = ['chr', 'pos', 'effect_allele', 'other_allele'] if is_singles else ['region', 'annot_group', 'max_maf']

# Iterate over input files to make long-form dataframe
summary_stats_dfs = []
for f in sumstats_list:
    print(f)
    temp = pd.read_table(f).drop_duplicates()
    cohort = f.replace('.munged_meta_singles.txt.gz', '').replace('.munged_meta_regions.txt.gz', '').replace(f'.{pheno}', '')
    temp['cohort'] = cohort
    temp['ss_filename'] = f
    test_counts = temp[info_cols].value_counts()
    summary_stats_dfs.append(temp)
summary_stats = pd.concat(summary_stats_dfs)

# Compute Effective sample size for binary phenotypes, keep N for quantitative
compute_n_eff = 2 / ((1 / summary_stats['n_control']) + (1 / summary_stats['n_case']))
compute_n_eff.update(summary_stats['n'])
summary_stats['n_eff'] = compute_n_eff

merge_results = [] # This test will collect results

# Always perform Stouffer tests on P-values
stouffer_results = perform_two_tailed_stouffers_meta(summary_stats, info_cols)
merge_results.append(stouffer_results)

# If effect sizes are available, do inverse variance weighted and heterogeneity
if 'beta_single' in summary_stats.columns or 'beta_burden' in summary_stats.columns:
    inv_var_results = perform_inv_var_weighted_meta(summary_stats, info_cols, is_singles)
    merge_results.append(inv_var_results)

    het_chi2_results = perform_chi_sq_heterogeneity_meta(summary_stats, info_cols, is_singles)
    merge_results.append(het_chi2_results)

# Combine results together
all_meta_results = pd.concat(merge_results, axis=1)

all_meta_results.insert(0, 'analysis', analysis)
all_meta_results.insert(1, 'phenotype', pheno)
print(all_meta_results)

# Insert the gene region information if needed,
# sort results by genomic position
if not is_singles:
    gene_file = pd.read_table(args.geneFile, index_col='gene_id')
    gene_id_order = all_meta_results.index.get_level_values('region')
    all_meta_results.insert(0, 'chr', gene_file.reindex(gene_id_order)['chromosome'].values)
    all_meta_results.insert(1, 'pos_start', gene_file.reindex(gene_id_order)['seq_region_start'].values)
    all_meta_results.insert(2, 'pos_stop', gene_file.reindex(gene_id_order)['seq_region_end'].values)
    all_meta_results.insert(3, 'gene_symbol', gene_file.reindex(gene_id_order)['gene_symbol'].values)

    all_meta_results = all_meta_results.sort_values(by=['chr', 'pos_start', 'annot_group', 'max_maf'])
else:
    all_meta_results = all_meta_results.sort_values(by=['chr', 'pos'])

# Write full output
all_meta_results.to_csv(output_name, sep='\t', na_rep='NA')

# Write filtered output
p_cols = [c for c in all_meta_results.columns if 'p_' in c and 'chi' not in c and 'heterogeneity' not in c]
n_study_cols = [c for c in all_meta_results.columns if 'N_studies' in c]

filtered_out = all_meta_results[all_meta_results[p_cols].fillna(1).min(axis=1) <= p_thresh].copy()
filtered_out = filtered_out[filtered_out[n_study_cols].max(axis=1) > 1].copy()
filtered_out.to_csv(output_filtered_name, sep='\t', na_rep='NA')

# Get max sample sizes for all cohort/phenotype tests
N_df = pd.Series(index=['PHENO', 'ANALYSIS', 'N_Samples'], dtype=object)
N_df['PHENO'] = pheno
N_df['ANALYSIS'] = analysis

n_cols = [c for c in all_meta_results.columns if 'N_eff' in c]
max_meta_N = all_meta_results[n_cols].max().max()

N_df['N_Samples'] = max_meta_N
N_df.index.name = 'col_name'
N_df.name = 'value'
N_df.to_csv(output_N_name, sep='\t')

# Keep a record of the individual study statistics used for each test
save_cols = ['p_single', 'beta_single', 'se_single', 'p_burden', 'beta_burden', 'se_burden', 'skat_p_col', 'skato_p_col']
save_cols = [c for c in save_cols if c in summary_stats.columns]

save_dfs = []
for col_name in save_cols:
    col_2D = summary_stats.pivot(index=info_cols, columns='cohort', values=col_name)
    col_2D.columns = [f'{col_name}_{c}' for c in col_2D.columns]
    save_dfs.append(col_2D)

save_df = pd.concat(save_dfs, axis=1)
save_df.index.names = info_cols

save_df.to_csv(output_effects_name, sep='\t', na_rep='NA')
