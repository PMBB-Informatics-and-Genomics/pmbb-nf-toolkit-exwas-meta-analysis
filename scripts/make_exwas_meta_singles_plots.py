from manhattan_plot import ManhattanPlot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse as ap
import sys
import os

def make_arg_parser():
    parser = ap.ArgumentParser(description=".")

    parser.add_argument('-p', '--phenotype', required=True, help='phenotype')
    parser.add_argument('-a', '--analysis', required=True, help='cohort')
    parser.add_argument('-s', '--sumstats', required=True)

    return parser

# Parse command line arguments
args = make_arg_parser().parse_args()
analysis = args.analysis
pheno = args.phenotype
singles_sumstats = args.sumstats

# Set up output file name and plot title
output_manhattan = f'{analysis}.{pheno}.singles.manhattan_vertical.png'
output_qq = f'{analysis}.{pheno}.singles.qq.png'
plot_title = f'SAIGE ExWAS Singles {analysis}: {pheno.replace("_", " ")}'

# Construct Manhattan plot object
mp = ManhattanPlot(singles_sumstats, test_rows=None, title=plot_title)
mp.load_data()

# Replace values less than min float precision with min value
mp.df['ID'] = 'chr' + mp.df['chr'].astype(str) + ':' + mp.df['pos'].astype(str) + ':' \
                + mp.df['other_allele'] + ':' + mp.df['effect_allele']

# Clean data for Manhattan plotting procedure
mp.clean_data(col_map={'chr': '#CHROM', 'p_single_stouffer_meta': 'P', 'pos': 'POS'})
# Protect from underflow errors
mp.df['P'] = np.where(mp.df['P'].astype(float) < sys.float_info.min, sys.float_info.min, mp.df['P'].astype(float))
mp.get_thinned_data()
mp.thinned = mp.thinned.dropna(subset='P')

# Dynamically figure out p-value threshold for plotting
num_ind_tests = len(mp.df.ID.unique())
if num_ind_tests == 0:
    print('All Single Assoc Tests were Rare')
    open(output_manhattan, 'w+').write('All Single Assoc Tests Were Ultrarare')
    open(output_qq, 'w+').write('All Single Assoc Tests Were Ultrarare')
    open(output_qq.replace('.png', '.csv'), 'w+').write('All Single Assoc Tests Were Ultrarare')
    sys.exit()

p_thresh = 0.05 / num_ind_tests
if ~np.any(mp.thinned['P'] < p_thresh):
    # If there are no p-values below the threshold, make it less stringent
    keep_signals = min(10, len(mp.thinned))
    p_thresh = 10 ** -np.nanquantile(mp.thinned['ROUNDED_Y'], 1 - keep_signals / len(mp.thinned))

print(p_thresh)
mp.sig_line = p_thresh

# Replace underscores with colons to avoid confusing the LaTeX formating
mp.thinned['ID'] = mp.thinned['ID'].str.replace('_', ':')

# Generate the Manhattan plot and QQ plot
mp.update_plotting_parameters(vertical=True, sig=p_thresh, annot_thresh=p_thresh, merge_genes=True)
mp.full_plot(save=output_manhattan, rep_boost=False, extra_cols={'beta_single_inv_var_meta': 'BETA'}, number_cols=['BETA'], keep_chr_pos=False)
plt.clf()

print(f"Saved Manhattan plot to: {output_manhattan}")

mp.qq_plot(save=output_qq)
print(f"Saved qq plot to: {output_qq}")

# plotting manifest
def plots_filename(row):
    manhattan_file = f'Singles_Plots/{output_manhattan}'
    qq_file = f'Singles_Plots/{output_qq}'
    return (manhattan_file, qq_file)

   
singles_sumstats_df = pd.read_csv(singles_sumstats,sep='\t')
singles_sumstats_df.insert(0,'cohort',analysis)
singles_sumstats_df = singles_sumstats_df.groupby(['cohort','phenotype']).first().reset_index()
singles_sumstats_df['singles_results'] = "Summary/exwas_singles_meta_all_suggestive.csv"
singles_sumstats_df['singles_sumstats'] = f"ExWAS_Meta/Sumstats/{singles_sumstats}" 
singles_sumstats_df[['singles_manhattan', 'singles_qq']] = singles_sumstats_df.apply(lambda row: plots_filename(row), axis=1,result_type='expand')
output_file = f"{analysis}.{pheno}.exwas_singles.plots_manifest.csv"
singles_sumstats_df.to_csv(output_file,index=False)