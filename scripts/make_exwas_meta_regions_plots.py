from manhattan_plot import ManhattanPlot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse as ap
import os
import sys

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
regions_sumstats = args.sumstats

# Read in full summary stats file
full_sumstats = pd.read_table(regions_sumstats)

# loop through all grouptest annotations and MAF value combinations
for ann in full_sumstats['annot_group'].unique():
    for maf in full_sumstats['max_maf'].unique():

        # Set up output file name and plot title
        maf = np.float64(maf)
        output_manhattan = f'{analysis}.{pheno}.{ann.replace(";","-")}.maf{str(maf).replace(".","-")}.regions.manhattan.png'
        output_qq = f'{analysis}.{pheno}.{ann.replace(";","-")}.maf{str(maf).replace(".","-")}.regions.qq.png'
        plot_title = f'SAIGE ExWAS Regions {analysis}: {pheno.replace("_", " ")}\nGroup = {ann}, MAF = {maf}'
        print(f'\n{plot_title}')

        # Instantiate Manhattan Plot object
        mp = ManhattanPlot(regions_sumstats, test_rows=None, title=plot_title)
        mp.load_data()

        # Replace values less than min float precision with min value
        mp.df = mp.df[(mp.df['annot_group'] == ann) & (mp.df['max_maf'] == maf)]
        # not every combination exists. Pass over iteration if DataFrame is empty

        if mp.df.empty:
            print(
                f"The combination of Group = {ann}, MAF = {maf} does not exist")
            continue
        
        # Clean data for Manhattan plot
        mp.clean_data(col_map={'chr': '#CHROM', 'pos_start': 'POS', 'gene_symbol': 'ID', 'p_burden_stouffer_meta': 'P'})
        # Protect against underflow
        mp.df['P'] = np.where(mp.df['P'].astype(float) < sys.float_info.min, sys.float_info.min, mp.df['P'].astype(float))
        mp.get_thinned_data()
        mp.thinned = mp.thinned.dropna(subset='P')

        # Dynamically compute the p-value threshold
        num_ind_tests = len(mp.df.ID.unique())
        p_thresh = 0.05 / num_ind_tests

        if len(mp.thinned) <= 1:
            continue

        # Set the p-value threshold to be less stringent if there are no significant hits
        if ~np.any(mp.thinned['P'] < p_thresh):
            keep_signals = min(10, len(mp.thinned))
            p_thresh = 10 ** -np.nanquantile(mp.thinned['ROUNDED_Y'], 1 - keep_signals / len(mp.thinned))

        print(p_thresh)
        mp.sig_line = p_thresh

        # Generate the Manhattan Plot
        mp.update_plotting_parameters(vertical=False, sig=p_thresh, annot_thresh=p_thresh, 
                                      merge_genes=False, ld_block=0)
        mp.full_plot(save=output_manhattan, rep_boost=False, legend_loc='top')
        print(f"Saved Manhattan plot to: {output_manhattan}")
        plt.clf()

        # Generate the QQ Plot
        mp.qq_plot(save=output_qq)
        print(f"Saved qq plot to: {output_qq}")
        plt.clf()
        
# plotting manifest
def plots_filename(row):
    # pheno = row['phenotype']
    ann = row['annot_group']
    maf = row['max_maf']
    # output_manhattan = f'Plots/{cohort}.{pheno}.{ann.replace(";","-")}.maf{str(maf).replace(".","-")}.regions.manhattan.png'
    output_manhattan = f'Regions_Plots/{analysis}.{pheno}.{ann.replace(";","-")}.maf{str(maf).replace(".","-")}.regions.manhattan.png'
    output_qq = f'Regions_Plots/{analysis}.{pheno}.{ann.replace(";","-")}.maf{str(maf).replace(".","-")}.regions.qq.png'
    return (output_manhattan, output_qq)
        
# regions_sumstats_df = pd.read_csv(regions_sumstats,sep='\t')
# regions_sumstats_df.rename(columns=columns_map_inv,inplace=True)
full_sumstats.insert(0,'cohort',analysis)
full_sumstats = full_sumstats.groupby(['cohort','phenotype', 'annot_group', 'max_maf']).first().reset_index()
full_sumstats['regions_sumstats'] = f"ExWAS_Meta/Sumstats/{regions_sumstats}"
full_sumstats['regions_results'] = "Summary/exwas_regions_meta_all_suggestive.csv"
full_sumstats[['regions_manhattan', 'regions_qq']] = full_sumstats.apply(lambda row: plots_filename(row), axis=1,result_type='expand')
output_file = f"{analysis}.{pheno}.exwas_regions.plots_manifest.csv"
full_sumstats.to_csv(output_file,index=False)