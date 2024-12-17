nextflow.enable.dsl = 2

// log info
log.info """\
NEXTFLOW - DSL2 - ExWAS Meta-Analysis - P I P E L I N E
============================================
run_as            :  ${workflow.commandLine}
run_location      :  ${launchDir}
started_at        :  ${workflow.start}
container         :  ${workflow.containerEngine}:${workflow.container}
python            :  ${params.my_python}

Phenotype Lists
============================================
binary phenos     :  ${params.bin_pheno_list}
quant. phenos     :  ${params.quant_pheno_list}

Input Information
============================================
regions_sumstats_suffix   :  ${params.regions_sumstats_suffix}
singles_sumstats_suffix   :  ${params.singles_sumstats_suffix}

Meta-Analysis Parameters
============================================
Stouffer P-Values will be used with
Inverse-Variance Weighted Effect Size Estimates

Post-Processing Parameters
============================================
p thresh          :  ${params.p_cutoff_summarize}

Meta-Analysis Cohort Scheme
============================================
${params.analyses.collect().join('\n').replace('=', '\t:\t')}

"""
.stripIndent()

// This is the main workflow block which will be executed
// by default with nextflow run exwas_meta_analysis.nf
workflow {
    cohort_pheno_sumstats = EXWAS_META_SETUP()
    exwas_meta_sumstats = EXWAS_META(cohort_pheno_sumstats)
}

// biofilter_wrapper.nf should be linked into the same directory
// as this .nf file
// we don't need this for now...TBD on how to annotate single variants with genes
// include { BIOFILTER_POSITIONS } from './biofilter_wrapper.nf'

// We do this a couple of times, so wrap in a function to avoid redundancy
// (redundancy is evil)
import groovyx.gpars.dataflow.DataflowBroadcast

DataflowBroadcast get_pheno_channel(Map params) {
    all_phenos = []
    all_phenos.addAll(params.bin_pheno_list)
    all_phenos.addAll(params.quant_pheno_list)
    pheno = Channel.fromList(all_phenos)
    return pheno
}

workflow EXWAS_META_SETUP {
    main:
        /*
        This is the code chunk for starting with The ExWAS Meta-Analysis
        If you have summary statistics from somewhere else,
        you can pass them directly
        */

        // Make phenotype Channel
        pheno = get_pheno_channel(params)

        // Make cohort Channel by iterating over analysis groups from config
        // to get unique cohorts
        all_unique_cohorts = []
        params.analyses.each { analysis, cohort_list -> all_unique_cohorts.addAll(cohort_list) }
        all_unique_cohorts = all_unique_cohorts.unique()

        // Make cohort x pheno combinations
        cohort_pheno = Channel.fromList(all_unique_cohorts).combine(pheno)

        // Map cohort x pheno combinations to file paths for regions AND singles
        cohort_pheno_regions_sumstats = cohort_pheno.map { c, p ->
            new Tuple(c, p, "${launchDir}/${c}/Sumstats/${c}.${p}${params.regions_sumstats_suffix}")
        }
        cohort_pheno_singles_sumstats = cohort_pheno.map { c, p ->
            new Tuple(c, p, "${launchDir}/${c}/Sumstats/${c}.${p}${params.singles_sumstats_suffix}")
        }
    emit:
        cohort_pheno_regions_sumstats // channel of tuples with 3 items: cohort, pheno, file
        cohort_pheno_singles_sumstats // channel of tuples with 3 items: cohort, pheno, file
}

workflow EXWAS_META {
    take:
        cohort_pheno_regions_sumstats // channel of tuples with 3 items: cohort, pheno, file
        cohort_pheno_singles_sumstats // channel of tuples with 3 items: cohort, pheno, file
    main:

        // This first section does a lot of Channel Logic to organize
        // our summary stats files and cohorts
        // -------------------------------------------------

        // Get all of the phenotypes in one list
        pheno = get_pheno_channel(params)

        // Get a List of Tuples of ALL (cohort, analysis) combinations from params
        // Examples: (AFR_ALL, AFR_EUR), (EUR_ALL, AFR_EUR)
        // Examples: (AFR_M, ALL_M), (AMR_M, ALL_M), (EUR_M, ALL_M)
        all_analyses_cohorts = []
        params.analyses.each { analysis, cohort_list ->
            cohort_list.each { cohort ->
                all_analyses_cohorts.add(new Tuple(cohort, analysis))
            }
        }

        /*
        Goal: get a channel to join with the munge output channel.
            create Channel of (cohort, analysis) pairs,
            combine (cohort, analysis) pairs with phenotype, -> (cohort, analysis, pheno)
            reorder Tuples to be (cohort, pheno, analysis),
            group tuples by cohort and pheno to get (cohort, pheno, analysis_list)
        */
        all_analyses_cohorts = Channel.fromList(all_analyses_cohorts).combine(pheno).map {
            cohort, analysis, pheno ->
            new Tuple(cohort, pheno, analysis)
            }.groupTuple(by: [0, 1])

        // -----------------------------------------------------------------
        // This next section does some organization to see which tests to perform
        // -----------------------------------------------------------------
        // Possible Meta-Analysis Tests

        // Defining script locations relative to the primary workflow file
        exwas_meta_munge_script = "${projectDir}/scripts/munge_sumstats_for_exwas_meta.py"
        meta_analysis_script = "${projectDir}/scripts/do_meta_statistical_tests.py"
        meta_singles_plot_script = "${projectDir}/scripts/make_exwas_meta_singles_plots.py"
        meta_regions_plot_script = "${projectDir}/scripts/make_exwas_meta_regions_plots.py"

        // -----------------------------------------------------------------
        // Now, we will actually handle calling processes for the meta-analyses
        // -----------------------------------------------------------------

        // Make sure we have regions info to do the meta-analyses
        singles_sumstats_drop = null
        regions_sumstats_drop = null

        if (!(params.burden_cols == null && params.skat_p_col == null && params.skato_p_col == null)) {
            // Iterate over input channel to see which cohort/phenotype files exist
            regions_sumstats_keep = cohort_pheno_regions_sumstats.filter { c, p, sumstats -> \
                    file(sumstats).exists() && (new File(sumstats)).length() > 100
            }
            regions_sumstats_drop = cohort_pheno_regions_sumstats.filter { c, p, sumstats -> \
                    ! (file(sumstats).exists() && (new File(sumstats)).length() > 100)
                }.map { c, p, sumstats -> \
                    new Tuple(c, p, 'regions')
            }

            // Pass each individual cohort/phenotype summary stats file to be munged
            regions_col_info = params.burden_cols + params.regions_info_cols
            regions_col_info['p_skat'] = params.skat_p_col
            regions_col_info['p_skato'] = params.skato_p_col
            regions_munge_output = munge_regions_sumstats_file(regions_sumstats_keep, regions_col_info, exwas_meta_munge_script)

            /*
            Goal: get a Channel with (cohort, pheno, analysis, sumstats file)
                join (cohort, pheno, analysis_list) with (cohort, pheno, sumstats_file)
                to get (cohort, pheno, analysis_list, sumstats_file)
                transpose to get new Channel with (cohort, pheno, analysis, file path)
            */
            analysis_cohort_regions_files = all_analyses_cohorts.join(regions_munge_output, by: [0, 1]).transpose()

            /*
            Take (cohort, pheno, analysis, file path),
                group by pheno and analysis to get (cohort_list, pheno, analysis, file_list),
                reorder to get new Channel with (analysis, pheno, cohort_list, file_list)
            */
            analysis_regions_inputs = analysis_cohort_regions_files.groupTuple(by: [1, 2]).map { cl, p, a, fl -> new Tuple(a, p, cl, fl) }

            gene_location_file = params['gene_location_file']

            // Call the meta-analysis tests for regions
            (meta_regions_sumstats, meta_regions_filtered, meta_regions_Ns, meta_regions_effects) = do_regions_tests(
                analysis_regions_inputs,
                meta_analysis_script,
                params['p_cutoff_summarize'],
                gene_location_file
            )

            regions_filtered_sumstats_list = meta_regions_filtered.map { analysis, pheno, sumstats -> sumstats }.collect()

            // Grab filtered summary stats for top hits table
            make_regions_summary_table(regions_filtered_sumstats_list)

            // Make the region-based Manhattan plots
            regions_plots = plot_regions_meta(meta_regions_sumstats, meta_regions_plot_script)

            // Here we collect all the individual manifest files and concatenate them
            // take the 4 input tuples of cohort, pheno, pngs and csvs, extract csvs, filter on manifest
            regions_csvs = regions_plots.map { cohort, pheno, pngs, csvs -> new Tuple(csvs) }.transpose().filter { it.name =~ /.*manifest.csv/ }.collect()
            exwas_regions = 'exwas_regions'
            regions_manifest = collect_exwas_regions_plots(exwas_regions, regions_csvs)
        }

        // Make sure we have singles info to do the meta-analyses
        if (params.singles_effect_cols != null) {
            // Iterate over input channel to see which cohort/phenotype files exist
            singles_sumstats_keep = cohort_pheno_singles_sumstats.filter { c, p, sumstats -> file(sumstats).exists() }
            singles_sumstats_drop = cohort_pheno_singles_sumstats.filter { c, p, sumstats -> \
                    !file(sumstats).exists()
                }.map { c, p, sumstats -> \
                    new Tuple(c, p, 'singles')
            }

            // Pass each individual cohort/phenotype summary stats file to be munged
            singles_col_info = params.singles_effect_cols + params.singles_info_cols
            singles_munge_output = munge_singles_sumstats_file(singles_sumstats_keep, singles_col_info, exwas_meta_munge_script)

            /*
            Goal: get a Channel with (cohort, pheno, analysis, sumstats file)
                join (cohort, pheno, analysis_list) with (cohort, pheno, sumstats_file)
                to get (cohort, pheno, analysis_list, sumstats_file)
                transpose to get new Channel with (cohort, pheno, analysis, file path)
            */
            analysis_cohort_singles_files = all_analyses_cohorts.join(singles_munge_output, by: [0, 1]).transpose()

            /*
            Take (cohort, pheno, analysis, file path),
                group by pheno and analysis to get (cohort_list, pheno, analysis, file_list),
                reorder to get new Channel with (analysis, pheno, cohort_list, file_list)
            */
            analysis_singles_inputs = analysis_cohort_singles_files.groupTuple(by: [1, 2]).map { cl, p, a, fl -> new Tuple(a, p, cl, fl) }

            // Call the meta-analysis tests for singles
            (meta_singles_sumstats, meta_singles_filtered, meta_singles_Ns, meta_singles_effects) = do_singles_tests(
                analysis_singles_inputs,
                meta_analysis_script,
                params['p_cutoff_summarize']
            )

            singles_filtered_sumstats_list = meta_singles_filtered.map { analysis, pheno, sumstats -> sumstats }.collect()

            // Grab filtered summary stats for top hits table
            make_singles_summary_table(singles_filtered_sumstats_list)

            // Generate Manhattan plots for the rare single-variant tests
            singles_plots = plot_singles_meta(meta_singles_sumstats, meta_singles_plot_script)

            // Here we collect all the individual manifest files and concatenate them
            // take the 4 input tuples of cohort, pheno, pngs and csvs, filter on manifest
            singles_csvs = singles_plots.map { cohort, pheno, pngs, csvs -> new Tuple(csvs) }.transpose().filter { it.name =~ /.*manifest.csv/ }.collect()
            exwas_singles = 'exwas_singles'
            singles_manifest = collect_exwas_singles_plots(exwas_singles, singles_csvs)
        }

        // If any files did not exist or were too small, they were dropped and will be logged here.
        if (singles_sumstats_drop != null && singles_sumstats_drop != null) {
            write_dropped_file(regions_sumstats_drop.concat(singles_sumstats_drop).collect(flat: false))
            sample_sizes_list = meta_regions_Ns.concat(meta_singles_Ns).map { analysis, pheno, info -> info }.collect()
        } else if (singles_sumstats_drop == null) {
            write_dropped_file(regions_sumstats_drop.collect(flat: false))
            sample_sizes_list = meta_regions_Ns.map { analysis, pheno, info -> info }.collect()
        } else {
            write_dropped_file(singles_sumstats_drop.collect(flat: false))
            sample_sizes_list = meta_singles_Ns.map { analysis, pheno, info -> info }.collect()
        }

        // Here we can collect all the final meta-analysis sizes (helpful later)
        sample_size_table = make_analysis_size_table(sample_sizes_list)
}

// If any of the expected summary stats files didn't exist, write a log
process write_dropped_file {
    publishDir "${launchDir}/ExWAS_Meta/"

    machineType 'n2-standard-4'

    input:
        val cohort_pheno_tuples
    output:
        path('sumstats_dropped.txt')
    shell:
        """
        echo "${cohort_pheno_tuples}" \
          | sed 's|], \\[|\\n|g' \
          | sed 's|\\[\\[||g' \
          | sed 's|]]||g' \
          | sort \
          > sumstats_dropped.txt
        """
}

process munge_regions_sumstats_file {
    publishDir "${launchDir}/ExWAS_Meta/Input_Munged/"

    machineType 'n2-standard-4'

    input:
        tuple val(cohort), val(pheno), path(sumstats_file)
        val column_params
        path munge_script
    output:
        tuple val(cohort), val(pheno), path("${cohort}.${pheno}.munged_meta_regions.txt.gz")
    shell:
        """
        echo "${column_params.collect().join('\n')}" > colnames.txt
        ${params.my_python} ${munge_script} \
          -c ${cohort} \
          -p ${pheno} \
          -s ${sumstats_file} \
          --colnames colnames.txt
        """
    stub:
        """
        touch ${cohort}.${pheno}.munged_meta_regions.txt.gz
        """
}

process munge_singles_sumstats_file {
    publishDir "${launchDir}/ExWAS_Meta/Input_Munged/"

    machineType 'n2-standard-4'

    input:
        tuple val(cohort), val(pheno), path(sumstats_file)
        val column_params
        path munge_script
    output:
        tuple val(cohort), val(pheno), path("${cohort}.${pheno}.munged_meta_singles.txt.gz")
    shell:
        """
        echo "${column_params.collect().join('\n')}" > colnames.txt
        ${params.my_python} ${munge_script} \
          --singles \
          -c ${cohort} \
          -p ${pheno} \
          -s ${sumstats_file} \
          --colnames colnames.txt
        """
    stub:
        """
        touch ${cohort}.${pheno}.munged_meta_singles.txt.gz
        """
}

process do_regions_tests {
    publishDir "${launchDir}/ExWAS_Meta/Sumstats/"

    machineType 'n2-standard-4'

    input:
        tuple val(analysis), val(pheno), val(cohort_list), path(sumstats_file_list)
        path meta_test_python_script
        val params_p_thresh
        path regions_info_table
    output:
        tuple val(analysis), val(pheno), path("${analysis}.${pheno}.meta_regions.txt.gz")
        tuple val(analysis), val(pheno), path("${analysis}.${pheno}.meta_regions.filtered.txt")
        tuple val(analysis), val(pheno), path("${analysis}.${pheno}.meta_regions.Ns.txt")
        tuple val(analysis), val(pheno), path("${analysis}.${pheno}.meta_regions.effects.txt.gz")
    shell:
        """
        ${params.my_python} ${meta_test_python_script} \
          --analysis ${analysis} \
          --pheno ${pheno} \
          --sumstats ${sumstats_file_list.join(' ')} \
          --pthresh ${params_p_thresh} \
          --geneFile ${regions_info_table}
        """
    stub:
        """
        touch ${analysis}.${pheno}.meta_regions.txt.gz
        touch ${analysis}.${pheno}.meta_regions.filtered.txt
        touch ${analysis}.${pheno}.meta_regions.Ns.txt
        touch ${analysis}.${pheno}.meta_regions.effects.txt.gz
        """
}

process do_singles_tests {
    publishDir "${launchDir}/ExWAS_Meta/Sumstats/"

    machineType 'n2-standard-4'

    input:
        tuple val(analysis), val(pheno), val(cohort_list), path(sumstats_file_list)
        path meta_test_python_script
        val params_p_thresh
    output:
        tuple val(analysis), val(pheno), path("${analysis}.${pheno}.meta_singles.txt.gz")
        tuple val(analysis), val(pheno), path("${analysis}.${pheno}.meta_singles.filtered.txt")
        tuple val(analysis), val(pheno), path("${analysis}.${pheno}.meta_singles.Ns.txt")
        tuple val(analysis), val(pheno), path("${analysis}.${pheno}.meta_singles.effects.txt.gz")
    shell:
        """
        ${params.my_python} ${meta_test_python_script} \
          --analysis ${analysis} \
          --pheno ${pheno} \
          --sumstats ${sumstats_file_list.join(' ')} \
          --pthresh ${params_p_thresh} \
          --singles
        """
    stub:
        """
        touch ${analysis}.${pheno}.meta_singles.txt.gz
        touch ${analysis}.${pheno}.meta_singles.filtered.txt
        touch ${analysis}.${pheno}.meta_singles.Ns.txt
        touch ${analysis}.${pheno}.meta_singles.effects.txt.gz
        """
}

process plot_regions_meta {
    publishDir "${launchDir}/Regions_Plots"

    machineType 'n2-standard-4'

    label 'safe_to_skip'

    input:
        tuple val(analysis), val(pheno), path(sumstats)
        path exwas_meta_regions_plot_script
    output:
        tuple val(analysis), val(pheno), path("${analysis}.${pheno}*.png"), path("${analysis}.${pheno}.exwas_regions.plots_manifest.csv")
    shell:
        """
        ${params.my_python} ${exwas_meta_regions_plot_script} \
          -a ${analysis} \
          -p ${pheno} \
          -s ${sumstats} \
        """
    stub:
        """
        touch ${analysis}.${pheno}.regions.stub.png
        touch ${analysis}.${pheno}.exwas_regions.plots_manifest.csv
        """
}

process plot_singles_meta {
    publishDir "${launchDir}/Singles_Plots"

    machineType 'n2-standard-4'

    label 'safe_to_skip'

    input:
        tuple val(analysis), val(pheno), path(sumstats)
        path exwas_meta_singles_plot_script
    output:
        tuple val(analysis), val(pheno), path("${analysis}.${pheno}*.png"), path("${analysis}.${pheno}.exwas_singles.plots_manifest.csv")
    shell:
        """
        ${params.my_python} ${exwas_meta_singles_plot_script} \
          -a ${analysis} \
          -p ${pheno} \
          -s ${sumstats} \
        """
    stub:
        """
        touch ${analysis}.${pheno}.singles.stub.png
        touch ${analysis}.${pheno}.exwas_singles.plots_manifest.csv
        """
}

process make_regions_summary_table {
    publishDir "${launchDir}/Summary/"

    machineType 'n2-standard-4'

    input:
        path all_filtered_sumstats
    output:
        path('exwas_regions_meta_all_suggestive.csv')
    script:
        """
        #! ${params.my_python}

        import pandas as pd
        dfs = []
        input_list = '${all_filtered_sumstats.join(' ')}'.split()
        output = "exwas_regions_meta_all_suggestive.csv"

        # concatenate all results files together
        for f in input_list:
            dfs.append(pd.read_table(f))
        df = pd.concat(dfs).sort_values(by=['region', 'annot_group', 'max_maf'])

        # only present results in the top hits table with at least two studies
        n_study_cols = [c for c in df if 'N_studies' in c]
        df = df[df[n_study_cols].max(axis=1) > 1]
        df.to_csv(output, index=False)
        """
    stub:
        '''
        touch exwas_regions_meta_all_suggestive.csv
        '''
}

process make_singles_summary_table {
    publishDir "${launchDir}/Summary/"

    machineType 'n2-standard-4'

    input:
        path all_filtered_sumstats
    output:
        path('exwas_singles_meta_all_suggestive.csv')
    script:
        """
        #! ${params.my_python}

        import pandas as pd
        dfs = []
        input_list = '${all_filtered_sumstats.join(' ')}'.split()
        output = "exwas_singles_meta_all_suggestive.csv"

        # concatenate all results files together
        for f in input_list:
            dfs.append(pd.read_table(f))
        df = pd.concat(dfs).sort_values(by=['chr', 'pos'])

        # only present results in the top hits table with at least two studies
        n_study_cols = [c for c in df if 'N_studies' in c]
        df = df[df[n_study_cols].max(axis=1) > 1]
        df.to_csv(output, index=False)
        """
    stub:
        '''
        touch exwas_singles_meta_all_suggestive.csv
        '''
}

// Combine the N Samples Files to get a Table of Sample Size
process make_analysis_size_table {
    publishDir "${launchDir}/Summary/"

    machineType 'n2-standard-4'

    input:
        path(all_sample_sizes)
    output:
        path('exwas_meta_analysis_sizes.csv')
    script:
        """
        #! ${params.my_python}

        import pandas as pd
        dfs = []
        input_list = '${all_sample_sizes.join(' ')}'.split()
        output = "exwas_meta_analysis_sizes.csv"

        # concatenate sample size files together
        for f in input_list:
            dfs.append(pd.read_table(f, index_col='col_name', dtype=str))
        df = pd.concat(dfs, axis=1).transpose()
        df.columns = ['PHENO', 'ANALYSIS', 'N_Samples']
        df = df[['ANALYSIS', 'PHENO', 'N_Samples']]

        # make sure sample sizes are integers
        df['N_Samples'] = df['N_Samples'].str.replace('.0', '')
        df['N_Samples'] = df['N_Samples'].mask(df['N_Samples'].isin(['0', '-9']))
        print(df)

        df.sort_values(by=['ANALYSIS', 'PHENO']).to_csv(output, index=False, na_rep='NA')
        """
    stub:
        '''
        touch exwas_meta_analysis_sizes.csv
        '''
}

process collect_exwas_regions_plots {
    // Takes as input, list of manifest files and concatenates them
    publishDir "${launchDir}/Summary/"

    machineType 'n2-standard-16'

    label 'safe_to_skip'

    input:
        val(saige_analysis)
        path(plots_manifests)
    output:
        path("${saige_analysis}.plots_manifest.csv")
    script:
        """
        #! ${params.my_python}

        import pandas as pd
        dfs = []
        input_list = [x.strip() for x in "${plots_manifests}".replace('[', '').replace(']', '').split()]
        output_file = "${saige_analysis}.plots_manifest.csv"

        for file_path in input_list:
            dfs.append(pd.read_csv(file_path))

        pd.concat(dfs,ignore_index=True).to_csv(output_file, index=False)
        """

    stub:
        """
        touch ${saige_analysis}.plots_manifest.csv
        """
}

process collect_exwas_singles_plots {
    // Takes as input, list of manifest files and concatenates them
    publishDir "${launchDir}/Summary/"

    machineType 'n2-standard-16'

    label 'safe_to_skip'

    input:
        val(saige_analysis)
        path(plots_manifests)
    output:
        path("${saige_analysis}.plots_manifest.csv")
    script:
        """
        #! ${params.my_python}

        import pandas as pd
        dfs = []
        input_list = [x.strip() for x in "${plots_manifests}".replace('[', '').replace(']', '').split()]
        output_file = "${saige_analysis}.plots_manifest.csv"

        for file_path in input_list:
            dfs.append(pd.read_csv(file_path))

        pd.concat(dfs,ignore_index=True).to_csv(output_file, index=False)
        """

    stub:
        """
        touch ${saige_analysis}.plots_manifest.csv
        """
}
