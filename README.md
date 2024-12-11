
Documentation for ExWAS Meta-Analysis
=====================================

# Module Overview


This module performs many different types of statistical tests for meta-analyzing effects and/or p-values from exome-wide association studies (ExWASs). It looks at both gene-burden region-based tests and single rare variants, performing analyses including Fisher’s tests and Inverse Variance-Weighted meta-analysis.
## Cloning Github Repository:


* Command: `git clone https://github.com/PMBB-Informatics-and-Genomics/geno_pheno_workbench.git`

* Navigate to relevant workflow directory run commands (our pipelines assume all of the nextflow files/scripts are in the current working directory)
## Software Requirements:


* [Nextflow version 23.04.1.5866](https://www.nextflow.io/docs/latest/cli.html)

* [Singularity version 3.8.3](https://sylabs.io/docs/) OR [Docker version 4.30.0](https://docs.docker.com/)
## Commands for Running the Workflow


* Singularity Command: `singularity build exwas_meta.sif docker://guarelin/exwas_meta:latest`

* Docker Command: `docker pull docker://guarelin/exwas_meta:latest`

* Run Command: `nextflow run exwas_meta_analysis.nf -profile cluster`

* Common `nextflow run` flags:

    * `-resume` flag picks up the workflow where it left off, otherwise, the workflow will rerun from the beginning

    * `-stub` performs a sort of dry run of the whole workflow, checks channels without executing any code

    * `-profile` selects the compute profiles we set up in nextflow.config (see nextflow.config file below)

    * `-profile` selects the compute profiles we set up in nextflow.config (see nextflow.config file below)

    * `-profile standard` uses the docker image to executes the processes

    * `-profile cluster` uses the singularity container and submits processes to a queue- optimal for HPC or LPC computing systems

    * `-profile all_of_us` uses the docker image to execute pipelines on the All of Us Researcher Workbench

* for more information visit the [Nextflow documentation](https://www.nextflow.io/docs/latest/cli.html)
# Configuration Parameters and Input File Descriptions

## Workflow


* `quant_pheno_list` (Type: List)

    * Quantitative phenotype list

* `analyses` (Type: Map (Dictionary))

    * Map of lists where keys are meta-analysis group nicknames and lists are groups of cohorts to include in that meta-analysis. This allows for multiple combinations of meta-analyses, for example all cohorts of one sex/ancestry, leave-one-biobank-out.

* `bin_pheno_list` (Type: List)

    * Binary phenotype list

* `my_python` (Type: File Path)

    * Path to the python executable to be used for python scripts - often it comes from the docker container (/opt/conda/bin/python)

    * Corresponding Input File: Python Executable

        * This is the /path/to/the/correct/python executable that has the right packages installed and right version. It usually comes from the docker/singularity container. If not, be sure to have most of the standard data packages like pandas, numpy, scipy, matplotlib, seaborn, etc

        * Type: Executable

        * Format: exe
## Pre-Processing


* `regions_sumstats_suffix` (Type: String)

    * Suffix for regions files from the ExWAS summary stats

    * Corresponding Input File: ExWAS Regions Summary Statistics

        * Input summary statistics of the regions tests. They are expected to be organized in the directory from which you’re running the pipeline like so: “COHORT/Sumstats/PHENO.SUFFIX”. This matches output of your other workflows, but if you’re starting here, you may use the python script “scripts/set_up_cohort_directory_structure.py” to create symlinks with the correct structure

        * Type: Summary Statistics

        * Format: tsv.gz

        * Input File Header:





        ```
        BETA    OR      SE      P       N       N_CASES N_CONTROLS      REGION  MAX_MAF ANNOT
        0.10749207338367353     1.113482034566681       0.10861255486849287     0.48893376459755766     16082   5316    10766   ENSG00000160256 0.01    pLOF
        0.2752070632715491      1.3168033082414594      1.1376977756875413      0.7748787400623763      16082   5316    10766   ENSG00000160256 0.01    damaging_missense
        0.09866625123248854     1.103697880275314       0.04777625246679377     0.09458445160078595     16082   5316    10766   ENSG00000160256 0.01    other_missense
        0.1869359835432516      1.2055501075333186      0.05549816240001925     0.0027432656092397415   16082   5316    10766   ENSG00000160256 0.01    synonymous
        
        ```

* `singles_sumstats_suffix` (Type: String)

    * Suffix for singles files from the ExWAS summary stats

    * Corresponding Input File: ExWAS Singles Summary Statistics

        * Input summary statistics of the singles tests. They are expected to be organized in the directory from which you’re running the pipeline like so: “COHORT/Sumstats/PHENO.SUFFIX”. This matches output of your other workflows, but if you’re starting here, you may use the python script “scripts/set_up_cohort_directory_structure.py” to create symlinks with the correct structure

        * Type: Summary Statistics

        * Format: tsv.gz

        * Input File Header:





        ```
        #CHROM  BP      A1      A2      BETA    OR      SE      P       A1_FREQ N       N_CASES N_CONTROLS
        21      46482292        G       T       0.06082811129594274     1.0627162295714525      0.061462174537242044    0.48893376459755766     0.10478298166200951     49575   12384      37191
        21      38975062        T       C       0.17065668690623242     1.1860834811261391      0.705489644020694       0.7748787400623763      0.026454221287319335    49575   12384      37191
        21      44118844        A       C       0.1934617601789459      1.2134429838217287      0.09367821095381963     0.09458445160078595     0.2519137577341014      49575   12384      37191
        21      24024049        G       A       0.2981114423114272      1.3473119270688356      0.08850429395727123     0.0027432656092397415   0.07047045978001559     49575   12384      37191
        
        ```

* `regions_info_cols` (Type: Map (Dictionary))

    * A map with six keys:  region, annot_group, max_maf, n, n_case, n_control, where the values are the corresponding test information columns in the regions files.

* `singles_info_cols` (Type: Map (Dictionary))

    * A map with three keys: beta, se, and p_value, where the values are the corresponding test statistic columns in the singles files.

* `burden_cols` (Type: Map (Dictionary))

    * A map with three keys: beta, se, and p_value, where the values are the corresponding gene burden test statistic columns in the regions files.

* `skat_p_col` (Type: String)

    * The name of the SKAT test p-value from the ExWAS summary stats. It can be set to null if you don’t want to meta-analyze these p-values.

* `skato_p_col` (Type: String)

    * The name of the SKAT-O test p-value from the ExWAS summary stats. It can be set to null if you don’t want to meta-analyze these p-values.

* `singles_effect_cols` (Type: Map (Dictionary))

    * A map with seven keys: chr, pos, effect_allele, other_allele, n, n_case, and n_control, where the values are the corresponding test information columns in the singles files.
## Post-Processing


* `gene_location_file ` (Type: File Path)

    * This file is used for getting gene-based coordinates for plotting .

    * Corresponding Input File: Gene Location File

        * CSV file of 

        * Type: Data Table

        * Format: tsv

        * Input File Header:





        ```
        gene_id chromosome  seq_region_start    seq_region_end  gene_symbol
        GENE1   1   1   90  GS1
        GENE2   2   91  100 GS2
        ```

* `p_cutoff_summarize` (Type: Float)

    * P-Value Threshold for Summarizing Results at the End, arbitrary p-value threshold for creating a table of results combined with low p-values 
# Output Files from ExWAS_Meta-Analysis


* Regions Meta-Analysis Manhattan Plots

    * A dot plot (manhattan plot) of significant gene regions associated with a phenotype. One plot will be created for each unique combination of phenotype, cohort, annotation group (pLof, etc.), and MAF threshold. 

    * Type: Manhattan Plot

    * Format: png

    * Parallel By: Cohort, Phenotype, Annot Group, MAF

* Regions Meta-Analysis QQ Plots

    * A QQ Plot of the Null Model vs Log10P results of the analysis for gene regions.  One plot will be created for each unique combination of phenotype, cohort, annotation group (pLof, etc.), and MAF threshold. 

    * Type: QQ Plot

    * Format: png

    * Parallel By: Cohort, Phenotype, Annot Group, MAF

* Regions Meta-Analysis Summary Statistics

    * A gzipped, unfiltered TSV (tab-separated) file of the results for the gene (regions) analysis if run. One file will be created for each unique Cohort, Phenotype, and analysis (regular, cauchy, rare, ultra rare) combination.

    * Type: Summary Statistics

    * Format: tsv.gz

    * Parallel By: Cohort, Phenotype

    * Output File Header:





    ```
    phenotype|gene           |annot            |max_maf|p_value           |p_value_burden    |p_value_skat      |beta_burden        |se_burden         |mac   |mac_case|mac_control|rare_var_count|ultrarare_var_count
    T2Diab   |ENSG00000141956|pLoF             |0.0001 |0.0479451461682565|0.0479451461682565|0.0479451461682565|0.0588652858997042 |0.0297621953331829|12.0  |5.0     |7.0        |0.0           |9.0
    T2Diab   |ENSG00000141956|pLoF             |0.001  |0.0479451461682565|0.0479451461682565|0.0479451461682565|0.0588652858997042 |0.0297621953331829|12.0  |5.0     |7.0        |0.0           |9.0
    T2Diab   |ENSG00000141956|pLoF             |0.01   |0.0479451461682565|0.0479451461682565|0.0479451461682565|0.0588652858997042 |0.0297621953331829|12.0  |5.0     |7.0        |0.0           |9.0
    T2Diab   |ENSG00000141956|damaging_missense|0.0001 |0.464219450219203 |0.464219450219203 |0.464219450219203 |-0.0110759683619445|0.0151328276810456|52.0  |7.0     |45.0       |0.0           |41.0
    ```

* Regions Meta-Analysis Top Hits Table

    * A FILTERED top hits csv summary file of results including cohort, phenotype, gene, group annotation, p-values, and other counts. One single summary file will be aggregated from all the “top hits” in each “Regions Summary Statistics” file. 

    * Type: Summary Table

    * Format: csv

    * Output File Header:





    ```
    region,annot_group,max_maf,chr,pos_start,pos_stop,gene_symbol,analysis,phenotype,p_burden_stouffer_meta,p_burden_stouffer_N_eff,p_burden_stouffer_N_studies,beta_burden_inv_var_meta,se_burden_inv_var_meta,p_burden_inv_var_meta,N_eff_inv_var_meta,N_studies_inv_var_meta,p_burden_chi2_stat,p_burden_heterogeneity
    ENSG00000000419,damaging_missense,0.0001,20,50934867,50959140,DPM1,Leave_EUR_Out,T2D,0.606368639793734,11702.0,3,0.0336482650736499,0.0653029773095249,0.6063686397937349,11702.0,3,,
    ENSG00000000419,damaging_missense,0.0001,20,50934867,50959140,DPM1,AFR_EUR,LDL_median,0.4348806366999341,24879.0,2,-0.2393382100454541,0.3751452681678022,0.5234814383226609,24879.0,2,38.57666608752835,
    ENSG00000000419,damaging_missense,0.0001,20,50934867,50959140,DPM1,ALL_F,LDL_median,0.0654224944953294,12515.0,2,-0.9459153671090688,0.5391298446427648,0.0793410420455428,12515.0,2,25.260017289413383,
    ENSG00000000419,damaging_missense,0.0001,20,50934867,50959140,DPM1,ALL_M,AAA,0.4794196005233119,15172.0,2,-0.0451279677642891,0.0638088900735573,0.4794196005233122,15172.0,2,,
    
    ```

* Singles Meta-Analysis Manhattan Plots

    * A dot plot (manhattan plot) of significant variants associated with a phenotype. One plot will be created for each unique combination of phenotype, cohort, annotation group (pLof, etc.), and MAF threshold. 

    * Type: Manhattan Plot

    * Format: png

    * Parallel By: Cohort, Phenotype

* Singles Meta-Analysis QQ Plots

    * A QQ Plot of the Null Model vs Log10P results of the analysis for variants.  One plot will be created for each unique combination of phenotype, cohort, annotation group (pLof, etc.), and MAF threshold. 

    * Type: QQ Plot

    * Format: png

    * Parallel By: Cohort, Phenotype

* Singles Meta-Analysis Summary Statistics

    * A gzipped, unfiltered TSV (tab-separated) file of the results for the variant (singles) analysis. One file will be created for each unique Cohort, Phenotype, and analysis (regular, cauchy, rare, ultra rare) combination.

    * Type: Summary Statistics

    * Format: tsv.gz

    * Parallel By: Cohort, Phenotype

    * Output File Header:





    ```
    phenotype|chromosome|base_pair_location|variant_id        |other_allele|effect_allele|effect_allele_count|effect_allele_frequency|missing_rate|beta      |standard_error|t_statistic|variance|p_value   |p_value_na|is_spa_test|allele_freq_case|allele_freq
    T2Diab   |21        |41801254          |21_41801254_TCTG_T|TCTG        |T            |277                |0.0046611              |0.0         |-0.099231 |0.167775      |-3.52526   |35.5258 |0.5542179 |0.5542179 |False      |0.00426841      |0.00474126
    T2Diab   |21        |41801360          |21_41801360_C_T   |C           |T            |41                 |0.00068991             |0.0         |-0.864441 |0.633121      |-3.98228   |5.38237 |0.08606924|0.08606924|False      |0.000297796     |0.000769948
    T2Diab   |21        |41801603          |21_41801603_C_T   |C           |T            |24                 |0.00040385             |0.0         |0.322923  |0.570593      |0.991852   |3.07148 |0.5714322 |0.5714322 |False      |0.000496327     |0.000384974
    T2Diab   |21        |41801645          |21_41801645_G_A   |G           |A            |58                 |0.000975971            |0.0         |0.0167811 |0.35132       |0.135962   |8.10206 |0.9619027 |0.9619027 |False      |0.00109192      |0.000952304
    ```

* Singles Meta-Analysis Top Hits Table

    * A FILTERED top hits csv summary file of results including cohort, phenotype, gene, group annotation, p-values, and other counts. One single summary file will be aggregated from all the “top hits” in each “Singles (Variant) Summary Statistics” file. 

    * Type: Summary Table

    * Format: csv

    * Output File Header:





    ```
    chr,pos,effect_allele,other_allele,analysis,phenotype,p_single_stouffer_meta,p_single_stouffer_N_eff,p_single_stouffer_N_studies,beta_single_inv_var_meta,se_single_inv_var_meta,p_single_inv_var_meta,N_eff_inv_var_meta,N_studies_inv_var_meta,p_single_chi2_stat,p_single_heterogeneity
    1,69745,T,C,ALL_M,AAA,0.579096,15172.0,2,-1.06221,1.91491,0.5790965097927716,15172.0,2,,
    1,930282,A,G,ALL_M,LDL_median,0.4106934,12364.0,2,-8.493,10.3236,0.410691052172855,12364.0,2,,
    1,935839,T,C,ALL,T2D,0.619276,39632.0,2,-0.334515,0.673235,0.6192757781492377,39632.0,2,,
    1,935849,C,G,ALL,T2D,0.1225267999999999,39632.0,2,-0.515683,0.333937,0.1225272089986704,39632.0,2,,
    ```

* Meta-Analysis Sample Sizes

    * A table containing the maximum sample size of each of the meta-analyses based on the input cohorts and phenotypes. The actual numbers for each test may vary if there is missingness for certain variants, but this captures the largest sample size.

    * Type: Summary Table

    * Format: csv

    * Output File Header:





    ```
    ANALYSIS,PHENO,N_Samples
    AFR_EUR,AAA,31265
    AFR_EUR,AAA,31265
    AFR_EUR,BMI_median,38134
    AFR_EUR,BMI_median,38134
    
    ```
# Example Config File Contents


```
params {
    // Map the overall "analyses" (meta-analysis combinations) to cohort/study population lists
    analyses = [
        'AFR_EUR': ['PMBB_AFR_ALL', 'PMBB_EUR_ALL'],
        'ALL': ['PMBB_AFR_ALL', 'PMBB_EUR_ALL', 'PMBB_EAS_ALL', 'PMBB_AMR_ALL', 'PMBB_SAS_ALL'],
        'ALL_M': ['PMBB_AFR_M', 'PMBB_EUR_M', 'PMBB_EAS_M', 'PMBB_AMR_M', 'PMBB_SAS_M'],
        'ALL_F': ['PMBB_AFR_F', 'PMBB_EUR_F', 'PMBB_EAS_F', 'PMBB_AMR_F', 'PMBB_SAS_F'],
        'Leave_EUR_Out': ['PMBB_AFR_ALL', 'PMBB_EAS_ALL', 'PMBB_AMR_ALL', 'PMBB_SAS_ALL']
    ]

    // Executables for python and GWAMA
    // my_python = '/home/guarelin/mambaforge/envs/py311/bin/python'
    my_python = '/opt/conda/bin/python'
    // my_Rscript = '/opt/conda/bin/Rscript'

    // Lists of phenotypes
    bin_pheno_list =  ['T2D', 'AAA']
    quant_pheno_list = ['LDL_median', 'BMI_median']

    // Pre- and Post-Processing Params (probably starts with .)
    regions_sumstats_suffix = '.exwas_regions.saige.gz'
    singles_sumstats_suffix = '.exwas_singles.saige.gz'

    // Top-Hits tables will be filtered to this p-value
    p_cutoff_summarize = 0.00001
    // this is for getting ENSEMBL gene symbols and coordinates for summary stats and plotting
    // tab-separated, columns include: gene_id, chromosome, seq_region_start, seq_region_end, gene_symbol
    gene_location_file = '/path/to/data/homo_sapiens_111_b38.txt'

    // Meta-Analysis Test Info
    regions_info_cols = [
        'region': 'gene',
        'annot_group': 'annot',
        'max_maf': 'max_maf',
        'n': 'N',
        'n_case': 'N_case',
        'n_control': 'N_ctrl'
    ]

    // Single-Variant Test Info
    singles_info_cols = [
        'chr': 'chromosome',
        'pos': 'base_pair_location',
        'effect_allele': 'effect_allele',
        'other_allele': 'other_allele',
        'n': 'n',
        'n_case': 'n_case',
        'n_control': 'n_ctrl'
    ]

    // set any of these column parameters to null 
    // if you don't want to meta-analyze those effects
    burden_cols = ['beta': 'beta_burden', 'se': 'se_burden', 'p_value' : 'p_value_burden']
    skat_p_col = 'p_value_skat'
    skato_p_col = null
    singles_effect_cols = ['beta': 'beta', 'se': 'standard_error', 'p_value': 'p_value']

    // Options for which tests to run
    tests_to_run = [
        'fishers' : true,
        'one_tailed_weighted_fishers': true,
        'two_tailed_weighted_fishers': false,
        'cauchy' : true,
        'one_tailed_stouffer': true,
        'two_tailed_stouffer' : true,
        'heterogeneity' : true,
        'inverse_variance_weighted' : true
    ]

    // Annotate single rare variants with nearest gene (should be the gene it's in)
    annotate = true

    // The following arguments go with annotate=true and will be used by the biofilter_wrapper sub-workflow
    biofilter_build = '38' // can be 19 or 38
    biofilter_loki = '/path/to/data/loki.db'
    biofilter_script = '/app/biofilter.py' // Must be an executable python file
    biofilter_close_dist = 5E4 // How "close" to a gene a variant should be to be considered "close" vs "far"
}
```
# Current Dockerfile for the Container/Image


```docker
FROM continuumio/miniconda3
WORKDIR /app

# biofilter version argument
ARG BIOFILTER_VERSION=2.4.3

RUN apt-get update \    
    # install packages needed to install biofilter and NEAT-plots
    && apt-get install -y --no-install-recommends libz-dev g++ gcc git wget tar unzip make \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    # install python packages needed for pipeline
    && conda install -y -n base -c conda-forge wget libtiff conda-build scipy pandas seaborn matplotlib numpy apsw sqlite \
    && conda clean --all --yes \
    # install NEAT-plots
    && git clone https://github.com/PMBB-Informatics-and-Genomics/NEAT-Plots.git \
    && mv NEAT-Plots/manhattan-plot/ /app/ \
    && conda develop /app/manhattan-plot/ \
    # install biofilter
    && wget https://github.com/RitchieLab/biofilter/releases/download/Biofilter-${BIOFILTER_VERSION}/biofilter-${BIOFILTER_VERSION}.tar.gz -O biofilter.tar.gz \
    && tar -zxvf biofilter.tar.gz --strip-components=1 -C /app \
    && /opt/conda/bin/python setup.py install \
    # make biofilter executable
    && chmod a+rx /app/biofilter.py \
    # remove biofilter tarball and NEAT-plots directory
    && rm -R biofilter.tar.gz NEAT-Plots

USER root
```
# Current `nextflow.config` contents


```
includeConfig 'exwas_meta_analysis.config'

profiles {
    non_docker_dev {
        process.executor = awsbatch-or-lsf-or-slurm-etc
        process.queue = 'epistasis_normal'
        process.memory = '15GB'
    }

    standard {
        process.executor = awsbatch-or-lsf-or-slurm-etc
        process.container = 'guarelin/exwas_meta:latest'
        docker.enabled = true
    }

    cluster {
        process.executor = awsbatch-or-lsf-or-slurm-etc
        process.queue = 'epistasis_normal'
        process.memory = '15GB'
    	process.container = 'exwas_meta.sif'
        singularity.enabled = true
        singularity.runOptions = '-B /root/,/directory/,/names/'
    }

    all_of_us {
        process.executor = awsbatch-or-lsf-or-slurm-etc
        process.memory = '15GB'
        process.container = 'gcr.io/ritchie-aou-psom-9015/exwas_meta:latest'
        docker.enabled = true
    }
}

```