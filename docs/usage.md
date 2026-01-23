# nf-core/eisca: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/eisca/usage](https://nf-co.re/eisca/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,merge,group
sampe_1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,,CONTROL
sampe_2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,,CONTROL
sampe_3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz,,CONTROL
sampe_4,AEG588A4_S4_L003_R1_001.fastq.gz,AEG588A4_S4_L003_R2_001.fastq.gz,,TREATMENT
sampe_5,AEG588A5_S5_L003_R1_001.fastq.gz,AEG588A5_S5_L003_R2_001.fastq.gz,,TREATMENT
sampe_6,AEG588A6_S6_L003_R1_001.fastq.gz,AEG588A6_S6_L003_R2_001.fastq.gz,sample_x,TREATMENT
sampe_7,AEG588A6_S6_L004_R1_001.fastq.gz,AEG588A6_S6_L004_R2_001.fastq.gz,sample_x,TREATMENT
```
If the input is smart-seq2 data, then you need to prepare a samplesheet as follows.
```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,plate,group
cell_1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,plate_1,CONTROL
cell_2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,plate_1,CONTROL
cell_3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz,plate_1,CONTROL
cell_4,AEG588A4_S4_L003_R1_001.fastq.gz,AEG588A4_S4_L003_R2_001.fastq.gz,plate_2,TREATMENT
cell_5,AEG588A5_S5_L003_R1_001.fastq.gz,AEG588A5_S5_L003_R2_001.fastq.gz,plate_2,TREATMENT
cell_6,AEG588A6_S6_L003_R1_001.fastq.gz,AEG588A6_S6_L003_R2_001.fastq.gz,plate_2,TREATMENT
cell_7,AEG588A6_S6_L004_R1_001.fastq.gz,AEG588A6_S6_L004_R2_001.fastq.gz,plate_2,TREATMENT
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `merge` | Optional: Users can specify closely related samples to be merged without the integration process by assigning the same new name to those merged samples. |
| `plate` | Optional: This is used for smart-seq2 data, indicating which Plate ID the cell belongs to. |
| `group` | Optional: Users can group a set of samples by assigning the same group name to those samples. Once the group column is added, all samples must be assigned to a group. |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.


## Parameters

The pipeline wutest has following parameters:

| Parameter   | Description |
| ----------- | ----------- |
| --input  \<samplesheet.csv> | Input samplesheet file in CSV format |
| --outdir \<directory> | Specify a output directory |
| --analyses | Specify analysis phases in the pipeline. There are three analysis phases (see detailed explanation below): primary, secondary, and tertiary. Multiple phases can be applied and should be separated by commas. By default all analysis phases are applied and set as "primary,secondary,tertiary".   |
| --skip \<string> | To skip one or multiple analyses, you can specify a set of analysis names separated by commas. e.g. "--skip fastqc" |
| --trim \<true/false> | A Boolean option, if set true the pipeline will apply tool Trimgalore to perform trimming to the input reads (false by default). |
| -profile \<config profile> | Specify a config profile to run the pipeline, which can be docker, singularity and conda |
| --protocol \<string> | The protocol that was used to generate the single cell data, e.g. 10x Genomics v2 Chemistry. Can be 'auto' (cellranger only), '10XV1', '10XV2', '10XV3', or any other protocol string that will get directly passed the respective aligner. Set to 'smartseq2' if the input is smart-seq2 data.|
| --aligner \<string> | Name of the tool to use for scRNA alignment, can choose from `alevin`, `kallisto` and `star`. |
| --h5ad \<file> | An AnnData object used when your analysis starting point is `secondary` or `tertiary`. If not specified then the one generated by specified aligner will be used. |
| --genome \<string> | Name of iGenomes reference, e.g. `--genome GRCh38` |
| --transcript_fasta \<file> | A cDNA FASTA file |
| --txp2gene \<string> | Path to transcript to gene mapping file. |
| --fasta \<file> | Reference FASTA gnome file |
| --gtf \<file> | Reference GTF annotation file |
| --parity \<string> | Indicating whether the input reads are paired-end (paired) or single-end (sinlge). It is paired by default |
| --args_qccellfilter \<string> | Flagged argument settings for the process of Cell filtering and QC, e.g. "--min_genes 50 --min_cells 1" |
| --args_clustering \<string> | Flagged argument settings for the process of clustering analysis, e.g. "--regress --scale" |
| --args_annotation \<string> | Flagged argument settings for the process of cell-type annotation analysis, e.g. "--model Immune_All_Low.pkl" |
| --args_annotation_scvi \<string> | Flagged argument settings for the process of cell-type annotation analysis with scvi-tools, e.g. "--h5ad_ref ${h5ad_ref} --label_key cellType --scvi_epochs 10 --batch_size 50 --scanvi_epochs 10" |
| --args_trainctmodel \<string> | Flagged argument settings for the process of training CellTypist models, e.g. "--feature_selection" |
| --args_dea \<string> | Flagged argument settings for the process of differential expression analysis, e.g. "--groupby leiden_res_0.50" |
| --args_dea_scvi \<string> | Flagged argument settings for the process of differential expression analysis with scvi-tools, e.g. "--groupby scanvi_label --epochs 10 --batch_size 50" |
| --save_reference \<true/false> | A Boolean option, if set true the pipeline will save all the intermediate output files apart from end results (true by default). |


## Analysis phases
The pipline has 3 analysis phases:
1. **primary phase** inculdes analyses:
   - Quality control of raw reads
   - Mapping reads and quantification
   - Convert count matrix into Anndata and Securat objects
2. **secondary phase** inculdes analyses:
   - Single-cell quality control
   - Cell filtering
   - Clustering analysis
   - Merging/integration of samples
3. **tertiary phase** inculdes analyses:    
   - Cell type annotation
   - Cell type annotation with scvi-tools
   - Differential expression analysis
   - Differential expression analysis with scvi-tools
   - Trajectory & pseudotime analysis (To be implemented)
   - Other downstream analyses (To be implemented)

Users can run each analysis phase separately by specifying the parameter, e.g., `--analyses secondary`. If this parameter is used with the `--skip`, the specified analyses within the analysis phases will be skipped. 

For example, the following command-line will run the pipeline for the secondary phase, skipping clustering analysis, and performing only cell filtering with the option min_genes set to 50.
```bash
nextflow run TGAC/eisca --analyses secondary --h5ad raw_matrix.h5ad --input samplesheet.csv --genome GRCh38 --outdir results -profile docker --aligner kallisto --skip clustering --protocol 10XV2 --args_qccellfilter "--min_genes 50" 
```

## Mapping & quantification
The primary analyses including mapping and quantification are based on pipeline `nf-core/scrnaseq', where we choose 3 aligners for our pipeline:
- [Salmon Alevin](https://salmon.readthedocs.io/en/latest/alevin.html) (default)
- [Kallisto](https://pachterlab.github.io/kallisto/about) & [Bustools](https://bustools.github.io/)
- [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md)

The corresponding parameters for Salmon Alevin are as follows:
| Parameter   | Description |
| ----------- | ----------- |
| --salmon_index  \<file> | This can be used to specify a precomputed Salmon index in the pipeline, in order to skip the generation of required indices by Salmon itself. |
| --txp2gene  \<file> | Path to transcript to gene mapping file. This allows the specification of a transcript to gene mapping file for Salmon Alevin and AlevinQC. |
| --simpleaf_rlen  \<int> | It is the target read length the index will be built for, using simpleaf. |

The corresponding parameters for Kqallisto are as follows:
| Parameter   | Description |
| ----------- | ----------- |
| --kallisto_index  \<file> | Specify a path to the precomputed Kallisto index. |
| --kb_t1c  \<file> | Specify a path to the cDNA transcripts-to-capture. |
| --kb_t2c  \<file> | Specify a path to the intron transcripts-to-capture. |
| --kb_workflow  \<string> | Type of workflow. Use nac for an index type that can quantify nascent and mature RNA. Use lamanno for RNA velocity based on La Manno et al. 2018 logic. (default: standard) |
| --kb_filter  \<file> | Activate Kallisto/BUStools filtering algorithm. |

The corresponding parameters for STAR are as follows:
| Parameter   | Description |
| ----------- | ----------- |
| --star_index  \<file> | Specify a path to the precomputed STAR index. |
| --star_ignore_asjbgtf  \<file> | Ignore the SJDB GTF file. |
| --seq_center  \<file> | Name of sequencing center for BAM read group tag. |
| --star_feature  \<file> | Quantification type of different transcriptomic feature. Use GeneFull on pre-mRNA count for single-nucleus RNA-seq reads. Use Gene Velocyto to generate RNA velocity matrix. |


## Cell filtering
Users can set the options for cell filtering in the parameter `--args_qccellfilter`, which are as follows. 
| Options   | Description |
| ----------- | ----------- |
| --min_genes  \<int> | Filter cells by minimum number of genes. (default=100)|
| --min_counts  \<int> | Filter cells by minimum number of counts. (default=1)|
| --max_genes  \<int> | Filter cells by maximum number of genes. (default=0 means not applied) |
| --max_counts  \<int> | Filter cells by maximum number of counts. (default=0 means not applied) |
| --min_cells  \<int> | Filter genes by number of cells expressed. (default=3)  |
| --min_gcounts  \<int> | Filter genes by the minimum counts expressed. (default=0) |
| --pct_mt  \<int> | Filter genes by the maximum percentage of mitochondrial counts. (default=20) |
| --quantile_upper  \<float> | Filter genes by upper limit of quantile on number of genes. (default=1) |
| --quantile_lower  \<float> | Filter genes by lower limit of quantile on number of genes. (default=0) |
| --iqr_coef  \<int> | Remove outliers which larger than iqr_coef*IQR in total_counts. (default=2) |
| --doublet_rate  \<int> | The expected fraction of transcriptomes that are doublets.  (default=0.1)|
| --find_doublets  \<int> | Whether to perform doublets prediction |
| --mt  \<string> | Specify a prefix of mitochondrial gene IDs. (default='MT-') |
| --fontsize  \<int> | Specify the font size for plots. (default=12) |
| --pdf | An switch of whether to generate figure files in PDF format. (false by default)|

For example, `--args_qccellfilter "--min_genes 50 --pct_mt 20"`


## Clustering analysis
Users can set the options for clustering analysis in the parameter `args_clustering`, which are as follows. 
| Options   | Description |
| ----------- | ----------- |
| --normalize | An switch of whether to apply normalization to the data (false by default)|
| --keep_doublets | An switch of whether to filter out the cells called as doublets. (false by default)|
| --regress | An switch of whether to regress out the variations from the total counts and the percentage of mitochondrial genes expressed. (false by default)|
| --scale | An switch of whether to scale the expression to have zero mean and unit variance. (false by default)|
| --resolutions \<string> | Resolution is used to control number of clusters. (default="0.02,0.05,0.1,0.5")|
| --integrate \<[bbknn, harmony, scanorama, scvi]> | Choose a method for data integration across samples. Currently four integration approaches can be choosen: 'bbknn' - a fast and intuitive batch effect removal method focus on local structure; 'harmony' - a popular global correction approach that iteratively adjusts the embedding of cells in lower-dimensional space, which is effective at correcting large batch effects, especially in datasets with complex batch structures. 'scanorama' - a method focuses on global merging with strong capability in retaining subtle variation across conditions. 'scvi' - uses a deep generative model (conditional variational autoencoder) to learn a latent representation of the data that accounts for batch variation while preserving biological variation, and good for highly heterogeneous datasets. scvi package requries CPUs support AVX instructions. (default=None)|
| --meta  \<[auto, sample, group]> | Choose a metadata column as the batch classes on which the clustering UMAPs will be displayed. By default, it is set to 'auto', which means it will use the 'group' column as the batch classes if 'group' is defined in the samplesheet file; otherwise, it will use the 'sample' column. |
| --covar_cat  \<[string]> | Specify the list of categorical covariates for scVI models. (default=[]) |
| --covar_con  \<[string]> | Specify the list of continuous covariates for scVI models. (default=[]) |
| --epochs  \<int> | Specify the number of epochs for scVI model training. (default=None) |
| --batch_size  \<int> | Specify the batch size for scVI model training. (default=None) |
| --fontsize  \<int> | Specify the font size for plots. (default=12) |
| --pdf | An switch of whether to generate figure files in PDF format. (false by default)|

For example, `--args_clustering "--integrate harmony"`


## Cell-type annotation analysis
Users can set the options for cell-type annotation analysis in the parameter `--args_annotation`, which are as follows. 
| Options   | Description |
| ----------- | ----------- |
| --model_file  \<string> | Specify a CellTypist model file for cell-type annotation analysis. |
| --model  \<string> | Specify a CellTypist model name, igored if a model file specified. (default='Immune_All_Low.pkl')|
| --mode  \<[best match, prob match]> | The way cell prediction is performed. 'best match' is to choose the cell type with the largest score/probability as the final prediction. Setting to 'prob match' will enable a multi-label classification, which assigns 0 (i.e., unassigned), 1, or >=2 cell type labels to each query cell. (default='best match')|
| --p_thres  \<float> | Probability threshold for the multi-label classification. Ignored if mode is 'best match'. (default=0.5) |
| --no_majority_voting | An switch of whether to disable the majority voting classifier after over-clustering. |
| --update_models | An switch of whether to update CellTypist models. |
| --meta  \<[auto, sample, group]> | Choose a metadata column as the batch classes on which the clustering UMAPs will be displayed. By default, it is set to 'auto', which means it will use the 'group' column as the batch classes if 'group' is defined in the samplesheet file; otherwise, it will use the 'sample' column. |
| --fontsize  \<int> | Specify the font size for plots. (default=12) |
| --pdf | An switch of whether to generate figure files in PDF format. (false by default)|

For example, `--args_annotation "--model_file path-to-file/my_model.pkl"`


## Training CellTypist models 
Users can set the options for training CellTypist models in the parameter `--args_trainctmodel`, which are as follows. 
| Options   | Description |
| ----------- | ----------- |
| --model_filename  \<string> | Specify a CellTypist model name. (default='celltypist_model.pkl')|
| --labels \<string> | Specify a column of cell-type from cell metadata of Anndata. (default='cell_type')|
| --l2c  \<float> | Inverse of L2 regularization strength for traditional logistic classifier. A smaller value can possibly improve model generalization while at the cost of decreased accuracy. This argument is ignored if SGD learning is enabled. (default=1.0) |
| --alpha  \<float> | L2 regularization strength for SGD logistic classifier. A larger value can possibly improve model generalization while at the cost of decreased accuracy. This argument is ignored if SGD learning is disabled. (default=0.0001) |
| --n_jobs  \<int> | Number of CPUs used, by default all CPUs are used. |
| --feature_selection | An switch of whether to perform two-pass data training where the first round is used for selecting important features/genes using SGD learning. If True, the training time will be longer. |
| --use_SGD | An switch of whether to implement SGD learning for the logistic classifier. |
| --use_GPU | An switch of whether to use GPU for logistic classifier. |

For example, `--args_trainctmodel "--model_filename test_model.pkl --labels majority_voting --feature_selection"`


## Cell-type annotation analysis with scvi-tools
Users can set the options for cell-type annotation analysis with scvi-tools in the parameter `--args_annotation_scvi`, which are as follows. 
| Options   | Description |
| ----------- | ----------- |
| --h5ad_ref  \<path> | Specify a reference anndata data file to generate scANVI model for annotation. |
| --model_path  \<path> | Specify a directory containing scANVI model file and corresponding reference anndata file. |
| --batch_key  \<string> | Specify a batch key for modelling reference data with scVI model. (default='sample')|
| --label_key  \<string> | Specify a label key of cell types for reference data. (default=None)|
| --n_top_genes  \<int> | Specify the number of highly-variable genes to keep. (default=2000)|
| --scvi_epochs  \<int> | Specify the number of epochs for training scVI model. (default=None)|
| --batch_size  \<int> | Specify the batch size for training scVI model. (default=None)|
| --scanvi_epochs  \<int> | Specify the number of epochs for training scANVI model. (default=None)|
| --n_samples_pl  \<int> | Specify the number of samples per label for training scANVI model. (default=None)|
| --meta  \<[auto, sample, group]> | Choose a metadata column as the batch classes on which the clustering UMAPs will be displayed. By default, it is set to 'auto', which means it will use the 'group' column as the batch classes if 'group' is defined in the samplesheet file; otherwise, it will use the 'sample' column. |
| --fontsize  \<int> | Specify the font size for plots. (default=12) |
| --pdf | An switch of whether to generate figure files in PDF format. (false by default)|

For example, `--args_annotation_scvi "--h5ad_ref ${h5ad_ref} --label_key cellType --scvi_epochs 10 --batch_size 50 --scanvi_epochs 10"`


## Differential expression analysis
Users can set the options for differential analysis in the parameter `--args_dea`, which are as follows. 
| Options   | Description |
| ----------- | ----------- |
| --groupby  \<string> | Specify a column of the observation table to define groups. (default='leiden') |
| --groups  \<string> | Specify a subset of groups, e.g. 'group1,group2'. By defualt, all groups are chosen. (default=None) |
| --reference  \<string> | Users can spcecify a group name as reference, and all other groups will be comapred against with this group. By default each group will be compared against rest of groups. (default='rest') |
| --method  \<['t-test', 'wilcoxon', 'logreg', 't-test_overestim_var']> | Choose a test method for differential expression anlaysis. The default method is 't-test', 't-test_overestim_var' overestimates variance of each group, 'wilcoxon' uses Wilcoxon rank-sum, 'logreg' uses logistic regression. (default='t-test')|
| --n_genes  \<int> | Number of top marker genes to show in plots. (default=20) |
| --celltype_col \<string> | Spcecify a column of the observation table to define cell-types, and DEA will be performed between groups for each cell-type. (default=None) |
| --celltypes \<string> | Spcecify a subset of cell-types for DEA between groups, e.g. 'celltype1,celltype2'. By default all cell-types are used. (default=None) |
| --meta  \<[auto, sample, group]> | Choose a metadata column as the batch classes on which the clustering UMAPs will be displayed. By default, it is set to 'auto', which means it will use the 'group' column as the batch classes if 'group' is defined in the samplesheet file; otherwise, it will use the 'sample' column. |
| --fontsize  \<int> | Specify the font size for plots. (default=12) |
| --pdf | An switch of whether to generate figure files in PDF format. (false by default)|

For example:  
`--args_dea "--groupby leiden_res_0.50"` - perform DEA to find marker genes for each cluster against the rest using clusters defined in column 'leiden_res_0.50' at group level if 'group' is defined in the metadata. Applying `--meta sample` to perform DEA at sample level.  
`--args_dea "--groupby group --reference control"` - perform DEA to find DE genes between each group against the group 'control', groups are defined in column 'group'.  
`--args_dea "--groupby group --reference control --celltype_col majority_voting"` - same as above but for each cell-type defined in column 'majority_voting'.


## Differential expression analysis with scvi-tools
Users can set the options for differential analysis with scvi-tools in the parameter `--args_dea_scvi`, which are as follows. 
| Options   | Description |
| ----------- | ----------- |
| --groupby  \<string> | Specify a column of the observation table to define groups. (default='None') |
| --group1  \<string> | Specify a subset of groups, e.g. 'group1,group2'. By defualt, all groups are chosen. (default=None) |
| --group2  \<string> | Specify a subset of groups, e.g. 'group1,group2'. By defualt, all groups are chosen. (default=None) |
| --covar_cat  \<[string]> | Specify the list of categorical covariates for scVI models. (default=[]) |
| --covar_con  \<[string]> | Specify the list of continuous covariates for scVI models. (default=[]) |
| --epochs  \<int> | Specify the number of epochs for scVI model training. (default=None) |
| --batch_size  \<int> | Specify the batch size for scVI model training. (default=None) |
| --celltype_col \<string> | Spcecify a column of the observation table to define cell-types, and DEA will be performed between groups for each cell-type. (default=None)|
| --celltypes \<string> | Spcecify a subset of cell-types for DEA between groups, e.g. 'celltype1,celltype2'. By default all cell-types are used. (default=None) |
| --n_markers  \<int> | Number of top marker genes to show in plots. (default=3) |
| --deg_lfc  \<int> | Set threshold of Log Folder Change for DEGs. (default=0) |
| --deg_bayes  \<int> | Set threshold of bayes factor for DEGs. (default=0) |
| --deg_nzerosprop  \<int> | Set threshold of proportion of non-zero expression cells in group1 for DEGs. (default=0) |
| --meta  \<[auto, sample, group]> | Choose a metadata column as the batch classes on which the clustering UMAPs will be displayed. By default, it is set to 'auto', which means it will use the 'group' column as the batch classes if 'group' is defined in the samplesheet file; otherwise, it will use the 'sample' column. |
| --fontsize  \<int> | Specify the font size for plots. (default=12) |
| --pdf | An switch of whether to generate figure files in PDF format. (false by default)|

For example:  
`--args_dea_scvi "--groupby scanvi_label --epochs 10 --batch_size 50"` - perform DEA to find marker genes for cells of each cell type defined in metadata column `scanvi_label` against the rest cells at group level if `group` is defined in the metadata. Applying `--meta sample` to perform DEA at sample level.  
`--args_dea_scvi "--groupby group --group1 Control --group2 Treatment --celltype_col scanvi_label --celltypes naive_B_cell,classical_monocyte --epochs 10 --batch_size 50"` - perform DEA to find DE genes by comparing group `Control` against the group 'Treatment' for 2 subpopulations of cell types defined in metadata column `scanvi_label`. Comparing groups are defined in metadata column `group`.  
`--args_dea_scvi "--groupby group --group1 Control --group2 Treatment --epochs 10 --batch_size 50"` - perform DEA to find DE genes by comparing group `Control` against the group 'Treatment' for for all cells.


## Cell-cell communication analysis
Users can set the options for cell-cell communication analysis in the parameter `--args_cellchat`, which are as follows. 
| Options   | Description |
| ----------- | ----------- |
| --group  \<string> | Specify the column from metadata for grouping the cells. (default=majority_voting)|
| --n_actgrps  \<int> | Specify the number of top active groups based on counts of the ligand–receptor interactions. (default=12)|
| --n_pairLR  \<int> | Specify the number of top ligand–receptor (LR) pairs based on the probabilities of LR pairs. (default=100)|
| --normalize | Indicates whether to normalize the counts.|
| --db  \<[human, mouse]> | Specify the species of CellChatDB. (default=human) |
| --dbc  \<string> | Specify the categories of CellChatDB, e.g. Secreted Signaling. |
| --dbenps | Indicates whether to use all CellChatDB excepting 'Non-protein Signaling'  |
| --threads  \<int> | Number of threads for parallel runs. (default=4) |
| --mean_method  \<[trimean, truncatedMean]> | Specify the method for calculating the average gene expression per cell group. Method trimean produces fewer interactions than method truncatedMean. (default=trimean) |
| --mincells  \<int> | The minimum number of cells required in each cell group. (default=10) |
| --meta  \<[auto, sample, group]> | Choose a metadata column to define separate subsets of cells for analysis. (default=auto)|
| --pdf | An switch of whether to generate figure files in PDF format. (false by default)|

For example, `--args_cellchat "--group cellType --mean_method truncatedMean"`


## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run TGAC/eisca --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/eisca -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/eisca
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/eisca releases page](https://github.com/nf-core/eisca/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
