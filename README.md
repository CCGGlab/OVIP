This repository contains the scripts to recreate the figures which are based on RNA-seq data in Demuytere et al. _Effects of hyperthermia on cisplatin Tissue penetration and gene expression in peritoneal metastases: results from a randomized trial in Ovarian Cancer_, BJS 2024.

# Environment

Analysis was performed in a [Conda](https://anaconda.org/) environment. See [ovip.yml](ovip.yml) for details.

`conda env create -f ovip.yml -n ovip`

# Data

* Raw RNA-seq data is available upon reasonable request.
* Processed RNA-seq data (DESeq2 differential expression analysis, cf. methods) is available in `data/`:
    * `coldata.rds`: sample information (patient, treatment, dose, temperature)
    * `counts_matrix.rds`: count matrix (64252 features x 37 samples) created from aligned BAM files. Features named by ENSEMBL identifier (ENSG*)
    * `dds.rds`: `DESeqDataSet` created with _DESeq2_ (cf. methods)
    * `vst_dds.rds`: `DESeqTransform` object with variable stabilizing transformation (VST) applied on the original `DESeqDataSet`
    * `dge_results.rds`: list of `DESeqResults` for all group comparisons of dose (3) and temperature (1) plus combination of 120mg/41°C group vs. 75mg/37°C group
    * `genesets.rds`: genesets obtained from MSigDB v7.4 (except TFT_ls which is from RegNetwork):
        * `GO_BP_ls`: Gene Ontology Biological Process
        * `GO_CC_ls`: Gene Ontology Cellular Component
        * `GO_MF_ls`: Gene Ontology Molecular Function
        * `CP_ls`: Canonical pathways
        * `Ha_ls`: Hallmark
        * `Kegg_ls`: Kyoto Encyclopedia of Genes and Genomes 
        * `Rea_ls`: Reactome
        * `TFT_ls`: transcription factor genesets from the regulatory network repository (RegNetwork)
        * `CGP_ls`: chemical and genetic perturbations collection
        * `Immune_ls`: Immunologic Signature collection (C7)
        * `Cisplatin_ls`: Cisplatin response subset obtained by searching all collections of MSigDB
        * `Xenobiotic_ls`: Xenobiotic metabolism subset obtained by searching all collections of MSigDB
        * `OC_ls`: Ovarian cancer subset obtained by searching all collections of MSigDB
    * `boxplot_cibersort_x_lm22_relative.RData` (created via `scripts/icd_boxplot.R`)
        * `icd_sum`: aggregated immune cell deconvolution (ICD) results
        * `icd_pwt`: statistics for `icd_sum` (pairwise wilcoxon test)
        * `icd_boxplot`: figure 11B
    * `CIBERSORTx_Job4_Adjusted.txt`: CIBERSORTx results obtained from the [web tool](https://cibersortx.stanford.edu/)

# Figures

## Fig. 11

`source("scripts/manuscript_fig_11.R", encoding = "UTF-8")`

## Supp. Fig. 2

`source("/scripts/manuscript_supp_fig_2.R", encoding = "UTF-8")`
