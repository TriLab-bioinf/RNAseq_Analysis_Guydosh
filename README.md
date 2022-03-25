# RNAseq_analysis_RNAseL
R workflow for the analysis of bulk RNAseq data

## This RNAseq analysis workflow make use of two R scripts:
*RNAseL_analysis.Rmd:* This is the main R notebook for performing overall RNAseq data QC with pcaExplorer, differential gene expression analysis with DESeq2 and GO overrepresentation/enrichment analysis with CLusterProfiler.
*01_aux_rnaseq_functions.R:* This is an R scrip containing a number of auxiliary functions used by RNAseL_analysis.Rmd.
