# Bulk-RNAseq-Analysis

## Overview

This project presents an R-based workflow for bulk RNA-seq data analysis, integrating differential gene expression analysis with functional enrichment and pathway-level interpretation.
The pipeline is designed to investigate transcriptional changes between control and gene knockout conditions and explore their potential biological significance.

## Analytical Workflow

**1. Differential Gene Expression Analysis**

**2. Data Visualization**

**3. Functional Enrichment Analysis**
- Gene Ontology (GO): Biological Process, Molecular Function, and Cellular Component
- Gene Set Enrichment Analysis (GSEA) using ranked gene expression statistics
- KEGG pathway enrichment
- Reactome pathway enrichment
- Separate enrichment analyses of upregulated and downregulated genes

**4. Gene and Pathway Investigation**
- Identification of differentially expressed genes associated with selected biological pathways
- Gene-centered exploration of enriched pathways
- Investigation of candidate genes for further biological interpretation

## Tools and Libraries

R, DESeq2, clusterProfiler, ReactomePA, org.Hs.eg.db, enrichplot, ggplot2, pheatmap, EnhancedVolcano, DEGreport.

## Data

The workflow was developed for a six-sample control-versus-knockout RNA-seq design with three biological replicates per condition.

The analysis requires a raw gene count matrix. Research data are not included in this repository.

## Purpose

This repository demonstrates the integration of transcriptomic data analysis, statistical methods, and functional pathway interpretation to investigate molecular mechanisms underlying biological phenotypes.
