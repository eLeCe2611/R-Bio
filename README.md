# R-Bio: Bioinformatics Analysis Pipeline 🧬📊

This repository contains a collection of bioinformatics laboratory practices developed in **R**. The project covers the entire data analysis pipeline, from raw data processing (Microarrays and RNA-seq) to functional enrichment and biological network inference.

## 🔬 Overview

The practices are focused on analyzing high-throughput biological data, utilizing industry-standard libraries (Bioconductor) and datasets from public repositories like GEO (Gene Expression Omnibus).

### Key Areas Covered:
* **Data Preprocessing**: Normalization (RMA), quality control, and matrix processing.
* **Differential Expression**: Analysis of Microarray and RNA-seq datasets.
* **Functional Analysis**: Gene Ontology (GO) enrichment across three domains: Biological Process (BP), Cellular Component (CC), and Molecular Function (MF).
* **Systems Biology**: Construction of biological networks using Pearson correlation and network file generation (.sif).

## 📁 Repository Structure

The project is organized into practical units (P5 to P10):

| Unit | Focus | Key Components |
| :--- | :--- | :--- |
| **P5 & P6** | **Data Normalization** | RMA normalization, Microarray/RNA-seq processing, and count matrices. |
| **P7** | **Gene Ontology (GO)** | Functional enrichment analysis (BP, CC, MF) for multiple datasets. |
| **P8** | **Series Matrix Analysis** | Handling and analyzing GEO series matrix files. |
| **P9** | **Network Analysis** | Pearson correlation networks and Cytoscape-compatible (.sif) exports. |
| **P10** | **GEO SOFT Processing** | Advanced parsing and analysis of .soft metadata files. |

## 🛠️ Tech Stack
* **Language**: R
* **Main Ecosystem**: [Bioconductor](https://www.bioconductor.org/)
* **Key Libraries**: `limma`, `DESeq2`, `clusterProfiler`, `GSEABase`, and `ggplot2`.
* **Data Formats**: .txt, .csv, .soft, .sif (Simple Interaction File).

## 📝 Documentation
Each folder includes a **PDF report** (`EPD_LuisCarmona.pdf`) detailing the methodology used, the statistical results obtained, and the biological interpretation of the data.

## 👤 Author
* **Luis Carmona** - [eLeCe2611](https://github.com/eLeCe2611)

---
*This repository was created for academic purposes as part of the Bioinformatics curriculum.*
