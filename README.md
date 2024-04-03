**OVERVIEW**

Generally, this GitHub repository contains a script that takes the original single-cell RNA sequencing (scRNA-seq) gene expression data from blast-induced traumatic brain injury (bTBI) and sham (control) mice hippocampal samples and produces a Uniform Manifold Approximation and Projection for Dimension Reduction (UMAP) plot figure that visualizes cell clustering (using leidan clustering at 0.15 resolution) and highlights the bTBI vs control samples within each cluster. This plot is useful to visualize the clustering of the scRNA seq data into distinct cell types, where both bTBI and control samples are present in each cell type cluster. Further analysis can then be done to assign specific identities of cell types to each clusters based off marker genes and find differentially expressed genes between blast and control samples within each cell type.

The code first reads in the gene expression matrix taken from the paper (with associated barcodes for genes and corresponding UMI counts) and the associative list of gene names and converts both data sets (from both conditions) into one comprehensive data frame organized as gene by cell matrix (where entries represent UMI counts). Preprocessing/quality control is then conducted on the data frames where UMI counts are normalized, apoptotic/dead cells are removed, and the table is converted into an AnnData object such that ScanPy analysis can be performed (a Python library for single cell RNA sequencing analysis). Finally, analysis can then proceed where: highly variable genes are identified; principal component analysis allowed for reduction of the dimensionality of the dataset (where 50 dimensions was observed to explain most of the variance); cell neighborhoods are computed based on gene expression profiles; and UMAP creates low-dimensional representation of the data.  Visualized in the output figure is UMAP plots denoting: 1) clusters of cell gene expression profiles based on Leiden clustering (at 0.15 resolution) and 2) blast vs control cells.


--------
**DATA**

Data represents gene expression matrices containing barcode and unique molecular identifier (UMI) counts generated from bTBI and control mice hippocampal samples (n = 6). To create bTBI group condition mice were subject to a 5.0 MPA blast shock wave created through rupturing an aluminum sheet with high pressure compressed gas [1]. 48 hours after injury (or not in control) the hippocampal region was isolated, treated, and purified for nuclei where single-cell RNA sequencing (scRNA-seq) could be performed with 10X Genomics Chromium platform [1]. Single cell transcripts were aligned to mouse reference genome (GRCm38 Ensembl: version 92) and provided feature-barcode matrices for each sample. Further the expression matrices are paired with a list of gene names that match the UMI for that cell to specific genes.

Publically available gene expression dataset uploaded to Gene Expression Omnibus (GEO) project number: GSE230253 (GSM7210822 for Sham and GSM7210821 for bTBI) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230253).

[1] Zhang, Lingxuan, et al. “Single-nucleus transcriptomic mapping of blast-induced traumatic brain injury in mice hippocampus.” Scientific Data, vol. 10, no. 1, 20 Sept. 2023, https://doi.org/10.1038/s41597-023-02552-x. 

--------
**FOLDER STRUCTURE**

The repository is organized into three main folders: 

*Code*:
  
  containing the python script (hw6_20440_ps.py) that is utilized to download the dataset, preprocess the data, and create an UMAP visualization of the hippocampal cell clusters between control and blast TBI samples (that saves the figure into the “figures” folder.

*Data*:
  
  this folder contains the unprocessed gene expression matrix and gene list utilized in the python script. There are two   
  subfolders within this folder:
  
  (1) blast_ctrl
     
    this subfolder contains the unprocessed gene expression matrix for control hippocampal samples under "GSM7210822_genes_c-TBI.tsv" and associated gene name list under "GSM7210822_matrix_c-TBI.mtx".
 
  (2) blast_test
      
    this subfolder contains the unprocessed gene expression matrix for blast TBI hippocampal samples under "GSM7210821_genes_b-TBI.tsv" and associated gene name list under "GSM7210821_matrix_b-TBI.mtx".

*figures*:
 
  this folder is used to store the figure output of the provided python script: UMAP visualization of the hippocampal cell clusters between control and blast TBI samples saved under "umap_visuliazation_clusters.png"


--------
**INSTALLATION**

Run code by cloning repository and running the hw6_20440_ps.py script after downloading the required packages (you will need to have Python installed on your system along with necessary dependencies). Navigate into code folder: "cd Code" in terminal and then run "python hw6_20440_ps.py". Alternatively after cloning repository you can navigate to code folder and open preferred text editor (Atom, Sublime, etc.) then run script in editor.

Required packages:
  
  - anndata==0.10.6
  
  - numpy==1.25.2
  
  - pandas==2.2.1
  
  - scanpy==1.10.0
