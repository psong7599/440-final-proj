**OVERVIEW**


--------
**DATA**
Data represents gene expression matrices containing barcode and unique molecular identifier (UMI) counts generated from bTBI and control mice hippocampal samples (n = 6). To create bTBI group condition mice were subject to a 5.0 MPA blast shock wave created through rupturing an aluminum sheet with high pressure compressed gas [1]. 48 hours after injury (or not in control) the hippocampal region was isolated, treated, and purified for nuclei where single-cell RNA sequencing (scRNA-seq) could be performed with 10X Genomics Chromium platform [1]. Single cell transcripts were aligned to mouse reference genome (GRCm38 Ensembl: version 92) and provided feature-barcode matrices for each sample. Further the expression matrices are paired with a list of gene names that match the UMI for that cell to specific genes.

Publically available gene expression dataset uploaded to Gene Expression Omnibus (GEO) project number: GSE230253 (GSM7210822 for Sham and GSM7210821 for bTBI) (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230253).

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
