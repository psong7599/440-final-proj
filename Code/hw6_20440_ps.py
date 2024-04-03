import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad


# Data Formatting: Read in blast and control mice scRNA seq hippocampal data and reformat into matrix where rows are
# cells and columns of genes (with associated UMI counts)

# BLAST SAMPLES
# read in the data from file and format into correct gene expression matrix shape

# Path to the .mtx file that contains the blast UMI data
blast_umi_path = '..\\Data\\blast_test\\GSM7210821_matrix_b-TBI.mtx'

# Read in the data and store in dataframe
blast_df = pd.read_csv(blast_umi_path, skiprows=3, names=["Col 1"])

# Split the values in the "Col 1" column into three separate columns
blast_df[['Gene', 'Cell #', 'UMI']] = blast_df['Col 1'].str.split(' ', expand=True)

# Drop the original "Col 1" column
blast_df.drop(columns=['Col 1'], inplace=True)

# Path to the file that contains the gene names associated with the UMI counts
blast_gene_path = '..\\Data\\blast_test\\GSM7210821_genes_b-TBI.tsv'

# Read in the data and store in dataframe
blast_gene_df = pd.read_csv(blast_gene_path, names=["Gene #", 'Gene Name', 'etc'], sep='\t')

# Assign values to each gene name for matching in the UMI matrix
blast_gene_df['Count'] = range(1, len(blast_gene_df) + 1)

# Convert 'Count' column to object type
blast_gene_df['Count'] = blast_gene_df['Count'].astype(str)

# Merge the two DataFrames based on the 'Gene' column from the first DataFrame and the 'Count' column from the second
# DataFrame
merged_df = blast_df.merge(blast_gene_df[['Gene Name', 'Count']], left_on='Gene', right_on='Count', how='left')

# Create a DataFrame filled with zeros that will be filled in with the desired UMI values
zeros_df = pd.DataFrame(np.zeros((9823, len(blast_gene_df['Gene Name']))))

# Define column labels
columns = blast_gene_df['Gene Name']
# Define row labels
index = range(1, 9824)

# Set column and row labels
zeros_df.columns = columns
zeros_df.index = index
# Add a label to the index
zeros_df = zeros_df.rename_axis('Cell #')

# Convert 'Cell #' column to int type
merged_df['Cell #'] = merged_df['Cell #'].astype(int)

# Pivot the second matrix DataFrame to get the desired format for the gene expression matrix
pivoted_df = merged_df.pivot_table(index='Cell #', columns='Gene Name', values='UMI', fill_value=0, aggfunc='sum')

# Reindex the pivoted DataFrame to match the index of the first matrix DataFrame (if needed)
pivoted_df = pivoted_df.reindex_like(zeros_df)

# Fill missing values with 0 (if there are cells in the first matrix DataFrame not present in the second)
pivoted_df = pivoted_df.fillna(0)

# Print the pivoted DataFrame
print(pivoted_df)



# CONTROL SAMPLES
# read in the data from file and format into correct gene expression matrix shape

# Path to the .mtx file that contains the control UMI data
b_ctrl_umi_path = '..\\Data\\blast_ctrl\\GSM7210822_matrix_c-TBI.mtx'

# Read in the data and store in dataframe
b_ctrl_df = pd.read_csv(b_ctrl_umi_path, skiprows=3, names=["Col 1"])

# Split the values in the "Col 1" column into three separate columns
b_ctrl_df[['Gene', 'Cell #', 'UMI']] = b_ctrl_df['Col 1'].str.split(' ', expand=True)

# Drop the original "Col 1" column
b_ctrl_df.drop(columns=['Col 1'], inplace=True)

# Path to the file that contains the gene names associated with the UMI counts
b_ctrl_gene_path = '..\\Data\\blast_ctrl\\GSM7210822_genes_c-TBI.tsv'

# Read in the data and store in dataframe
b_ctrl_gene_df = pd.read_csv(b_ctrl_gene_path, names=["Gene #", 'Gene Name', 'etc'], sep='\t')

# Assign values to each gene name for matching in the UMI matrix
b_ctrl_gene_df['Count'] = range(1, len(b_ctrl_gene_df) + 1)

# Convert 'Count' column to object type
b_ctrl_gene_df['Count'] = b_ctrl_gene_df['Count'].astype(str)

# Merge the two DataFrames based on the 'Gene' column from the first DataFrame and the 'Count' column from the second
# DataFrame
merged_df_ctrl = b_ctrl_df.merge(b_ctrl_gene_df[['Gene Name', 'Count']], left_on='Gene', right_on='Count', how='left')

# Create a DataFrame filled with zeros
zeros_df_1 = pd.DataFrame(np.zeros((7834, len(b_ctrl_gene_df['Gene Name']))))

# Define column labels
colum = b_ctrl_gene_df['Gene Name']
# Define row labels
ind = range(1, 7835)

# Set column and row labels
zeros_df_1.columns = colum
zeros_df_1.index = ind
# Add a label to the index
zeros_df_1 = zeros_df_1.rename_axis('Cell #')

# Convert 'Cell #' column to int type
merged_df_ctrl['Cell #'] = merged_df_ctrl['Cell #'].astype(int)

# Pivot the second matrix DataFrame to get the desired format
pivoted_df_ctrl = merged_df_ctrl.pivot_table(index='Cell #', columns='Gene Name', values='UMI', fill_value=0,
                                             aggfunc='sum')

# Reindex the pivoted DataFrame to match the index of the first matrix DataFrame (if needed)
pivoted_df_ctrl = pivoted_df_ctrl.reindex_like(zeros_df_1)

# Fill missing values with 0 (if there are cells in the first matrix DataFrame not present in the second)
pivoted_df_ctrl = pivoted_df_ctrl.fillna(0)

# Print the pivoted DataFrame
print(pivoted_df_ctrl)

# **Data Quality Control and Normalization:**
# With the read-in and formatted data tables, run quality control and normaliztion for cell cluster visualization




# set parameters for scanpy figures generated
sc.settings.set_figure_params(dpi=150, facecolor="white", figsize=(10, 8))

# In[5]:


# Convert all values in the DataFrame to integers
pivoted_df = pivoted_df.astype(int)

# Convert all values in the DataFrame to integers
pivoted_df_ctrl = pivoted_df_ctrl.astype(int)

# empty array to combine blast and control samples
adatas = {}

# Create an AnnData object for blast object
blast_adata = sc.AnnData(pivoted_df)
blast_adata.var_names_make_unique()
adatas["blast"] = blast_adata  # label those samples blast for further analysis

# Create an AnnData object for blast ctrl object
blast_ctrl_adata = sc.AnnData(pivoted_df_ctrl)
blast_ctrl_adata.var_names_make_unique()
adatas["control"] = blast_ctrl_adata  # label those samples control for further analysis

# concatenate dataframes to get one large anndata onject for scanpy analysis
adata = ad.concat(adatas, label="sample")
adata.obs_names_make_unique()

# Print the AnnData object
print(adata)
print(adata.obs["sample"])

# In[6]:


# REMOVE MITOCHONDRIAL GENES AND DEAD CELLS SO ACCURATE ANALYSIS

# mitochondrial genes, "mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("mt-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

# filter cells with less than 100 genes expressed and genes that are detected in less than 3 cells
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)

# filter for 10% mitochondrial genes
threshold_mt = 10
mito_filter = adata.obs['pct_counts_mt'] < threshold_mt
adata_blast_filtered = adata[mito_filter]

# observe qc metrics jointly - can remove cells that have too many mitochondrial genes expressed or total counts
sc.pl.violin(
    adata_blast_filtered,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)

# In[7]:


# DATA NORMALIZATION

# Saving count data
adata_blast_filtered.layers["counts"] = adata_blast_filtered.X.copy()
# Normalizing to median total counts
sc.pp.normalize_total(adata_blast_filtered)
# Logarithmize the data
sc.pp.log1p(adata_blast_filtered)

# **Analysis:**
# Dimensional reduction, clustering, and visulization of expression data.

# In[36]:


# reduce the dimensionality need the most informative genes - annotate highly variable genes
sc.pp.highly_variable_genes(adata_blast_filtered, n_top_genes=2000, batch_key="sample")

# reduce the dimensionality of the data by running PCA
sc.tl.pca(adata_blast_filtered, n_comps=50)

# compute neighborhood graph based on gene expression data
sc.pp.neighbors(adata_blast_filtered)

# graph can be embedded in two dimensions for UMAP visualization according to blast vs control sample type
sc.tl.umap(adata_blast_filtered)

# plot UMAP visualization in three different resolutions
for res in [0.02, 0.15, 2.0]:
    sc.tl.leiden(
        adata_blast_filtered, key_added=f"leiden_res_{res:4.2f}", resolution=res
    )

# In[43]:


# plot umap embedding of scRNA data at different resolutions and highlighted by TBI vs control sample and dave png to figures folder
sc.pl.umap(
    adata_blast_filtered,
    color=["sample", "leiden_res_0.15"],
    title=["UMAP Plot Highlighted by Experimental Condition", "UMAP Plot Highlighted by Leiden Clustering"],
    save='..\\figures\\umap_visuliazation_clusters.png'
)