#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_squared_error, r2_score

currentpath = os.getcwd()
print("Current path:", currentpath)


# In[ ]:





# In[5]:


print("CPM Normalization")


# In[6]:


# we DO normalize the data to the number of counts per column (CPM)


# In[7]:


file_path = 'A549_featureCounts_table.txt'

# Read the file using tab as the separator
df = pd.read_csv(file_path, sep='\t')

# Take a quick look at the DataFrame
print(df.head(3))


# In[8]:


rows, columns = df.shape
print("Number of genes :", rows)
print("Number of samples:", columns - 6 )


# In[9]:


print(df.columns)


# In[10]:


i = df 
# i['Chrom'] = i['Chr'].apply(lambda s: list(set(filter(None, s.split(';'))))[0])
# i.head(3)
i['Chrom'] = i['Chr'].apply(lambda s: set(filter(None, s.split(';'))).pop())
i.head(3)

# Check the unique chromosome values in the new 'Chrom' column
print("Unique Chrom values:", i['Chrom'].unique())
print("Number of rows:", i.shape[0])
print("Number of columns:", i.shape[1])


# In[11]:


print("correlation matrix")


# In[52]:


correlation_matrix = i[['AC1_IP', 'AC1_NP', 'AC2_IP', 'AC2_NP', 'AP1_IP', 'AP1_NP',
                        'AP2_IP', 'AP2_NP', 'AV1_IP', 'AV1_NP', 'AV2_IP', 'AV2_NP']].corr()
print(correlation_matrix)

# Display correlation matrix as heatmap
plt.figure(figsize=(8, 8))
ax = sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f", linewidths=0.5)
plt.title("Correlation Matrix Heatmap")

# Rotate y-axis labels to 90 degrees
ax.set_yticklabels(ax.get_yticklabels(), rotation=90)

plt.show()


# In[ ]:





# In[53]:


# Compute the number of NA and NaN :

total_nans = i.isna().sum().sum()
print("Total NaN values:", total_nans)

nans_per_column = i.isna().sum()
# print("NaNs per column:\n", nans_per_column)
print("Shape before removing NaNs:", i.shape)


# In[54]:


i_clean = i.dropna()
print("Shape after removing NaNs:", i_clean.shape)
i = i_clean


# In[ ]:





# In[55]:


# I preferred to normaliz after excluding these chromosomes (chrY, and chrM)


# In[56]:


fi = i[~i["Chrom"].isin(["chrY", "chrM"])]


# In[57]:


print("Unique Chrom values:", fi['Chrom'].unique())
print("Number of rows:", fi.shape[0])
print("Number of columns:", fi.shape[1])

print("Name of columns: \n")
print(fi.columns.tolist())


# In[ ]:





# In[58]:


fi2 = fi.drop(columns=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'Chrom'])


# In[59]:


# Get summary statistics (mean, min, max) for numeric columns
summary_stats = fi2.describe().loc[['mean', 'min', 'max']]

# Compute the median for numeric columns and add it to the summary statistics
summary_stats.loc['median'] = fi2.median()

# Optional: reorder the rows so that the median appears after the mean
summary_stats = summary_stats.loc[['mean', 'median', 'min', 'max']]

print(summary_stats.round(2))

print("Name of columns: \n")
print(fi2.columns.tolist())


# In[60]:


print("keep only rows with NO zeros in these columns")


# In[61]:


# Extract column names
cols = fi2.columns.tolist()
print("Column names:", cols)

# Filter rows: keep only rows with no zeros in these columns
fi2_filtered = fi2[(fi2[cols] != 0).all(axis=1)]
fi2 = fi2_filtered

# Verify the results
print("Filtered DataFrame shape:", fi2_filtered.shape)


# In[62]:


print("CPM : the normalization per TOTAL COUNTS per each COLUMN")


# In[63]:


# Let's do the normalization per TOTAL COUNTS per each COLUMN

# Normalize selected columns by dividing by their own column sums (total counts normalization)
fi2n = fi2[cols].div(fi2[cols].sum(), axis=1) * 1e6  # Multiply by 1 million (CPM normalization)

# Display the normalized DataFrame
print(fi2n.head())

column_sums_after_norm = fi2n.sum()
print("Column sums after normalization:")
print(column_sums_after_norm)

# Verify the data frame
print("Filtered DataFrame shape:", fi2n.shape)


# In[64]:


print(cols)


# In[65]:


get_ipython().run_line_magic('matplotlib', 'inline')

import seaborn as sns
import matplotlib.pyplot as plt

# Set up the plot
plt.figure(figsize=(6, 3))

# Plot density curves for each sample
for col in cols:
    sns.kdeplot(fi2n[col], label=col, fill=True, alpha=0.4)

# Customize plot with smaller font sizes
plt.xlabel("Expression Level", fontsize=10)
plt.ylabel("Density", fontsize=10)
plt.title("Density Plots for all the Samples", fontsize=10)  
plt.legend(fontsize=8)
plt.xlim(0, 2000)
plt.grid(False)

# Display the plot
plt.show()


# In[66]:


# Verify the data frame
print("Filtered DataFrame shape:", fi2n.shape)

# Calculate summary statistics
summary_stats = fi2n.agg(['max', 'min', 'median']).T

# Display summary clearly
print("Summary statistics for each column:\n")
print(summary_stats.round(4))


# In[67]:


import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(8, 2))

# Create detailed boxplots with customized appearance
ax = sns.boxplot(
    data=fi2n,
    showmeans=False,
    medianprops={"color": "red", "linewidth": 2},
    boxprops={"facecolor": "lightblue", "edgecolor": "black", "linewidth": 1.5},
    whiskerprops={"color": "black", "linewidth": 1.5},
    capprops={"color": "black", "linewidth": 1.5},
    flierprops={"marker": "x", "markerfacecolor": "gray", "markeredgecolor": "gray", "markersize": 5}
)

# Annotate medians clearly
medians = fi2n.median()
for tick, median in enumerate(medians):
    ax.text(tick, median + 0.5, f'{median:.2f}', 
            horizontalalignment='center', 
            color='red', 
            weight='bold')

# Customize plot appearance
plt.xticks(rotation=45, fontsize=10)
plt.ylabel("Normalized Values", fontsize=10)
plt.title("Boxplots with Annotated Median Values", fontsize=10)
plt.ylim(0, 60)  # Adjust y-limit as needed

# Display plot
plt.show()


# In[68]:


rows, columns = fi2n.shape
print(f"Number of rows: {rows}")
print(f"Number of columns: {columns}")

print(fi2n.columns)


# In[ ]:





# In[69]:


# Filtering (different filtering criteria from the previous analysis, on raw values)
# we removed the rows where the sum was lower thn a threshold


# In[70]:


# Filtering genes with CPM >= 1 in at least 2 samples
fi2n_cpm = fi2n[(fi2n >= 1).sum(axis=1) >= 2]

print("Previous shape:", fi2n.shape)
print("Filtered shape:", fi2n_cpm.shape)


# In[71]:


print("working with the data frame fi2n_cpm")


# In[72]:


import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(8, 2))

# Create detailed boxplots with customized appearance
ax = sns.boxplot(
    data=fi2n_cpm,
    showmeans=False,
    medianprops={"color": "red", "linewidth": 2},
    boxprops={"facecolor": "lightblue", "edgecolor": "black", "linewidth": 1.5},
    whiskerprops={"color": "black", "linewidth": 1.5},
    capprops={"color": "black", "linewidth": 1.5},
    flierprops={"marker": "x", "markerfacecolor": "gray", "markeredgecolor": "gray", "markersize": 5}
)

# Annotate medians clearly
medians = fi2n_cpm.median()
for tick, median in enumerate(medians):
    ax.text(tick, median + 0.5, f'{median:.2f}', 
            horizontalalignment='center', 
            color='red', 
            weight='bold')

# Customize plot appearance
plt.xticks(rotation=45, fontsize=10)
plt.ylabel("Normalized CPM Values", fontsize=10)
plt.title("Boxplots with Annotated Median Values", fontsize=10)
plt.ylim(0, 60) 

# Display plot
plt.show()


# In[102]:


import seaborn as sns
import matplotlib.pyplot as plt

# Set plot size
plt.figure(figsize=(12, 4))

# Create violin plots
ax = sns.violinplot(
    data=fi2n_cpm,
    inner='box',  # clearly shows median, quartiles
    linewidth=1.5,
    palette="Pastel1"
)

# Annotate medians clearly
medians = fi2n_cpm.median()
for tick, median in enumerate(medians):
    plt.text(tick, median, f'{median:.2f}',
             horizontalalignment='center',
             color='red',
             weight='bold')

# Customize axes labels and title
plt.xlabel("Samples", fontsize=12)
plt.ylabel("CPM", fontsize=12)
plt.title("Violin Plots with Annotated Median Values", fontsize=12)
plt.xticks(rotation=45, fontsize=10)

# Optional y-axis limits
plt.ylim(bottom=0, top=1000)  # Adjust as needed

# Display plot
plt.grid(False)
plt.show()


# In[ ]:





# In[74]:


fi3 = fi2n_cpm


# In[75]:


print("the matrix fi3 contains the CPM values")


# In[76]:


stats3 = fi3[['AC1_IP', 'AC1_NP', 'AC2_IP', 'AC2_NP', 'AP1_IP', 'AP1_NP', 
              'AP2_IP', 'AP2_NP', 'AV1_IP', 'AV1_NP', 'AV2_IP', 'AV2_NP']].agg(['min', 'max', 'median'])

print(stats3)

print("Number of rows:", fi3.shape[0])
print("Number of columns:", fi3.shape[1])

print("Name of columns: \n")
print(fi3.columns.tolist())

fi3.to_csv("A549_featureCounts_table_cpm_values_afterf1cpm2samples.tsv", sep="\t", index=False)


# In[ ]:





# In[77]:


import seaborn as sns
import matplotlib.pyplot as plt

def plot_all_visualizations(df, col1, col2):
    """Generate Violin, ECDF, and Density plots for selected columns in a single row."""
    
    fig, axes = plt.subplots(1, 3, figsize=(9, 3))  # Creating a row with 3 subplots
    
    # Violin Plot
    melted_df = df[[col1, col2]].melt(var_name="Condition", value_name="Value")
    sns.violinplot(x="Condition", y="Value", data=melted_df, ax=axes[0])
    axes[0].set_ylim(0, 3000)
    axes[0].set_title(f"Violin Plot", fontsize=8)
    axes[0].set_xlabel("Condition", fontsize=8)
    axes[0].set_ylabel("Values", fontsize=8)

    # ECDF Plot
    sns.ecdfplot(data=df, x=col1, label=col1, ax=axes[1])
    sns.ecdfplot(data=df, x=col2, label=col2, ax=axes[1])
    axes[1].set_xlim(0, 300)
    axes[1].set_title(f"ECDF Plot", fontsize=8)
    axes[1].set_xlabel("Value", fontsize=8)
    axes[1].set_ylabel("ECDF", fontsize=8)
    axes[1].legend(fontsize=6)

    # Density Plot
    sns.kdeplot(data=df, x=col1, label=col1, fill=True, bw_adjust=1, alpha=0.5, ax=axes[2])
    sns.kdeplot(data=df, x=col2, label=col2, fill=True, bw_adjust=1, alpha=0.5, ax=axes[2])
    axes[2].set_xlim(0, 1000)
    axes[2].set_title(f"Density Plot", fontsize=8)
    axes[2].set_xlabel("Value", fontsize=8)
    axes[2].set_ylabel("Density", fontsize=8)
    axes[2].legend(fontsize=6)

    plt.tight_layout()
    plt.show()


# In[78]:


plot_all_visualizations(fi3, 'AC1_IP', 'AC2_NP')


# In[79]:


plot_all_visualizations(fi3, 'AC2_IP', 'AC2_NP')


# In[80]:


plot_all_visualizations(fi3, 'AV1_IP', 'AV1_NP')


# In[81]:


plot_all_visualizations(fi3, 'AV2_IP', 'AV2_NP')


# In[82]:


plot_all_visualizations(fi3, 'AP1_IP', 'AP2_NP')


# In[83]:


plot_all_visualizations(fi3, 'AP2_IP', 'AP2_NP')


# In[84]:


print(fi3.columns.tolist())


# In[85]:


# Plotting the correlation plots between the replicates


# In[46]:


import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, r2_score
import numpy as np


# In[47]:


def plot_correlation(df, col1, col2):
    """Generate a correlation scatter plot with regression line, R² value, and RMSE."""
    
    # Increase figure size for clarity
    plt.figure(figsize=(2, 2))
    
    # Drop NaN values from the selected columns
    df_clean = df[[col1, col2]].dropna()
    
    # Create scatter plot with regression line
    ax = sns.regplot(x=df_clean[col1], y=df_clean[col2],
                     scatter_kws={'alpha': 0.5}, line_kws={'color': 'red'})
    
    # Compute R² and RMSE
    r2 = r2_score(df_clean[col1], df_clean[col2])
    rmse = np.sqrt(mean_squared_error(df_clean[col1], df_clean[col2]))
    
    # Set labels, title, and tick font sizes
    plt.xlabel(col1, fontsize=10)
    plt.ylabel(col2, fontsize=10)
    plt.title(f"Correlation Plot: {col1} vs {col2}", fontsize=10)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    
    # Set axis limits if needed (adjust or remove if your data range is different)
    plt.xlim(0, 10000)
    plt.ylim(0, 10000)
    
    # Display R² and RMSE on the plot using axes coordinates
    plt.text(0.05, 0.95, f'R²: {r2:.4f}\nRMSE: {rmse:.2f}', 
             transform=plt.gca().transAxes, fontsize=10, 
             verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
    
    plt.show()


# In[48]:


plot_correlation(fi3, 'AC1_IP', 'AC2_IP')


# In[49]:


plot_correlation(fi3, 'AP1_IP', 'AP2_IP')


# In[50]:


plot_correlation(fi3, 'AV1_IP', 'AV2_IP')


# In[ ]:





# In[93]:





# In[115]:


print(fi3.head(3))


# In[116]:


print(fi3.tail(3))


# In[117]:


# Get the column names as a list
columns_list = fi3.columns.tolist()
print(columns_list)


# In[199]:


import pandas as pd
from sklearn.decomposition import PCA
# Standardize the data (recommended for PCA)
from sklearn.preprocessing import StandardScaler


# In[120]:


import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

# Handle missing values (example: remove rows with NaNs)
rna_data = fi3.dropna()

# Scale the data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(rna_data.T)

# 2. Perform PCA
pca = PCA(n_components=2)
principal_components = pca.fit_transform(scaled_data)

# 3. Visualize Results
# Create DataFrame for PCA results
pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])

# Add sample names to the PCA DataFrame
pca_df['Sample'] = rna_data.columns

# Scatter plot
plt.figure(figsize=(8, 8))
sns.scatterplot(x='PC1', y='PC2', hue='Sample', data=pca_df, s=100)

# Add sample labels as annotations
for i, sample_name in enumerate(pca_df['Sample']):
    plt.annotate(sample_name, (pca_df['PC1'][i], pca_df['PC2'][i]), fontsize=8, ha='center', va='bottom')

plt.xlabel('Principal Component 1 (PC1)')
plt.ylabel('Principal Component 2 (PC2)')
plt.title('PCA of CPM')
plt.grid(True)
plt.show()

# Explained variance ratio
print("Explained Variance Ratio:", pca.explained_variance_ratio_)


# In[86]:


rna_data = fi3.dropna()

# Standardizing the data (genes as rows, samples as columns)
# scaler = StandardScaler()
# scaled_data = scaler.fit_transform(rna_data)  # Keep genes as rows


# In[87]:


print("Heatmap")


# In[88]:


import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
import seaborn as sns
import matplotlib.pyplot as plt

# Handle missing values
rna_data = fi3.dropna()

# Standardize the data (Z-score normalization)
scaler = StandardScaler()
standardized_data = scaler.fit_transform(rna_data)

# Standardize the data (Z-score normalization)
scaler = StandardScaler()
standardized_data = scaler.fit_transform(rna_data)

# Clip values between -2 and 2
clipped_data = np.clip(standardized_data, -2, 2)

# Convert to DataFrame
clipped_df = pd.DataFrame(clipped_data, index=rna_data.index, columns=rna_data.columns)

### **1. Create Heatmap with Sample Dendrogram**
# Perform hierarchical clustering on samples (columns)
dist_matrix = pdist(clipped_df.T, metric='euclidean')   # Pairwise distances between columns
col_linkage = linkage(dist_matrix, method='ward')       # Use Ward's method for clustering

# Create the clustered heatmap
# plt.figure(figsize=(8, 8))

g = sns.clustermap(
    clipped_df,
    cmap='RdYlGn',
    col_linkage=col_linkage,
    figsize=(6, 6),
    vmin=-2,
    vmax=2,
    yticklabels=False,  # Turn off row labels
    # xticklabels=False  # (Optional) Turn off column labels if you wish
)

# Optionally set a title on the heatmap itself (use the heatmap axis)
g.ax_heatmap.set_title("Clipped Standardized Heatmap (Ward Linkage)", fontsize=10)

plt.show()

### **2. Create Separate Dendrogram for Samples**
# plt.figure(figsize=(5, 5))
# dendrogram(col_linkage, labels=list(rna_data.columns), leaf_rotation=90, leaf_font_size=6)

# plt.title('Hierarchical Clustering Dendrogram (Samples - Ward Linkage)', fontsize=10)  # Updated title
# plt.xlabel('Samples', fontsize=8)
# plt.ylabel('Distance', fontsize=8)
# plt.xticks(fontsize=6)
# plt.grid(False)
# plt.show()

# Perform hierarchical clustering on genes (rows)
# dist_matrix_rows = pdist(clipped_df, metric='euclidean')  # Pairwise distances between rows
# row_linkage = linkage(dist_matrix_rows, method='ward')  # Ward's method for rows

### **3. Create Separate Dendrogram for Genes**
# plt.figure(figsize=(5, 5))
# dendrogram(row_linkage, labels=list(rna_data.index), leaf_rotation=90, leaf_font_size=6)

# plt.title('Hierarchical Clustering Dendrogram (Genes - Ward Linkage)', fontsize=10)
# plt.xlabel('Genes', fontsize=8)
# plt.ylabel('Distance', fontsize=8)
# plt.xticks(fontsize=6)
# plt.show()


# In[89]:


print("MDS")


# In[90]:


get_ipython().run_line_magic('matplotlib', 'inline')

import pandas as pd
import numpy as np
from sklearn.manifold import MDS
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

# Handle missing values (example: remove rows with NaNs)
rna_data = fi3.dropna()

# Scale the data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(rna_data.T)

# 1. Perform MDS
mds = MDS(n_components=2, random_state=42)  # Set random_state for reproducibility
mds_coordinates = mds.fit_transform(scaled_data)

# 2. Visualize Results
# Create DataFrame for MDS results
mds_df = pd.DataFrame(data=mds_coordinates, columns=['MDS1', 'MDS2'])

# Add sample names to the MDS DataFrame
mds_df['Sample'] = rna_data.columns

# Scatter plot
plt.figure(figsize=(6, 6))
sns.scatterplot(x='MDS1', y='MDS2', hue='Sample', data=mds_df, s=100)

# Add sample labels as annotations
for i, sample_name in enumerate(mds_df['Sample']):
    plt.annotate(sample_name, (mds_df['MDS1'][i], mds_df['MDS2'][i]), fontsize=8, ha='center', va='bottom')

plt.xlabel('MDS Dimension 1')
plt.ylabel('MDS Dimension 2')
plt.title('MDS Plot of CPM')
plt.grid(True)
plt.show()


# In[145]:


print("UMAP")


# In[91]:


import pandas as pd
from sklearn.preprocessing import StandardScaler
import umap
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm

# Step 1: Remove rows with missing data (NaNs)
rna_data = fi3.dropna()

# Step 2: Scale the data (transpose so samples are rows)
scaler = StandardScaler()
scaled_array = scaler.fit_transform(rna_data.T)

# Convert scaled_array back to DataFrame (for easy labeling)
scaled_df = pd.DataFrame(scaled_array, index=rna_data.columns)

# Run UMAP
umap_model = umap.UMAP(n_components=2, random_state=42)
umap_components = umap_model.fit_transform(scaled_df)

# Generate a distinct color for each sample
num_samples = len(scaled_df.index)
colors = cm.tab20(np.linspace(0, 1, num_samples))  
# You can also try 'cm.rainbow', 'cm.tab10', etc.

# Plot each sample with a unique color
plt.figure(figsize=(6, 6))
for i, sample_name in enumerate(scaled_df.index):
    plt.scatter(umap_components[i, 0], umap_components[i, 1], 
                color=colors[i], alpha=0.7)
    plt.text(umap_components[i, 0], umap_components[i, 1], 
             sample_name, fontsize=9, alpha=0.8)

plt.title('UMAP of CPM data')
plt.xlabel('UMAP Dimension 1')
plt.ylabel('UMAP Dimension 2')
plt.grid(True)
plt.show()


# In[ ]:





# In[92]:


print("NMF")


# In[93]:


# NMF

import pandas as pd
import numpy as np
from sklearn.decomposition import NMF
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

# Handle missing values (example: remove rows with NaNs)
rna_data = fi3.dropna()

# Scale the data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(rna_data.T)

# 2. Perform NMF
# Note: NMF requires all input data to be non-negative, so ensure your data meets this requirement.
# If your scaled data has negative values, apply a transformation to make it non-negative.
non_negative_data = scaled_data - np.min(scaled_data)  # Shift values to be non-negative
nmf = NMF(n_components=2, random_state=42)
nmf_components = nmf.fit_transform(non_negative_data)

# 3. Visualize Results
# Create DataFrame for NMF results
nmf_df = pd.DataFrame(data=nmf_components, columns=['Component1', 'Component2'])

# Add sample names to the NMF DataFrame
nmf_df['Sample'] = rna_data.columns

# Scatter plot
plt.figure(figsize=(8, 8))
sns.scatterplot(x='Component1', y='Component2', hue='Sample', data=nmf_df, s=100)

# Add sample labels as annotations
for i, sample_name in enumerate(nmf_df['Sample']):
    plt.annotate(sample_name, (nmf_df['Component1'][i], nmf_df['Component2'][i]), fontsize=8, ha='center', va='bottom')

plt.xlabel('NMF Component 1')
plt.ylabel('NMF Component 2')
plt.title('NMF of CPM data')
plt.grid(True)
plt.show()


# In[ ]:


print("ICA")


# In[ ]:


# ICA


# In[94]:


get_ipython().run_line_magic('matplotlib', 'inline')

import pandas as pd
import numpy as np
from sklearn.decomposition import FastICA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

# Handle missing values (example: remove rows with NaNs)
rna_data = fi3.dropna()

# Scale the data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(rna_data.T)

# 2. Perform ICA
n_components = 2  # Choose the number of components
ica = FastICA(n_components=n_components, random_state=42)
ica_components = ica.fit_transform(scaled_data)

# 3. Visualize Results
# Create DataFrame for ICA results
ica_df = pd.DataFrame(data=ica_components, columns=[f'IC{i+1}' for i in range(n_components)])

# Add sample names to the ICA DataFrame
ica_df['Sample'] = rna_data.columns

# Scatter plot
plt.figure(figsize=(8, 6))
sns.scatterplot(x=f'IC1', y=f'IC2', hue='Sample', data=ica_df, s=100)
plt.legend(title='Sample Names')

# Add sample labels as annotations
for i, sample_name in enumerate(ica_df['Sample']):
    plt.annotate(sample_name, (ica_df[f'IC1'][i], ica_df[f'IC2'][i]), fontsize=8, ha='center', va='bottom', xytext=(0, 5), textcoords='offset points')

plt.xlabel('Independent Component 1 (IC1)')
plt.ylabel('Independent Component 2 (IC2)')
plt.title('ICA of CPM data')
plt.grid(False)
plt.show()

# If you want to see the components (genes' contributions):
# components_df = pd.DataFrame(ica.components_, columns=rna_data.index)
# print(components_df)


# In[95]:


print("LDA : in the previous documents")


# In[96]:


# using LDA


# In[97]:


# Chat's code

import pandas as pd
import numpy as np
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

# Handle missing values (example: remove rows with NaNs)
rna_data = fi3.dropna()

# Scale the data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(rna_data.T)

# 1. Prepare data for LDA
# For LDA, you need class labels for each sample.
# Assuming you have a list or Series called 'labels' with the same length as the number of samples:
# labels = [...]  # Replace with your actual labels

# If you don't have labels and just want to visualize the data with LDA, you can skip this part and
# set n_components to the number of features minus 1 (or fewer if desired).

# Example with simulated labels (replace with your actual labels if you have them)
# Here I'm assuming 3 different classes

n_samples = scaled_data.shape[0]
labels = np.random.choice(['Class 1', 'Class 2', 'Class 3'], size=n_samples)

# 2. Perform LDA
n_components = 2  # Choose the number of components (should be less than the number of classes)
lda = LinearDiscriminantAnalysis(n_components=n_components)
lda_components = lda.fit_transform(scaled_data, labels)  # Note: LDA requires labels for supervised dimensionality reduction

# 3. Visualize Results
# Create DataFrame for LDA results
lda_df = pd.DataFrame(data=lda_components, columns=[f'LD{i+1}' for i in range(n_components)])

# Add sample names and labels to the LDA DataFrame
lda_df['Sample'] = rna_data.columns
lda_df['Label'] = labels  # Add the labels

# Scatter plot
plt.figure(figsize=(8, 6))
sns.scatterplot(x=f'LD1', y=f'LD2', hue='Label', data=lda_df, s=100, style='Label')  # Use 'Label' for hue and style
plt.legend(title='Class Labels')

# Add sample labels as annotations
for i, sample_name in enumerate(lda_df['Sample']):
    plt.annotate(sample_name, (lda_df[f'LD1'][i], lda_df[f'LD2'][i]), fontsize=8, ha='center', va='bottom', xytext=(0, 5), textcoords='offset points')

plt.xlabel('Linear Discriminant 1 (LD1)')
plt.ylabel('Linear Discriminant 2 (LD2)')
plt.title('LDA of CPM Data')
plt.grid(True)
plt.show()

# Explained variance ratio
print("Explained Variance Ratio:", lda.explained_variance_ratio_)


# In[ ]:





# In[98]:


# use GMM = gaussian mixture models


# In[99]:


# use RBM = restricted Boltzman Machines


# In[100]:


print("RBM")


# In[165]:


import pandas as pd
import numpy as np
from sklearn.neural_network import BernoulliRBM
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import matplotlib.pyplot as plt
import seaborn as sns

# Handle missing values (example: remove rows with NaNs)
rna_data = fi3.dropna()

# Scale the data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(rna_data.T)

# 2. Perform RBM
n_components = 2  # Choose the number of components
rbm = BernoulliRBM(n_components=n_components, random_state=42, verbose=0)

# Create a pipeline to combine scaling and RBM
pipeline = Pipeline(steps=[('scaler', scaler), ('rbm', rbm)])

# Fit the pipeline to the data
pipeline.fit(scaled_data)

# Transform the data to get the hidden layer representation
rbm_components = pipeline.transform(scaled_data)

# 3. Visualize Results
# Create DataFrame for RBM results
rbm_df = pd.DataFrame(data=rbm_components, columns=[f'RBM{i+1}' for i in range(n_components)])

# Add sample names to the RBM DataFrame
rbm_df['Sample'] = rna_data.columns

# Scatter plot
plt.figure(figsize=(6, 6))
sns.scatterplot(x=f'RBM1', y=f'RBM2', hue='Sample', data=rbm_df, s=100)
plt.legend(title='Sample Names')

# Add sample labels as annotations
for i, sample_name in enumerate(rbm_df['Sample']):
    plt.annotate(sample_name, (rbm_df[f'RBM1'][i], rbm_df[f'RBM2'][i]), fontsize=8, ha='center', va='bottom', xytext=(0, 5), textcoords='offset points')

plt.xlabel('Restricted Boltzmann Machine Component 1 (RBM1)')
plt.ylabel('Restricted Boltzmann Machine Component 2 (RBM2)')
plt.title('RBM of CPM Data')
plt.grid(True)
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




