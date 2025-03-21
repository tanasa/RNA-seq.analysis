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





# In[2]:


# we did NOT normalize the data to the total number of counts per column


# In[ ]:





# In[3]:


file_path = 'A549_featureCounts_table.txt'

# Read the file using tab as the separator
df = pd.read_csv(file_path, sep='\t')

# Take a quick look at the DataFrame
print(df.head(3))


# In[4]:


rows, columns = df.shape
print("Number of genes :", rows)
print("Number of samples:", columns - 6 )


# In[5]:


print(df.columns)


# In[ ]:





# In[6]:


i = df
print(i.head(3))


# In[7]:


i = df 
# i['Chrom'] = i['Chr'].apply(lambda s: list(set(filter(None, s.split(';'))))[0])
# i.head(3)
i['Chrom'] = i['Chr'].apply(lambda s: set(filter(None, s.split(';'))).pop())
i.head(3)


# In[8]:


# Check the unique chromosome values in the new 'Chrom' column
print("Unique Chrom values:", i['Chrom'].unique())
print("Number of rows:", i.shape[0])
print("Number of columns:", i.shape[1])


# In[9]:


# Exclude the following chromosomes from the analysis : chrX, chrY, chrM


# In[10]:


fi = i[~i["Chrom"].isin(["chrX","chrY", "chrM"])]


# In[11]:


print("Unique Chrom values:", fi['Chrom'].unique())
print("Number of rows:", fi.shape[0])
print("Number of columns:", fi.shape[1])

print("Name of columns: \n")
print(fi.columns.tolist())


# In[12]:


fi2 = fi.drop(columns=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'Chrom'])


# In[13]:


# Get summary statistics (mean, min, max) for numeric columns
summary_stats = fi2.describe().loc[['mean', 'min', 'max']]
print(summary_stats.round(2))

print("Name of columns: \n")
print(fi2.columns.tolist())


# In[14]:


#  Filtering based on the sum of rows is lower than 1000


# In[15]:


# Compute the sum of numerical values in each row
fi2['Row_Sum'] = fi2.select_dtypes(include=['number']).sum(axis=1)

# Remove rows where the sum is lower than 1000
fi3 = fi2[fi2['Row_Sum'] >= 1000].drop(columns=['Row_Sum'])


# In[16]:


# Remove rows where the sum is lower than 1000


# In[17]:


# Get summary statistics (mean, min, max) for numeric columns
summary_stats3 = fi3.describe().loc[['mean', 'min', 'max']]
print(summary_stats3.round(2))

print("Name of columns: \n")
print(fi3.columns.tolist())


# In[18]:


stats3 = fi3[['AC1_IP', 'AC1_NP', 'AC2_IP', 'AC2_NP', 'AP1_IP', 'AP1_NP', 
              'AP2_IP', 'AP2_NP', 'AV1_IP', 'AV1_NP', 'AV2_IP', 'AV2_NP']].agg(['min', 'max', 'median'])

print(stats3)


# In[19]:


print("Number of rows:", fi3.shape[0])
print("Number of columns:", fi3.shape[1])

print("Name of columns: \n")
print(fi3.columns.tolist())

fi3.to_excel('A549_featureCounts_table_raw_counts_filter_sum_1000.xlsx', index=False)
fi3.to_csv('A549_featureCounts_table_raw_counts_filter_sum_1000.xls', sep="\t", index=False)


# In[ ]:





# In[20]:


import seaborn as sns
import matplotlib.pyplot as plt

def plot_all_visualizations(df, col1, col2):
    """Generate Violin, ECDF, and Density plots for selected columns in a single row."""
    
    fig, axes = plt.subplots(1, 3, figsize=(8, 2))  # Creating a row with 3 subplots
    
    # Violin Plot
    melted_df = df[[col1, col2]].melt(var_name="Condition", value_name="Value")
    sns.violinplot(x="Condition", y="Value", data=melted_df, ax=axes[0])
    axes[0].set_ylim(0, 100000)
    axes[0].set_title(f"Violin Plot", fontsize=8)
    axes[0].set_xlabel("Condition", fontsize=8)
    axes[0].set_ylabel("Values", fontsize=8)

    # ECDF Plot
    sns.ecdfplot(data=df, x=col1, label=col1, ax=axes[1])
    sns.ecdfplot(data=df, x=col2, label=col2, ax=axes[1])
    axes[1].set_xlim(0, 30000)
    axes[1].set_title(f"ECDF Plot", fontsize=8)
    axes[1].set_xlabel("Value", fontsize=8)
    axes[1].set_ylabel("ECDF", fontsize=8)
    axes[1].legend(fontsize=6)

    # Density Plot
    sns.kdeplot(data=df, x=col1, label=col1, fill=True, bw_adjust=1, alpha=0.5, ax=axes[2])
    sns.kdeplot(data=df, x=col2, label=col2, fill=True, bw_adjust=1, alpha=0.5, ax=axes[2])
    axes[2].set_xlim(0, 30000)
    axes[2].set_title(f"Density Plot", fontsize=8)
    axes[2].set_xlabel("Value", fontsize=8)
    axes[2].set_ylabel("Density", fontsize=8)
    axes[2].legend(fontsize=6)

    plt.tight_layout()
    plt.show()


# In[21]:


plot_all_visualizations(fi3, 'AC1_IP', 'AC2_NP')


# In[22]:


plot_all_visualizations(fi3, 'AC2_IP', 'AC2_NP')


# In[23]:


plot_all_visualizations(fi3, 'AV1_IP', 'AV1_NP')


# In[24]:


plot_all_visualizations(fi3, 'AV2_IP', 'AV2_NP')


# In[25]:


plot_all_visualizations(fi3, 'AP1_IP', 'AP2_NP')


# In[26]:


plot_all_visualizations(fi3, 'AP2_IP', 'AP2_NP')


# In[27]:


plt.figure(figsize=(12, 2))

# Reshaping data for seaborn
melted_df = fi3.melt(var_name="Condition", value_name="Value")

# Plot violin plot
sns.violinplot(x="Condition", y="Value", data=melted_df)

# Adding labels and title
plt.title("Experimental Conditions : Raw Values")
plt.ylabel("Raw Values")
plt.xlabel("Conditions")

# Setting y-axis limit
plt.ylim(0, 0.05e6)

# Display the plot
plt.show()


# In[28]:


print(fi3.columns.tolist())


# In[ ]:





# In[29]:


# Plotting the correlation plots between the replicates


# In[30]:


import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, r2_score
import numpy as np


# In[31]:


def plot_correlation(df, col1, col2):
    """Generate a correlation scatter plot with regression line, R² value, and RMSE."""

    plt.figure(figsize=(2, 2))

    # Drop NaN values from the selected columns
    df_clean = df[[col1, col2]].dropna()

    plt.figure(figsize=(2, 2))
    
    # Scatter plot with regression line
    sns.regplot(x=df_clean[col1], y=df_clean[col2], scatter_kws={'alpha': 0.5}, line_kws={'color': 'red'})
    
    # Compute R² and RMSE
    r2 = r2_score(df_clean[col1], df_clean[col2])
    rmse = np.sqrt(mean_squared_error(df_clean[col1], df_clean[col2]))
    
    # Set axis limits
    plt.xlim(0, 0.1 * 1e6)
    plt.ylim(0, 0.1 * 1e6)
    
    # Labels and formatting
    plt.xlabel(col1, fontsize=8)
    plt.ylabel(col2, fontsize=8)
    plt.title(f"Correlation Plot: {col1} vs {col2}", fontsize=8)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    
    # Display R² and RMSE on the plot
    plt.text(0.05 * 1e6, 0.045 * 1e6, f'R²: {r2:.4f}\nRMSE: {rmse:.2f}', 
             fontsize=8, bbox=dict(facecolor='white', alpha=0.5))

    plt.show()


# In[95]:


get_ipython().run_line_magic('matplotlib', 'inline')
plot_correlation(fi3, 'AC1_IP', 'AC2_IP')


# In[96]:


get_ipython().run_line_magic('matplotlib', 'inline')
plot_correlation(fi3, 'AP1_IP', 'AP2_IP')


# In[97]:


get_ipython().run_line_magic('matplotlib', 'inline')
plot_correlation(fi3, 'AV1_IP', 'AV2_IP')


# In[ ]:





# In[98]:


# PCA plot


# In[99]:


print(fi3.head(3))
print(fi3.tail(3))


# In[ ]:





# In[100]:


# Get the column names as a list
columns_list = fi3.columns.tolist()
print(columns_list)


# In[ ]:





# In[101]:


import pandas as pd
from sklearn.decomposition import PCA
# Standardize the data (recommended for PCA)
from sklearn.preprocessing import StandardScaler


# In[103]:


get_ipython().run_line_magic('matplotlib', 'inline')

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

# Handle missing values (example: remove rows with NaNs)
rna_data = fi3.dropna()

# 1. Scale the data
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
plt.figure(figsize=(7, 7))
sns.scatterplot(x='PC1', y='PC2', hue='Sample', data=pca_df, s=100)

# Add sample labels as annotations
for i, sample_name in enumerate(pca_df['Sample']):
    plt.annotate(sample_name, (pca_df['PC1'][i], pca_df['PC2'][i]), fontsize=8, ha='center', va='bottom')

plt.xlabel('Principal Component 1 (PC1)')
plt.ylabel('Principal Component 2 (PC2)')
plt.title('PCA')
plt.grid(False)
plt.show()

# Explained variance ratio
print("Explained Variance Ratio:", pca.explained_variance_ratio_)


# In[42]:


# from sklearn.preprocessing import StandardScaler
# from sklearn.decomposition import PCA

# Strategies for normalization :
# 1. Standardize the data (Z-score normalization)
# scaler = StandardScaler()
# standardized_data = scaler.fit_transform(rna_data)
# clipped_data = np.clip(standardized_data, -2, 2)
# 2. Min-Max Scaling (Alternative): not successful
# scaler = MinMaxScaler(feature_range=(-1,1))
# scaled_data = scaler.fit_transform(rna_data)


# In[ ]:





# In[49]:


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.preprocessing import StandardScaler

# Ensure rna_data is properly defined
rna_data = fi3.dropna()

# Standardize the data (Z-score normalization)
scaler = StandardScaler()
standardized_data = scaler.fit_transform(rna_data)

# Clip values between -2 and 2
clipped_data = np.clip(standardized_data, -2, 2)

# Convert to DataFrame
clipped_df = pd.DataFrame(clipped_data, index=rna_data.index, columns=rna_data.columns)

### **1. Create Heatmap with Sample Dendrogram**
# Perform hierarchical clustering on samples (columns)
dist_matrix = pdist(clipped_df.T, metric='euclidean')
col_linkage = linkage(dist_matrix, method='ward')  # pdist returns a condensed distance matrix

# Create a clustered heatmap without y-axis tick labels (names/numbers)
g = sns.clustermap(clipped_df, cmap='RdYlGn', col_linkage=col_linkage,
                   figsize=(4, 6), vmin=-2, vmax=2, yticklabels=False)
plt.suptitle('Euclidean Distance, Ward Linkage', fontsize=10)
plt.show()

### **2. Create Separate Dendrogram for Samples**

plt.figure(figsize=(4, 4))
dendrogram(col_linkage, labels=list(rna_data.columns), leaf_rotation=90, leaf_font_size=6)
plt.title('Dendrogram', fontsize=10)
plt.xlabel('Samples', fontsize=8)
plt.ylabel('Distance', fontsize=8)
plt.xticks(fontsize=6)
plt.show()


# In[ ]:





# In[104]:


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
get_ipython().run_line_magic('matplotlib', 'inline')

plt.figure(figsize=(6, 6))
sns.scatterplot(x='MDS1', y='MDS2', hue='Sample', data=mds_df, s=100)

# Add sample labels as annotations
for i, sample_name in enumerate(mds_df['Sample']):
    plt.annotate(sample_name, (mds_df['MDS1'][i], mds_df['MDS2'][i]), fontsize=8, ha='center', va='bottom')

plt.xlabel('MDS Dimension 1')
plt.ylabel('MDS Dimension 2')
plt.title('MDS')
plt.grid(True)
plt.show()


# In[113]:


print("UMAP")


# In[114]:


import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import umap
import matplotlib.pyplot as plt

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

# Generate a unique color for each sample using a colormap
num_points = len(scaled_df.index)
colors = plt.cm.hsv(np.linspace(0, 1, num_points))

# Plot UMAP with sample names as labels
plt.figure(figsize=(4, 4))
plt.scatter(umap_components[:, 0], umap_components[:, 1], c=colors, alpha=0.7)

# Annotate each sample
for i, sample_name in enumerate(scaled_df.index):
    plt.text(umap_components[i, 0], umap_components[i, 1], sample_name, fontsize=9, alpha=0.8)

plt.title('UMAP')
plt.xlabel('UMAP Dimension 1')
plt.ylabel('UMAP Dimension 2')
plt.grid(True)
plt.show()


# In[112]:


print("NMF")


# In[66]:


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
plt.figure(figsize=(10, 8))
sns.scatterplot(x='Component1', y='Component2', hue='Sample', data=nmf_df, s=100)

# Add sample labels as annotations
for i, sample_name in enumerate(nmf_df['Sample']):
    plt.annotate(sample_name, (nmf_df['Component1'][i], nmf_df['Component2'][i]), fontsize=8, ha='center', va='bottom')

plt.xlabel('NMF Component 1')
plt.ylabel('NMF Component 2')
plt.title('NMF')
plt.grid(True)
plt.show()


# In[ ]:





# In[67]:


print("ICA")


# In[93]:


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
plt.figure(figsize=(6, 6))
sns.scatterplot(x=f'IC1', y=f'IC2', hue='Sample', data=ica_df, s=100)
plt.legend(title='Sample Names')

# Add sample labels as annotations
for i, sample_name in enumerate(ica_df['Sample']):
    plt.annotate(sample_name, (ica_df[f'IC1'][i], ica_df[f'IC2'][i]), fontsize=8, ha='center', va='bottom', xytext=(0, 5), textcoords='offset points')

plt.xlabel('Independent Component 1 (IC1)')
plt.ylabel('Independent Component 2 (IC2)')
plt.title('ICA')
plt.grid(False)
plt.show()

# If you want to see the components (genes' contributions):
# components_df = pd.DataFrame(ica.components_, columns=rna_data.index)
# print(components_df)


# In[ ]:





# In[72]:


print("LDA")


# In[94]:


# Chat code
# Gemini's code

get_ipython().run_line_magic('matplotlib', 'inline')

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
plt.figure(figsize=(6, 6))
sns.scatterplot(x=f'LD1', y=f'LD2', hue='Label', data=lda_df, s=100, style='Label')  # Use 'Label' for hue and style
plt.legend(title='Class Labels')

# Add sample labels as annotations
for i, sample_name in enumerate(lda_df['Sample']):
    plt.annotate(sample_name, (lda_df[f'LD1'][i], lda_df[f'LD2'][i]), fontsize=8, ha='center', va='bottom', xytext=(0, 5), textcoords='offset points')

plt.xlabel('Linear Discriminant 1 (LD1)')
plt.ylabel('Linear Discriminant 2 (LD2)')
plt.title('LDA')
plt.grid(True)
plt.show()

# Explained variance ratio
print("Explained Variance Ratio:", lda.explained_variance_ratio_)


# In[ ]:





# In[ ]:


# use GMM = gaussian mixture models


# In[ ]:


# use RBM = restricted Boltzman Machines


# In[80]:


print("Restricted Boltzman Machines")


# In[92]:


get_ipython().run_line_magic('matplotlib', 'inline')

from adjustText import adjust_text
import pandas as pd
import numpy as np
from sklearn.neural_network import BernoulliRBM
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import matplotlib.pyplot as plt
import seaborn as sns

# Example: assume fi3 is your DataFrame of shape (genes, samples)
# Drop rows with NaNs
rna_data = fi3.dropna()

# Transpose so rows = samples, columns = features (as scikit-learn expects)
data_for_rbm = rna_data.T

# Define the pipeline with a scaler and RBM
n_components = 2  # we want a 2D representation
pipeline = Pipeline([
    ('scaler', StandardScaler()), 
    ('rbm', BernoulliRBM(n_components=n_components, random_state=42, verbose=0))
])

# Fit the pipeline to the data (this fits both scaler and RBM)
pipeline.fit(data_for_rbm)

# Transform the data to get the 2D RBM components
rbm_components = pipeline.transform(data_for_rbm)

# Build a DataFrame for easy plotting
rbm_df = pd.DataFrame(
    data=rbm_components, 
    columns=[f'RBM{i+1}' for i in range(n_components)]
)

# Add sample names to the DataFrame (these come from the original columns of fi3)
rbm_df['Sample'] = rna_data.columns

# Plot the 2D RBM representation
plt.figure(figsize=(10, 8))

fig, ax = plt.subplots(figsize=(10, 8))

# Create the scatter plot
sns.scatterplot(x='RBM1', y='RBM2', hue='Sample', data=rbm_df, s=100, ax=ax)
ax.set_xlabel('RBM Component 1')
ax.set_ylabel('RBM Component 2')
ax.set_title('RBM with Sample Labels Table')

# Create table data (you can customize which columns you want)
table_data = rbm_df[['Sample', 'RBM1', 'RBM2']]
# Add the table to the plot; here we place it to the right of the axes
table = plt.table(cellText=table_data.values,
                  colLabels=table_data.columns,
                  cellLoc='center',
                  loc='right')

# Optionally, adjust the table's appearance and scaling
table.scale(1, 1.5)  # Increase row height for readability
plt.subplots_adjust(right=0.7)  # Make space on the right for the table

plt.show()


# In[ ]:




