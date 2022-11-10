##############################################################################
##### SUPERVISED PROJECT BioStat M1 SSD
##############################################################################

#%% IMPORTING PACKAGES

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

#%%###########################################################################
##### DATA IMPORT
##############################################################################

path = "C:/Users/33638/Documents/MasterSSD/M1/S7/Projet_biostat/"

#%% DATA HEAD
df_head = pd.read_csv(path + "expression_data_tsg3_3686_samples_20982_genes__head.csv", sep = ";")
print(df_head)

#%% DATA TARGET
df_target = pd.read_csv(path + "expression_data_tsg3_3686_samples_20982_genes__targets.csv", sep = ";")
print(df_target)
print(df_target.columns)
#%% DATA CANCER
df_tissue_group = pd.read_csv(path + "expression_data_tcga_brca_TCGA-BRCA_log_fpkm_1250_samples_41779_genes.csv", sep = ";")
print(df_tissue_group)

#%% DATA TISSUE
df_tissue = pd.read_csv(path + "Tissue_specific_genes.csv", sep = ";")
df_tissue["id_gene"].round(0)
print(df_tissue)
print(df_tissue.columns)
print(df_tissue.loc[:, "tissue_L1"])
print(df_tissue.loc[:, "SNR_L1"])

#%% COMPLETE DATA
df_tsg3 = pd.read_csv(path + "expression_data_tsg3_3686_samples_20982_genes.csv", sep = ";")
print(df_tsg3)
print(df_tsg3.shape)

#%% ADD TISSU GROUP LEVEL 1 TO THE DATAFRAME
df_tsg3 = df_tsg3[sorted(df_tsg3)]
df_target = df_target.sort_values(by='id_sample')

tissue_lvl1 = list(df_target['tissue_group_level1'])
tissue_lvl1.extend([np.nan, np.nan])

df_tsg3.loc[len(df_tsg3)] = tissue_lvl1

#%% TEST ON THE FIRST GENE
test_gene1 = df_tsg3.iloc[[0,-1]]
gene = test_gene1['gene_symbols'][0]
test_gene1 = test_gene1.transpose()
test_gene1 = test_gene1.drop(['gene_symbols', 'id_gene'])
test_gene1.set_axis(['values', 'tissue'], axis=1, inplace = True)
test_gene1 = test_gene1.astype({'values': 'float'})
test_gene1_group = test_gene1.groupby('tissue')['values'].mean()

plt.bar(test_gene1_group.index, test_gene1_group.values)
plt.xticks(range(len(test_gene1_group.index)), test_gene1_group.index, rotation=90)
plt.title(gene)








