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

#%% CONCATENATE DATA
df = df_tsg3.merge(df_tissue)

#%% GRAPH TO UNDERSTAND DATA
nb_tissue_L1 = df["tissue_L1"].value_counts()
print(nb_tissue_L1)