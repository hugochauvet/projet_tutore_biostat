##############################################################################
##### SUPERVISED PROJECT BioStat M1 SSD
##############################################################################

#%% IMPORTING PACKAGES

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

#%%###########################################################################
##### DATA IMPORT
##############################################################################

#%% DATA TARGET
def load_data_target(path):
    return pd.read_csv(path + "expression_data_tsg3_3686_samples_20982_genes__targets.csv", sep = ";")


#%% COMPLETE DATA
def load_data_tsg3(path):
    return  pd.read_csv(path + "expression_data_tsg3_3686_samples_20982_genes.csv", sep = ";")

#%% ADD TISSU GROUP LEVEL 1 TO THE DATAFRAME
def merge_tsg3_target(df_tsg3, df_target):
    df_tsg3 = df_tsg3[sorted(df_tsg3)]
    df_target = df_target.sort_values(by='id_sample')
    tissue_lvl1 = list(df_target['tissue_group_level1'])
    tissue_lvl1.extend([np.nan, np.nan])
    df_tsg3.loc[len(df_tsg3)] = tissue_lvl1
    return df_tsg3

#%% TEST ON THE FIRST GENE
def gene_expr_tissue(gene, df_merge):
    idx_gene = df_merge[df_merge['gene_symbols'] == gene].index
    test_gene1 = df_merge.iloc[[idx_gene[0],-1]]
    test_gene1 = test_gene1.transpose()
    test_gene1 = test_gene1.drop(['gene_symbols', 'id_gene'])
    test_gene1.set_axis(['values', 'tissue'], axis=1, inplace = True)
    test_gene1 = test_gene1.astype({'values': 'float'})
    test_gene1_group = test_gene1.groupby('tissue')['values'].mean() # mean normal expression
    max_value = test_gene1_group.nlargest(2)[0]
    second_value = test_gene1_group.nlargest(2)[1]
    SNR = round(max_value / second_value, 3)
    df_gene_plot = test_gene1_group.reset_index(level=0)
    # BARPLOT
    col_bar = ['grey' if x < max_value else 'red' if SNR > 3 else 'blue' for x in df_gene_plot['values']]
    ax = sns.barplot(data = df_gene_plot, x = "tissue", y = "values", palette = col_bar)
    ax.tick_params(axis = 'x', rotation = 90)
    ax.set_title(gene)
    col_SNR = 'red' if SNR >= 3 else 'black'
    ax.text(0.05, max_value*0.95, 'SNR = ' + str(SNR), fontsize = 16, color = col_SNR)

#%%###########################################################################
##### MAIN
##############################################################################

#%% LOAD DATA
path = "C:/Users/33638/Documents/MasterSSD/M1/S7/Projet_biostat/"
df_target = load_data_target(path)
df_tsg3 = load_data_tsg3(path)
df_merge = merge_tsg3_target(df_tsg3, df_target)

#%% EX4.1
print(df_merge['gene_symbols'])
gene_expr_tissue('CCDC158', df_merge)








