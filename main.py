##############################################################################
##### SUPERVISED PROJECT BioStat M1 SSD
##############################################################################

#%% IMPORTING PACKAGES

import numpy as np
import pandas as pd
import os

#%% DATA IMPORT

path = "C:/Users/33638/Documents/MasterSSD/M1/S7/Projet_biostat/"
df = pd.read_csv(path + "expression_data_tsg3_3686_samples_20982_genes__head.csv", sep = ";", index_col = 0)
print(df)
print(df.shape)
