# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 17:23:47 2024

@author: teesh
"""

import numpy as np
import pandas as pd
from scipy.io import loadmat

CancerTypes = ['ACC','BLCA','BRCA','CESC',
        'CHOL','COAD','DLBC','ESCA',
        'GBM','HNSC','KICH','KIRC',
        'KIRP','LAML','LGG','LIHC',
        'LUAD','LUSC','MESO','OV',
        'PAAD','PCPG','PRAD','READ',
        'SARC','SKCM','STAD','TGCT',
        'THCA','THYM','UCEC','UCS','UVM'] # 33 cancer types

groups = ['WHITE','BLACK','NAT_A']

path_to_data = 'C:/Users/teesh/OneDrive - Indian Institute of Technology Guwahati/Dataset/EssentialData/'
path_to_GA = path_to_data + 'Genetic_Ancestry.xlsx'
path_to_Protein = path_to_data + 'ProteinData/Protein.txt'
path_to_mRNA = path_to_data + 'mRNAData/mRNA.mat'
path_to_Methylation = path_to_data + 'MethylationData/Methylation.mat'
path_to_MicroRNA = path_to_data + 'MicroRNAData/MicroRNA-Expression.mat'

# Protein data
def get_protein(CancerTypes, path_to_Protein):
    df = pd.read_csv(path_to_Protein, sep='\t', index_col='SampleID')
    df = df.dropna(axis=1)
    # tumorTypes = tumor_types(cancer_type)
    df = df[df['TumorType'].isin(CancerTypes)]
    df = df.drop(columns=['TumorType'])
    index = df.index.values
    index_new = [row[:12] for row in index]
    df.index = index_new
    
    return df

# mRNA - Expression data
def get_mRNA(CancerTypes, path_to_mRNA):
    A = loadmat(path_to_mRNA)
    X, Y, GeneName, SampleName = A['X'].astype('float32'), A['Y'], A['GeneName'][0], A['SampleName']
    GeneName = [row[0] for row in GeneName]
    SampleName = [row[0][0] for row in SampleName]
    Y = [row[0][0] for row in Y]
    df_X = pd.DataFrame(X, columns=GeneName, index=SampleName)
    df_Y = pd.DataFrame(Y, index=SampleName, columns=['Disease'])
    df_Y = df_Y[df_Y['Disease'].isin(CancerTypes)]
    df = df_X.join(df_Y, how='inner')
    df = df.drop(columns=['Disease'])
    index = df.index.values
    index_new = [row[:12] for row in index]
    df.index = index_new
    df = df.reset_index().drop_duplicates(subset='index', keep='first').set_index('index')
    
    return df

# DNA Methylation data
def get_Methylation(CancerTypes, path_to_Methylation,groups):
    MethylationData = loadmat(path_to_Methylation)
    # extracting input combinations data...
    X, Y, GeneName, SampleName = MethylationData['X'].astype('float32'), MethylationData['CancerType'], MethylationData['FeatureName'][0], MethylationData['SampleName']
    GeneName = [row[0] for row in GeneName]
    SampleName = [row[0][0] for row in SampleName]
    Y = [row[0][0] for row in Y]
    MethylationData_X = pd.DataFrame(X, columns=GeneName, index=SampleName)
    MethylationData_Y = pd.DataFrame(Y, index=SampleName, columns=['Disease'])
    MethylationData_Y = MethylationData_Y[MethylationData_Y['Disease'].isin(CancerTypes)]
    MethylationData_in = MethylationData_X.join(MethylationData_Y, how='inner')
    MethylationData_in = MethylationData_in.drop(columns=['Disease'])
    index = MethylationData_in.index.values
    index_new = [row[:12] for row in index]
    MethylationData_in.index = index_new
    MethylationData_in = MethylationData_in.reset_index().drop_duplicates(subset='index', keep='first').set_index('index')
    # adding ethnic group information...
    MethyAncsDataPath = path_to_data + 'MethylationData/MethylationGenetic.xlsx'
    #print(MethyAncsDataPath)
    # fetching ethnic_group info from MethylationGenetic.xlsx
    MethyAncsData = [pd.read_excel(MethyAncsDataPath,
                         disease, usecols='A,B',
                         index_col='bcr_patient_barcode',
                         keep_default_na=False)
           for disease in CancerTypes]
    MethyAncsData_ethnic_group = pd.concat(MethyAncsData)
    ethnic_group_groups = ['WHITE',
              'BLACK OR AFRICAN AMERICAN',
              'ASIAN',
              'AMERICAN INDIAN OR ALASKA NATIVE',
              'NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER']
    MethyAncsData_ethnic_group['ethnic_group'] = MethyAncsData_ethnic_group['race']
    MethyAncsData_ethnic_group = MethyAncsData_ethnic_group[MethyAncsData_ethnic_group['ethnic_group'].isin(ethnic_group_groups)]
    MethyAncsData_ethnic_group.loc[MethyAncsData_ethnic_group['ethnic_group'] == 'WHITE', 'ethnic_group'] = 'WHITE'
    MethyAncsData_ethnic_group.loc[MethyAncsData_ethnic_group['ethnic_group'] == 'BLACK OR AFRICAN AMERICAN', 'ethnic_group'] = 'BLACK'
    MethyAncsData_ethnic_group.loc[MethyAncsData_ethnic_group['ethnic_group'] == 'ASIAN', 'ethnic_group'] = 'ASIAN'
    MethyAncsData_ethnic_group.loc[MethyAncsData_ethnic_group['ethnic_group'] == 'AMERICAN INDIAN OR ALASKA NATIVE', 'ethnic_group'] = 'NAT_A'
    MethyAncsData_ethnic_group.loc[MethyAncsData_ethnic_group['ethnic_group'] == 'NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER', 'ethnic_group'] = 'OTHER'
    MethyAncsData_ethnic_group = MethyAncsData_ethnic_group[MethyAncsData_ethnic_group['ethnic_group'].isin(groups)]
    # Keep patients with ethnic group information
    MethylationData_in = MethylationData_in.join(MethyAncsData_ethnic_group, how='inner')
    MethylationData_in = MethylationData_in.dropna(axis='columns')
    
    return MethylationData_in

# MicroRNA Expression data
def get_MicroRNA(CancerTypes, path_to_MicroRNA):
    A = loadmat(path_to_MicroRNA)
    X, Y, GeneName, SampleName = A['X'].astype('float32'), A['CancerType'], A['FeatureName'][0], A['SampleName']
    GeneName = [row[0] for row in GeneName]
    SampleName = [row[0][0] for row in SampleName]
    Y = [row[0][0] for row in Y]
    df_X = pd.DataFrame(X, columns=GeneName, index=SampleName)
    df_Y = pd.DataFrame(Y, index=SampleName, columns=['Disease'])
    df_Y = df_Y[df_Y['Disease'].isin(CancerTypes)]
    df = df_X.join(df_Y, how='inner')
    df = df.drop(columns=['Disease'])
    index = df.index.values
    index_new = [row[:12] for row in index]
    df.index = index_new
    df = df.reset_index().drop_duplicates(subset='index', keep='first').set_index('index')
    
    return df

def get_ethnic_group(CancerTypes, path_to_GA, df, groups):
    df_list = [pd.read_excel(path_to_GA, disease, usecols='A,E', 
                    index_col='Patient_ID', keep_default_na=False)
                   for disease in CancerTypes]
    df_group = pd.concat(df_list)
    df_group = df_group[df_group['EIGENSTRAT'].isin(['EA', 'AA', 'EAA', 'NA', 'OA'])]
    df_group['ethnic_group'] = df_group['EIGENSTRAT']
    
    df_group.loc[df_group['EIGENSTRAT'] == 'EA', 'ethnic_group'] = 'WHITE'
    df_group.loc[df_group['EIGENSTRAT'] == 'AA', 'ethnic_group'] = 'BLACK'
    df_group.loc[df_group['EIGENSTRAT'] == 'EAA', 'ethnic_group'] = 'ASIAN'
    df_group.loc[df_group['EIGENSTRAT'] == 'NA', 'ethnic_group'] = 'NAT_A'
    df_group.loc[df_group['EIGENSTRAT'] == 'OA', 'ethnic_group'] = 'OTHER'
    df_group = df_group.drop(columns=['EIGENSTRAT'])
    
    df_group = df_group[df_group['ethnic_group'].isin(groups)]
    
    # Keep patients with ethnic_group information
    df = df.join(df_group, how='inner')
    #print(df.shape)
    df = df.dropna(axis='columns')
    
    return df

df = get_protein(CancerTypes, path_to_Protein)
df = get_ethnic_group(CancerTypes, path_to_GA, df, groups)
print('Protein')
print(np.shape(df))
print(np.unique(df['ethnic_group'],return_counts=True))

df = get_mRNA(CancerTypes, path_to_mRNA)
df = get_ethnic_group(CancerTypes, path_to_GA, df, groups)
print('mRNA')
print(np.shape(df))
print(np.unique(df['ethnic_group'],return_counts=True))

df = get_MicroRNA(CancerTypes, path_to_MicroRNA)
df = get_ethnic_group(CancerTypes, path_to_GA, df, groups)
print('MicroRNA')
print(np.shape(df))
print(np.unique(df['ethnic_group'],return_counts=True))

df = get_Methylation(CancerTypes, path_to_Methylation, groups)
print('Methylation')
print(np.shape(df))
print(np.unique(df['ethnic_group'],return_counts=True))

df = get_mRNA(CancerTypes, path_to_mRNA)
df = get_ethnic_group(CancerTypes, path_to_GA, df, groups)
print('mRNA')
print(np.shape(df))
unique_values, counts = np.unique(df['ethnic_group'], return_counts=True)
percentages = counts / np.sum(counts) * 100
print("Unique values and their percentages in mRNA:")
for value, count, percentage in zip(unique_values, counts, percentages):
    print(f"{value}: {count} ({percentage:.1f}%)")

df = get_MicroRNA(CancerTypes, path_to_MicroRNA)
df = get_ethnic_group(CancerTypes, path_to_GA, df, groups)
print('MicroRNA')
print(np.shape(df))
unique_values, counts = np.unique(df['ethnic_group'], return_counts=True)
percentages = counts / np.sum(counts) * 100
print("Unique values and their percentages in MicroRNA:")
for value, count, percentage in zip(unique_values, counts, percentages):
    print(f"{value}: {count} ({percentage:.1f}%)")

df = get_Methylation(CancerTypes, path_to_Methylation, groups)
print('Methylation')
print(np.shape(df))
unique_values, counts = np.unique(df['ethnic_group'], return_counts=True)
percentages = counts / np.sum(counts) * 100
print("Unique values and their percentages in Methylation:")
for value, count, percentage in zip(unique_values, counts, percentages):
    print(f"{value}: {count} ({percentage:.1f}%)")








