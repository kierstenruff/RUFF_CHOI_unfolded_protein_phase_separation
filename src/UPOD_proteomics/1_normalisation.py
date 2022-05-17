
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

def normalisation(input_folder, output_folder, sample_name):
    # Read in compiled ratios, including pre-labelled channels

    compiled = pd.read_excel(
        f'{input_folder}{sample_name}_Compiled.xlsx', sheet_name=None)
    compiled.update({key: value.drop([col for col in value.columns.tolist(
    ) if 'Unnamed: ' in col], axis=1) for key, value in compiled.items()})

    # Plot histograms to visualise variation per channel in total peptide abundance
    # collect individual data types
    peptides = compiled['Peptides'].copy()
    proteins = compiled['Proteins'].copy()

    # # melt dataframe to ease normalisation calculation
    # for df in [peptides, proteins]:

    #     info_cols = [col for col in df.columns.tolist() if col in ['Protein IDs', 'Sequence', 'Proteins']]
    #     sample_cols = [col for col in df.columns.tolist() if col not in info_cols]

    #     df = pd.melt(
    #         df, 
    #         id_vars=info_cols, 
    #         value_vars=sample_cols,
    #         var_name=sample_name, 
    #         value_name='abundance')

    # Visualise sum of abundances per plex --> based on equal loading theory should be equivalent
    fig, ax = plt.subplots(figsize=(15, 6))
    sns.barplot(data=peptides.groupby('sample_name').sum(
    ).reset_index(), x='sample_name', y='abundance')
    plt.xticks(rotation=45)
    plt.savefig(f'{output_folder}{sample_name}_pre-normalisation.png')
    plt.show()

    # Sum-normalisation according to peptide abundance(i.e. should have loaded same amount of total protein?) or cell number normalisation(total amount of protein should have been proportional to the number of cells to start with?)
    norm_factor = peptides.groupby(
        'sample_name').sum().reset_index().replace(0, np.nan)
    norm_factor['abundance'] = norm_factor['abundance'] / \
        norm_factor['abundance'][33] #normalising to highest abundance, to avoid compression of values downwards
    norm_factor = dict(norm_factor[['sample_name', 'abundance']].values)

    # Alternative normalisation via cell count
    norm_pooled_factor = peptides.groupby(
        'sample_name').sum().reset_index().replace(0, np.nan)
    norm_pooled_factor['replicate'] = norm_pooled_factor['sample_name'].str.split('_').str[1]
    references = norm_pooled_factor[norm_pooled_factor['sample_name'].str.contains('ref')]
    references = dict(references[['replicate', 'abundance']].values)
    norm_pooled_factor['reference'] = norm_pooled_factor['replicate'].map(references)
    norm_pooled_factor['abundance'] = norm_pooled_factor['abundance'] / \
        norm_pooled_factor['reference'] 
    norm_pooled_factor = dict(norm_pooled_factor[['sample_name', 'abundance']].values)
    
    # Apply normalisation factor to both peptide- and protein-level quantitation
    for df in [peptides, proteins]:
        df['sum_norm_factor'] = df['sample_name'].map(norm_factor)
        df['sum_norm'] = df['abundance'] / \
            df['sum_norm_factor']

        df['ref_norm_factor'] = df['sample_name'].map(norm_pooled_factor)
        df['ref_norm'] = df['abundance'] / \
            df['ref_norm_factor']


    # Assign sample, replicate descriptors
    peptides[['sample_type', 'replicate']] = peptides['sample_name'].str.split('_', expand=True)
    proteins[['sample_type', 'replicate']] = proteins['sample_name'].str.split('_', expand=True)

    # visualise abundance after normalisation
    fig, ax = plt.subplots(figsize=(15, 6))
    sns.barplot(data=peptides.groupby('sample_name').sum(
    ).reset_index(), x='sample_name', y='sum_norm')
    plt.xticks(rotation=45)
    plt.savefig(f'{output_folder}{sample_name}_post-normalisation.png')
    plt.show()

    # Save normalised results to excel
    FileHandling.df_to_excel(
        output_path=f'{output_folder}{sample_name}_normalised.xlsx',
        sheetnames=['Proteins', 'Peptides'],
        data_frames=[proteins, peptides]
    )


if __name__ == '__main__':
    input_folder = 'python_results/initial_cleanup/'
    output_folder = 'python_results/normalisation/'
    sample_names = ['AGG', 'SUP']

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    for sample in sample_names:
        normalisation(input_folder, output_folder, sample)