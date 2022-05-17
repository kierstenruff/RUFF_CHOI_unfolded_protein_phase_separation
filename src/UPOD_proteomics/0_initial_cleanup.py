"""Initial cleanup script to parse MaxQuant results for peptide and protein quantitation, collect columns of interest and assign sample labels to individual TMT channels"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')


def tmt_peptides(input_path, sample_names=None, pep_cols=[]):

    peptides = pd.read_table(input_path, sep='\t')
    logger.info(
        f'Imported peptides from {input_path}. {peptides.shape[0]} entries found.')
    cleaned_dfs = {}
    standard_cols = ['Sequence', 'Proteins', 'Gene names',
                     'Protein names', 'Unique (Groups)', 'Unique (Proteins)', ] + pep_cols

    # filter out non-unique, reverse and contaminant peptides
    filtered_pep = peptides[(peptides['Reverse'] != '+') & (
        peptides['Potential contaminant'] != '+')]

    if sample_names is None:
        logger.info(f'Sample names not set. Collecting all samples.')
        logger.debug(f'Columns found: {peptides.columns.tolist()}')
        sample_names = [x.replace('Experiment ', '').split(
            '_')[:-1] for x in peptides.columns.tolist() if 'Experiment ' in x]
        sample_names = list(set([('_').join(x) for x in sample_names]))
        logger.info(f'Samples detected: {sample_names}')

    cleaned_dfs = {}
    standard_cols = ['Sequence', 'Proteins', 'Gene names',
                     'Protein names', 'Unique (Groups)', 'Unique (Proteins)', ] + pep_cols
    for sample in sample_names:
        sample_df = filtered_pep[standard_cols +
                                 [x for x in filtered_pep.columns.tolist() if f' {sample}_' in x]]
        sample_cols = [col for col in sample_df.columns.tolist(
        ) if 'Reporter intensity corrected' in col]
        # filter those without any values for value in variable:
        # In theory (and in previous MQ cleanup), drop any peptides with any missing values here? To start with better to only drop ones that are all missing
        sample_df = sample_df.replace(0, np.nan).dropna(
            axis=0, how='all', subset=sample_cols)
        cleaned_dfs[sample] = sample_df[standard_cols + sample_cols]
    logger.info(f'Successfully cleaned peptide dataframe.')

    return cleaned_dfs


def tmt_proteins(input_path, sample_names=None, prot_cols=[]):

    logger.info(f'Collecting proteins')
    proteins = pd.read_table(input_path, sep='\t')
    logger.info(
        f'Imported proteins from {input_path}. {proteins.shape[0]} entries found.')

    # remove contaminant and reverse proteins
    proteins = proteins[(proteins['Reverse'] != '+') &
                        (proteins['Potential contaminant'] != '+')]
    logger.info(
        f'Removed contaminant and reverse proteins: {proteins.shape[0]} entries remain.')

    cleaned_prot_dfs = {}
    standard_cols = ['Protein IDs', 'Gene names',
                     'Protein names', 'Number of proteins'] + prot_cols

    for sample in sample_names:
        sample_cols = standard_cols + \
            [x for x in proteins.columns.tolist() if f' {sample}_' in x]
        sample_df = proteins[sample_cols]
        ratio_cols = [col for col in sample_df.columns.tolist(
        ) if 'Reporter intensity corrected' in col]
        logger.debug(f'Ratio cols: {ratio_cols}')
        #collect columns of interest
        sample_vals = proteins[sample_cols + ratio_cols]
        #collect only proteins with at least one quantification in that sample
        sample_reps = sample_df.replace(0, np.nan).dropna(
            axis=0, how='all', subset=ratio_cols)
        logger.debug(f'Sample reps: {sample_reps.head(10)}')
        # collect only proteins which are master proteins
        master_proteins = sample_reps
        logger.debug(f'Master proteins: {master_proteins.head(10)}')
        cleaned_prot_dfs[sample] = master_proteins

    logger.info(f'Successfully cleaned proteins dataframe.')

    return cleaned_prot_dfs


def tmt_cleaner(input_folder, output_path, sample_names=None, proteins_file='proteinGroups.txt', peptides_file='peptides.txt', prot_cols=[], pep_cols=[], channel_labels=False):

    cleaned_dfs = {}

    cleaned_peptides = tmt_peptides(
        input_path=f'{input_folder}{peptides_file}', sample_names=sample_names, pep_cols=pep_cols)
    cleaned_proteins = tmt_proteins(
        input_path=f'{input_folder}{proteins_file}', sample_names=sample_names, prot_cols=prot_cols)

    logger.info(f'Sorting cleaned data per sample...')
    ## Collecting specific results for each set of samples for further processing
    for sample in sample_names:
        #collect peptide dataframe, rename relevant columns
        pep_dataframe = cleaned_peptides[sample]
        prot_dataframe = cleaned_proteins[sample]
        prot_dataframe.rename(columns={'Protein IDs': 'Proteins'}, inplace=True)

        if channel_labels == True:

            pep_dataframe = tmt_channel_labels(
                sample_name=sample, dataframe= pep_dataframe,input_path=f"{input_folder.split('/')[0]}/channel_description.xlsx")
            prot_dataframe = tmt_channel_labels(
                sample_name=sample, dataframe=prot_dataframe, input_path=f"{input_folder.split('/')[0]}/channel_description.xlsx")

        #save to individual excel spreadsheets
        FileHandling.df_to_excel(f'{output_folder}{sample}_Compiled.xlsx', sheetnames=['Proteins', 'Peptides'], data_frames=[prot_dataframe, pep_dataframe])

    logger.info(
        f'Proteins and peptides successfully cleaned. Dataframes save to {output_folder}.')

    return cleaned_peptides, cleaned_proteins


def tmt_channel_labels(sample_name, dataframe, input_path=f'experimental_data/channel_description.xlsx'):

    """Melts quantitative information from either peptides or proteins dataframe, then assigns sample label derived from template file located at input_path.
    Template file should have standard layout, where sheetname matches sample_name and each sheet contains a dataframe of labels (rows) vs plexes (columns) in which each cell then contains the corresponding sample label. An example file can be located here: https://github.com/dezeraecox-experiments/EXP121_Proteome-abundance-Dnmt3-stemcells/blob/master/raw_data/channel_description.xlsx.
    Resultant dataframe is then unmelted to produce wideform table where each quantitative column is labelled by the sample name (as opposed to channel_replicate)."""

    # read in channel descriptions, format for labels
    channel_layout = pd.read_excel(input_path, sheet_name=f'{sample_name}')
    channel_layout = pd.melt(channel_layout,
                            id_vars=['Label', 'Channel'],
                            value_vars=[col for col in channel_layout.columns.tolist() if 'Plex' in col],
                            var_name='Plex',
                            value_name='sample_name',
                            )
    channel_layout = channel_layout[channel_layout['sample_name'] != 'Empty']
    channel_layout['Plex'] = channel_layout['Plex'].str.strip('Plex ')
    channel_layout['channel_plex'] = [f'{channel}_{plex}' for channel, plex in channel_layout[['Channel', 'Plex']].values]
    channel_layout = dict(channel_layout[['channel_plex', 'sample_name']].values)


    # Label TMT channels using channel description file
    info_cols = [col for col in dataframe.columns.tolist(
    ) if col in ['Protein IDs', 'Sequence', 'Proteins']]
    sample_cols = [col for col in dataframe.columns.tolist(
    ) if 'Reporter intensity corrected ' in col]

    # Assign channels to sample-name_replicates
    quant_info = dataframe[info_cols + sample_cols].copy()
    quant_info = pd.melt(quant_info,
                         id_vars=info_cols,
                         value_vars=sample_cols,
                         var_name='channel_plex',
                         value_name='abundance',
                         )
    quant_info['channel'] = quant_info['channel_plex'].str.strip(
        'Reporter intensity corrected ').str.split(' ').str[0]
    quant_info['plex'] = quant_info['channel_plex'].str.split('_').str[1]
    quant_info['channel_plex'] = quant_info['channel'] + \
        '_' + quant_info['plex']
    quant_info['sample_name'] = quant_info['channel_plex'].map(channel_layout)
    quant_info.dropna(subset=['sample_name'], inplace=True)

    return quant_info



if __name__ == "__main__":

    input_folder = 'data/'
    output_folder = 'python_results/initial_cleanup/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    tmt_cleaner(
        input_folder, 
        output_folder,
        sample_names=['AGG', 'SUP'],
        proteins_file='proteinGroups.txt',
        peptides_file='peptides.txt',
        prot_cols=[],
        pep_cols=[],
        channel_labels=True
        )

