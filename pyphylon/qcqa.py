"""
Functions for quality control and quality assurance of genomic data.
"""

import os
import numpy as np
import pandas as pd
from kneebow.rotor import Rotor


# Main filtration workflow methods

def filter_by_species(summary, species_name):
    """
    Filter the summary DataFrame for a specific species.
    
    Parameters:
    - summary (pd.DataFrame): DataFrame containing genome summary data.
    - species_name (str): The name of the species to filter by.

    Returns:
    - pd.DataFrame: Filtered DataFrame containing only the specified species.
    """
    species_summary = summary[summary["genome_name"].str.contains(species_name)]
    species_summary = species_summary.dropna(subset=['genome_length'])
    species_summary = species_summary.dropna(subset=['patric_cds'])
    
    # Ensure genome_length and patric_cds are ints
    species_summary['genome_length'] = species_summary['genome_length'].astype('int')
    species_summary['patric_cds'] = species_summary['patric_cds'].astype('int')

    return species_summary


def filter_by_genome_quality(
    species_summary,
    min_thresh_n50=None,
    max_contig=None,
    contamination_cutoff=None,
    completeness_cutoff=None,
    checkm_filter_statuses=('WGS',),
    checkm_missing='keep',
    return_stats=True,
):
    """
    Filter the genome data based on various quality metrics.

    Parameters:
    - species_summary (pd.DataFrame): DataFrame containing genome data.
    - min_thresh_n50 (int, optional): Minimum threshold for N50 score.
    - max_contig (int, optional): Maximum allowed number of contigs.
    - contamination_cutoff (float, optional): Cutoff for CheckM contamination score.
    - completeness_cutoff (float, optional): Cutoff for CheckM completeness score.
    - checkm_filter_statuses (tuple, optional): Genome statuses to apply CheckM
        filtering to. Default is ('WGS',) for backward compatibility. Set to
        ('Complete', 'WGS') to filter all genome types, or None to filter all.
    - checkm_missing (str, optional): How to handle genomes with null CheckM data.
        'keep' retains them (default), 'drop' discards them.
    - return_stats (bool, optional): Whether to return filtration statistics.

    Returns:
    - pd.DataFrame: Filtered DataFrame containing high-quality genomes.
    - pd.DataFrame (optional): Filtration statistics if return_stats=True.
    """

    species_complete_summary = species_summary[species_summary.genome_status == 'Complete']
    species_wgs_summary = species_summary[species_summary.genome_status == 'WGS']

    # Record initial lengths (for filtration statistics)
    filtration_metrics_list = ['prefiltration', 'L50/N50', 'contig_count', 'CheckM_completeness_contamination']
    filtration_columns = ['initial', 'num_filtered', 'remaining']
    df_filtration = pd.DataFrame(index=filtration_metrics_list, columns=filtration_columns)

    df_filtration.loc['prefiltration', 'initial'] = species_summary.shape[0]
    df_filtration.loc['prefiltration', 'remaining'] = df_filtration.loc['prefiltration', 'initial']
    df_filtration.loc['prefiltration', 'num_filtered'] = 0

    # Filter complete sequences by L50 & N50 score metrics
    species_complete_summary = _filter_l50(species_complete_summary)
    species_complete_summary = _filter_n50(species_complete_summary, min_thresh_n50)

    # Record L50 & N50 filtration metrics
    df_filtration.loc['L50/N50', 'initial'] = df_filtration.loc['prefiltration', 'remaining']
    df_filtration.loc['L50/N50', 'remaining'] = species_complete_summary.shape[0] + species_wgs_summary.shape[0]
    df_filtration.loc['L50/N50', 'num_filtered'] = df_filtration.loc['L50/N50', 'initial'] - df_filtration.loc['L50/N50', 'remaining']

    # Filter other WGS sequences by contig count
    species_wgs_summary = _filter_by_contig(species_wgs_summary, max_contig)

    # Record contig count filtration metrics
    df_filtration.loc['contig_count', 'initial'] = df_filtration.loc['L50/N50', 'remaining']
    df_filtration.loc['contig_count', 'remaining'] = species_complete_summary.shape[0] + species_wgs_summary.shape[0]
    df_filtration.loc['contig_count', 'num_filtered'] = df_filtration.loc['contig_count', 'initial'] - df_filtration.loc['contig_count', 'remaining']

    # Merge Complete + WGS, then apply CheckM filtering to selected statuses
    filtered = pd.concat([species_complete_summary, species_wgs_summary])

    if checkm_filter_statuses is None:
        checkm_subset = filtered
        skip_subset = filtered.iloc[0:0]  # empty DataFrame with same columns
    else:
        checkm_mask = filtered['genome_status'].isin(checkm_filter_statuses)
        checkm_subset = filtered[checkm_mask]
        skip_subset = filtered[~checkm_mask]

    checkm_subset = _filter_checkM_contamination(checkm_subset, contamination_cutoff, checkm_missing)
    checkm_subset = _filter_checkM_completeness(checkm_subset, completeness_cutoff, checkm_missing)

    filtered_species_summary = pd.concat([checkm_subset, skip_subset])

    df_filtration.loc['CheckM_completeness_contamination', 'initial'] = df_filtration.loc['contig_count', 'remaining']
    df_filtration.loc['CheckM_completeness_contamination', 'remaining'] = filtered_species_summary.shape[0]
    df_filtration.loc['CheckM_completeness_contamination', 'num_filtered'] = df_filtration.loc['CheckM_completeness_contamination', 'initial'] - df_filtration.loc['CheckM_completeness_contamination', 'remaining']

    # Typecast relevant columns as numeric
    filtered_species_summary['contig_l50'] = filtered_species_summary['contig_l50'].astype('int')
    filtered_species_summary['contig_n50'] = filtered_species_summary['contig_n50'].astype('int')
    filtered_species_summary['contigs'] = filtered_species_summary['contigs'].astype('int')
    filtered_species_summary['checkm_contamination'] = pd.to_numeric(filtered_species_summary['checkm_contamination'], errors='coerce')
    filtered_species_summary['checkm_completeness'] = pd.to_numeric(filtered_species_summary['checkm_completeness'], errors='coerce')
    filtered_species_summary['gc_content'] = filtered_species_summary['gc_content'].astype('float')
    
    if return_stats:
        return filtered_species_summary, df_filtration
    else:
        return filtered_species_summary


# Individual filtration functions

def _filter_l50(species_complete_summary, l50_score=1):
    """
    Filter genomes by L50 score.

    Parameters:
    - species_complete_summary (pd.DataFrame): DataFrame containing complete genomes.
    - l50_score (int, optional): L50 score to filter by. Default is 1.

    Returns:
    - pd.DataFrame: Filtered DataFrame with genomes having the specified L50 score.
    """
    species_complete_summary = species_complete_summary.dropna(subset=['contig_l50'])
    species_complete_summary['contig_l50'] = species_complete_summary['contig_l50'].astype('int')
    
    good_l50 = species_complete_summary['contig_l50'] == l50_score
    species_complete_summary = species_complete_summary[good_l50]

    return species_complete_summary


def _filter_n50(species_complete_summary, min_thresh_n50):
    """
    Filter genomes by N50 score.

    Parameters:
    - species_complete_summary (pd.DataFrame): DataFrame containing complete genomes.
    - min_thresh_n50 (int): Minimum N50 score threshold.

    Returns:
    - pd.DataFrame: Filtered DataFrame with genomes having N50 scores above the threshold.
    """
    species_complete_summary = species_complete_summary.dropna(subset=['contig_n50'])
    species_complete_summary['contig_n50'] = species_complete_summary['contig_n50'].astype('int')
    
    if min_thresh_n50:
        cond = species_complete_summary['contig_n50'] > min_thresh_n50
        species_complete_summary = species_complete_summary[cond]
    
    return species_complete_summary


def _filter_by_contig(species_wgs_summary, max_contig):
    """
    Filter genomes by contig count.

    Parameters:
    - species_wgs_summary (pd.DataFrame): DataFrame containing WGS genomes.
    - max_contig (int, optional): Maximum allowed number of contigs.

    Returns:
    - pd.DataFrame: Filtered DataFrame with genomes having contig counts below the threshold.
    """
    if species_wgs_summary.empty:
        return species_wgs_summary
    species_wgs_summary = species_wgs_summary.dropna(subset=['contigs'])
    species_wgs_summary['contigs'] = species_wgs_summary['contigs'].astype('int')
    
    if max_contig:
        species_wgs_summary = species_wgs_summary[species_wgs_summary['contigs'] <= max_contig]
    else:
        species_wgs_summary = _remove_contig_outliers(species_wgs_summary)

    return species_wgs_summary


def _filter_checkM_contamination(species_wgs_summary, contamination_cutoff, checkm_missing='keep'):
    """
    Filter genomes by CheckM contamination score.

    Parameters:
    - species_wgs_summary (pd.DataFrame): DataFrame containing WGS genomes.
    - contamination_cutoff (float, optional): Maximum allowed CheckM contamination score.
    - checkm_missing (str, optional): How to handle genomes with null CheckM data.
        'keep' retains them (default), 'drop' discards them.

    Returns:
    - pd.DataFrame: Filtered DataFrame with genomes having contamination scores below the cutoff.
    """
    if species_wgs_summary.empty:
        return species_wgs_summary

    has_data = species_wgs_summary.dropna(subset=['checkm_contamination']).copy()
    missing_data = species_wgs_summary[species_wgs_summary['checkm_contamination'].isna()]

    has_data['checkm_contamination'] = has_data['checkm_contamination'].astype('float')

    if contamination_cutoff:
        cond = has_data['checkm_contamination'] < contamination_cutoff
    else:
        kneebow_cutoff = _get_kneebow_cutoff(has_data, column='checkm_contamination', curve='elbow')
        cond = has_data['checkm_contamination'] < kneebow_cutoff

    has_data = has_data[cond]

    if checkm_missing == 'keep':
        return pd.concat([has_data, missing_data])
    else:
        return has_data


def _filter_checkM_completeness(species_wgs_summary, completeness_cutoff, checkm_missing='keep'):
    """
    Filter genomes by CheckM completeness score.

    Parameters:
    - species_wgs_summary (pd.DataFrame): DataFrame containing WGS genomes.
    - completeness_cutoff (float, optional): Minimum allowed CheckM completeness score.
    - checkm_missing (str, optional): How to handle genomes with null CheckM data.
        'keep' retains them (default), 'drop' discards them.

    Returns:
    - pd.DataFrame: Filtered DataFrame with genomes having completeness scores above the cutoff.
    """
    if species_wgs_summary.empty:
        return species_wgs_summary

    has_data = species_wgs_summary.dropna(subset=['checkm_completeness']).copy()
    missing_data = species_wgs_summary[species_wgs_summary['checkm_completeness'].isna()]

    has_data['checkm_completeness'] = has_data['checkm_completeness'].astype('float')

    if completeness_cutoff:
        cond = has_data['checkm_completeness'] > completeness_cutoff
    else:
        kneebow_cutoff = _get_kneebow_cutoff(has_data, column='checkm_completeness', curve='knee')
        cond = has_data['checkm_completeness'] > kneebow_cutoff

    has_data = has_data[cond]

    if checkm_missing == 'keep':
        return pd.concat([has_data, missing_data])
    else:
        return has_data


# Helper functions

def _remove_contig_outliers(species_wgs_summary):
    """
    Remove outliers in contig counts based on IQR and median thresholds.

    Parameters:
    - species_wgs_summary (pd.DataFrame): DataFrame containing WGS genomes.

    Returns:
    - pd.DataFrame: Filtered DataFrame with outliers removed.
    """
    df = species_wgs_summary.copy()
    
    # Step 1: Calculate the IQR
    Q1 = df['contigs'].quantile(0.25)
    Q3 = df['contigs'].quantile(0.75)
    IQR = Q3 - Q1
    
    # Step 2: Find the upper fence
    upper_fence = Q3 + 1.5 * IQR

    # Step 3: Remove all entries above upper fence limit (clear outliers)
    outlier_cond = species_wgs_summary['contigs'] <= upper_fence
    species_wgs_summary = species_wgs_summary[outlier_cond]
    
    # Step 4: Remove all entries > 2.5 * median of remaining values
    upper_limit_median = 2.5 * species_wgs_summary['contigs'].median()
    outlier_cond2 = species_wgs_summary['contigs'] <= upper_limit_median
    species_wgs_summary = species_wgs_summary[outlier_cond2]
    
    return species_wgs_summary


def _get_kneebow_cutoff(species_wgs_summary, column, curve):
    """
    Determine the cutoff value using the Kneebow package.

    Parameters:
    - species_wgs_summary (pd.DataFrame): DataFrame containing WGS genomes.
    - column (str): Column name to determine the cutoff for.
    - curve (str): Type of curve ('elbow' or 'knee').

    Returns:
    - float: The cutoff value.
    """
    if species_wgs_summary.empty:
        return 0.0

    df = species_wgs_summary.copy()

    if curve.lower() == 'elbow':
        df = df.sort_values(by=column, ascending=True)
    elif curve.lower() == 'knee':
        df = df.sort_values(by=column, ascending=False)
    else:
        raise ValueError(f'curve must be either "elbow" or "knee". {curve} was provided instead.')
    
    df = df.reset_index()
    
    results_itr = zip(list(df.index), list(df[column]))
    data = list(results_itr)
    
    rotor = Rotor()
    rotor.fit_rotate(data)
    
    if curve.lower() == 'elbow':
        elbow_idx = rotor.get_elbow_index()
        kneebow_cutoff = df[column][elbow_idx]
    elif curve.lower() == 'knee':
        knee_idx = rotor.get_knee_index()
        kneebow_cutoff = df[column][knee_idx]
    
    return kneebow_cutoff


def append_entry(df_filtration, entry_name, entry_remaining):
    """
    Append a new entry to the filtration statistics DataFrame.

    Parameters:
    - df_filtration (pd.DataFrame): DataFrame containing filtration statistics.
    - entry_name (str): Name of the new entry.
    - entry_remaining (int): Number of remaining genomes after filtration.

    Returns:
    - pd.DataFrame: Updated filtration statistics DataFrame.
    """
    df_temp = pd.DataFrame(index=[entry_name], columns=df_filtration.columns)
    
    df_temp.loc[entry_name, 'initial'] = df_filtration.iloc[-1]['remaining']
    df_temp.loc[entry_name, 'remaining'] = entry_remaining
    df_temp.loc[entry_name, 'num_filtered'] = df_temp.loc[entry_name, 'initial'] - df_temp.loc[entry_name, 'remaining']
    
    df_filtration = pd.concat([df_filtration, df_temp])

    return df_filtration
