#!/usr/bin/env python
from __future__ import print_function
import sys
import scipy
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import glob
import argparse


"""
select_pseudo_irts_from_lib.py: 
This script selects and exports a specified number of high-intensity peptides spanning 
the entire RT range to serve as iRT standards for the RT alignment between library and DIA MS run.
"""

__author__      = "Leon Bichmann"


# Select iRT standards based on the library intensity spanning the RT range
def get_pseudo_irts(lib, n_irts, min_rt, max_rt):

    # select irts based on dda Intensity
    df_pre = pd.read_csv(lib, sep='\t')
    # sum up individual transition intensities
    df_sum = df_pre.groupby(['ModifiedPeptideSequence', 'PrecursorCharge'])['LibraryIntensity'].apply(sum).reset_index()
    df_merged = df_pre.merge(df_sum, on=['ModifiedPeptideSequence', 'PrecursorCharge'])[['ModifiedPeptideSequence', 'NormalizedRetentionTime', 'LibraryIntensity_y']]

    rt_sample_space = np.linspace(min_rt, max_rt, n_irts)
    rt_sub_df = []
    for rt in rt_sample_space:
        try:
            best_pep = df_merged[(df_pre['NormalizedRetentionTime'] < rt + 0.5) & (df_merged['NormalizedRetentionTime'] > rt)].sort_values('LibraryIntensity_y', ascending=False).iloc[0]['ModifiedPeptideSequence']
            rt_sub_df.append(best_pep)
        except:
            pass

    irts = list(set(rt_sub_df))

    df_sub = df_pre[(df_pre['ModifiedPeptideSequence'].isin(irts))]

    return df_sub


def main():
    model = argparse.ArgumentParser(description='Postprocess Neoepitopes predicted by MHCNuggets')

    model.add_argument(
        '-i', '--input_libraries',
        type=str,
        help='library file serving as reference'
    )

    model.add_argument(
        '-n', '--n_irts',
        type=int,
        help='number of irts to select for alignment'
    )

    model.add_argument(
        '-ri', '--min_rt',
        type=int,
        help='minimum rt of irts to select for alignment'
    )

    model.add_argument(
        '-rn', '--max_rt',
        type=int,
        help='maximum rt of irts to select for alignment'
    )

    model.add_argument(
        '-o', '--output',
        type=str,
        help='output aligned library file'
    )

    args = model.parse_args()
    print(args)

    lib=args.input_libraries
    n_irts=args.n_irts
    min_rt=args.min_rt
    max_rt=args.max_rt

    df_sub=get_pseudo_irts(lib, n_irts, min_rt, max_rt)
    df_sub.to_csv(args.output, sep='\t')


if __name__=="__main__":
   main()
