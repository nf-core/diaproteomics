#!/usr/bin/env python
import argparse
import pandas as pd


"""
reformat_output_for_msstats.py: 
Reformat the peptide quantity output of the DIAlignR chromatogram alignment to pass it onto msstats
"""

__author__      = "Leon Bichmann"


def main():
    model = argparse.ArgumentParser(description='Reformat dialignR output for MSstats')

    model.add_argument(
        '-i', '--input',
        type=str,
        help='dialignR file'
    )

    model.add_argument(
        '-e', '--exp_design',
        type=str,
        help='experimental design file'
    )

    model.add_argument(
        '-l', '--library',
        type=str,
        help='library file'
    )

    model.add_argument(
        '-f', '--fdr_level',
        type=str,
        help='fdr level'
    )

    model.add_argument(
        '-o', '--output',
        type=str,
        help='output reformatted file'
    )

    args = model.parse_args()

    exp_design=args.exp_design
    dialignR=args.input
    lib=args.library
    fdr_level=args.fdr_level

    #parse library and rename columns
    df_lib=pd.read_csv(lib,sep='\t')
    df_lib=df_lib[['ModifiedPeptideSequence', 'ProteinId']]
    df_lib.columns=['PeptideSequence', 'ProteinName']

    if fdr_level=='peptide':
       df_lib['ProteinName']=df_lib['PeptideSequence']

    #parse experimental design and rename columns
    df_I=pd.read_csv(exp_design, sep='\t')
    if any(['/' in f for f in df_I['Spectra_Filepath'].values.tolist()]):
       files=[f.replace('.mzML','').replace('.raw','').replace('.Raw','').replace('.RAW','').split('/')[-1] for f in df_I['Spectra_Filepath'].values.tolist()]
    else:
       files=[f.replace('.mzML','').replace('.raw','').replace('.Raw','').replace('.RAW','') for f in df_I['Spectra_Filepath'].values.tolist()]
    print(files)
    df_I['run']=files
    df_I=df_I[['Sample', 'Sample_Group', 'MSstats_Condition', 'run']]
    df_I.colums=['Sample', 'Sample_Group','MSstats_Condition', 'run']

    #parse dialignr output and merge with experimental design and library
    df=pd.read_csv(dialignR, sep=',')
    df=df.merge(df_I, how='outer', on='run')

    df=df[['sequence', 'charge', 'MSstats_Condition', 'Sample_Group', 'Sample', 'intensity']]
    df.columns=['PeptideSequence', 'PrecursorCharge', 'Condition', 'BioReplicate', 'Run', 'Intensity']
    df['IsotopeLabelType']='L'
    df['ProductCharge']=0
    df['FragmentIon']='NA'

    df=df.merge(df_lib, how='inner', on='PeptideSequence')
    df=df.drop_duplicates()

    #output reformatted
    df.to_csv(args.output, sep=',', index=False)


if __name__=="__main__":
   main()
