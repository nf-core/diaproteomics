#!/usr/bin/env python
import argparse
import pandas as pd

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
        '-o', '--output',
        type=str,
        help='output reformatted file'
    )

    args = model.parse_args()

    exp_design=args.exp_design
    dialignR=args.input
    lib=args.library

    #parse library and rename columns
    df_lib=pd.read_csv(lib,sep='\t')
    df_lib=df_lib[['ModifiedPeptideSequence', 'ProteinId']]
    df_lib.columns=['PeptideSequence', 'ProteinName']

    #parse experimental design and rename columns
    df_I=pd.read_csv(exp_design, sep='\t')
    files=[f.replace('.mzML','').split('/')[-1] for f in df_I['Spectra_Filepath'].values.tolist()]
    df_I['run']=files
    df_I=df_I[['Fraction_Group', 'Sample', 'Fraction', 'run']]
    df_I.colums=['Fraction_Group', 'Sample', 'Fraction', 'run']

    #parse dialignr output and merge with experimental design and library
    df=pd.read_csv(dialignR, sep=',')
    df=df.merge(df_I, how='outer', on='run')

    df=df[['sequence', 'charge', 'Fraction', 'Fraction_Group', 'intensity', 'run']]
    df.columns=['PeptideSequence', 'PrecursorCharge', 'Condition', 'Run', 'Intensity', 'Reference']
    df['BioReplicate']=df['Run']
    df['IsotopeLabelType']='L'
    df['ProductCharge']=0
    df['FragmentIon']='NA'
    files=[f+'.mzML' for f in df['Reference'].values.tolist()]
    df['Reference']=files

    df=df.merge(df_lib, how='inner', on='PeptideSequence')
    df=df.drop_duplicates()

    #output reformatted
    df.to_csv(args.output,sep='\t', index=False)


if __name__=="__main__":
   main()
