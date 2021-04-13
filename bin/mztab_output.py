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
    df_lib=df_lib[['ModifiedPeptideSequence', 'PrecursorMz', 'ProteinId']]
    df_lib.columns=['PeptideSequence', 'PrecursorMz', 'ProteinName']

    #parse experimental design and rename columns
    df_I=pd.read_csv(exp_design, sep='\t')
    if any(['/' in f for f in df_I['Spectra_Filepath'].values.tolist()]):
       files=[f.replace('.mzML','').replace('.raw','').replace('.Raw','').replace('.RAW','').split('/')[-1] for f in df_I['Spectra_Filepath'].values.tolist()]
    else:
       files=[f.replace('.mzML','').replace('.raw','').replace('.Raw','').replace('.RAW','') for f in df_I['Spectra_Filepath'].values.tolist()]
    print(files)
    df_I['run']=files
    df_I=df_I[['Sample', 'BatchID', 'MSstats_Condition', 'run']]
    df_I.columns=['Sample', 'BatchID','MSstats_Condition', 'run']

    #parse dialignr output and merge with experimental design and library
    df=pd.read_csv(dialignR, sep='\t')
    df=df.merge(df_I, how='outer', on='run')

    df=df[['sequence', 'RT', 'charge', 'Sample', 'intensity']]
    df.columns=['PeptideSequence', 'retention_time', 'PrecursorCharge', 'Run', 'Intensity']

    df=df.merge(df_lib, how='inner', on='PeptideSequence')
    df=df.drop_duplicates()

    #reformat MTD header
    header='\t'.join(['MTD','mzTab-version','1.0'])+'\n\n'

    #reformat as PRT protein mzTab
    remaining_prh_cols='unit_id,description,taxid,species,database,database_version,search_engine,search_engine_score,reliability,num_peptides_unambiguous,ambiguity_members,modifications,uri,go_terms,protein_coverage'.split(',')

    num_peptides=df.groupby(['ProteinName'])['PeptideSequence'].apply(list).apply(len).reset_index()['PeptideSequence']
    num_peptides_distinct=df.groupby(['ProteinName'])['PeptideSequence'].apply(list).apply(set).apply(len).reset_index()['PeptideSequence']
    df_prot=df.groupby(['ProteinName','Run'])['Intensity'].apply(sum).reset_index().pivot(index='ProteinName', columns='Run', values='Intensity').reset_index()
    col_idxs=df_prot.columns[1:]
    df_prot.columns=['accession']+['protein_abundance_sub['+str(i)+']' for i in col_idxs]
    for i in col_idxs:
        df_prot['protein_abundance_stdev_sub['+str(i)+']']='null'
        df_prot['protein_abundance_std_error_sub['+str(i)+']']='null'

    cols=df_prot.columns.tolist()
    df_prot['PRH']='PRT'
    df_prot=df_prot[['PRH']+cols]
    df_prot['num_peptides']=num_peptides
    df_prot['num_peptides_distinct']=num_peptides_distinct
    for c in remaining_prh_cols:
        df_prot[c]='null'

    #reformat as PEP peptide mzTab
    remaining_pep_cols='unit_id,unique,database,database_version,search_engine_score,reliability,modifications,uri,spectra_ref'.split(',')

    retention_time=df.groupby(['PeptideSequence','ProteinName','PrecursorCharge','PrecursorMz'])['retention_time'].apply(np.median).reset_index()['retention_time']
    df_pep=df.groupby(['PeptideSequence','ProteinName','PrecursorCharge','PrecursorMz','Run'])['Intensity'].apply(sum).unstack('Run').reset_index()
    col_idxs=df_pep.columns[4:]
    df_pep.columns=['sequence','accession','charge','mass_to_charge']+['peptide_abundance_sub['+str(i)+']' for i in col_idxs]
    for i in col_idxs:
        df_pep['peptide_abundance_stdev_sub['+str(i)+']']='null'
        df_pep['peptide_abundance_std_error_sub['+str(i)+']']='null'

    cols=df_pep.columns.tolist()
    df_pep['PEH']='PEP'
    df_pep=df_pep[['PEH']+cols]
    df_pep['retention_time']=retention_time
    df_pep['search_engine']='OpenSwathWorkflow'
    for c in remaining_pep_cols:
        df_pep[c]='null'

    #output mzTab
    

if __name__=="__main__":
    main()
