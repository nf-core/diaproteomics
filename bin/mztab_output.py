#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np


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
        '-t', '--fdr_threshold_pep',
        type=float,
        help='fdr treshold'
    )

    model.add_argument(
        '-p', '--fdr_threshold_prot',
        type=float,
        help='fdr treshold'
    )

    model.add_argument(
        '-ms1', '--ms1_scoring',
        type=str,
        help='ms1_scoring'
    )

    model.add_argument(
        '-r', '--rt_extraction_window',
        type=float,
        help='rt_extraction_window'
    )

    model.add_argument(
        '-m', '--mz_extraction_window',
        type=float,
        help='mz_extraction_window'
    )

    model.add_argument(
        '-m1','--mz_extraction_window_ms1',
        type=float,
        help='mz_extraction_window_ms1'
    )

    model.add_argument(
        '-mu', '--mz_extraction_unit',
        type=str,
        help='mz_extraction_unit'
    )

    model.add_argument(
        '-mu1', '--mz_extraction_unit_ms1',
        type=str,
        help='mz_extraction_unit_ms1'
    )

    model.add_argument(
        '-dg', '--dialignr_global_align_fdr',
        type=float,
        help='dialignr_global_align_fdr'
    )

    model.add_argument(
        '-dn','--dialignr_align_fdr',
        type=float,
        help='dialignr_align_fdr'
    )

    model.add_argument(
        '-du', '--dialignr_unalign_fdr',
        type=float,
        help='dialignr_unalign_fdr'
    )

    model.add_argument(
        '-da', '--dialignr_analyte_fdr',
        type=float,
        help='dialignr_analyte_fdr'
    )

    model.add_argument(
        '-dq', '--dialignr_query_fdr',
        type=float,
        help='dialignr_query_fdr'
    )

    model.add_argument(
        '-w', '--workflow_version',
        type=str,
        help='workflow version'
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
    if fdr_level=='protein':
       fdr_short='prot'
       fdr_threshold=str(args.fdr_threshold_prot)
    else:
       fdr_short='pep'
       fdr_threshold=str(args.fdr_threshold_pep) 

    ms1_scoring=str(args.ms1_scoring)
    rt_extraction_window=str(args.rt_extraction_window)
    mz_extraction_window=str(args.mz_extraction_window)
    mz_extraction_window_ms1=str(args.mz_extraction_window_ms1)
    mz_extraction_unit=str(args.mz_extraction_unit)
    mz_extraction_unit_ms1=str(args.mz_extraction_unit_ms1)
    dialignr_global_align_fdr=str(args.dialignr_global_align_fdr)
    dialignr_align_fdr=str(args.dialignr_align_fdr)
    dialignr_unalign_fdr=str(args.dialignr_unalign_fdr)
    dialignr_analyte_fdr=str(args.dialignr_analyte_fdr)
    dialignr_query_fdr=str(args.dialignr_query_fdr)
    diaproteomics_version=str(args.workflow_version)

    #parse library and rename columns
    df_lib=pd.read_csv(lib,sep='\t')
    df_lib=df_lib[['ModifiedPeptideSequence', 'PrecursorMz', 'ProteinId']]
    df_lib.columns=['PeptideSequence', 'PrecursorMz', 'ProteinName']

    #parse experimental design and rename columns
    df_I=pd.read_csv(exp_design, sep='\t')
    if any(['/' in f for f in df_I['Spectra_Filepath'].values.tolist()]):
       files=[f.replace('.mzml','').replace('.mzML','').replace('.raw','').replace('.Raw','').replace('.RAW','').split('/')[-1] for f in df_I['Spectra_Filepath'].values.tolist()]
    else:
       files=[f.replace('.mzml','').replace('.mzML','').replace('.raw','').replace('.Raw','').replace('.RAW','') for f in df_I['Spectra_Filepath'].values.tolist()]
    print(files)
    df_I['run']=files
    df_I=df_I[['Sample', 'BatchID', 'MSstats_Condition', 'run']]
    df_I.columns=['Sample', 'BatchID','MSstats_Condition', 'run']

    assays=list(set(df_I['MSstats_Condition'].tolist()))
    assayDict={}
    for a in enumerate(assays):
        assayDict[a[1]]=a[0]+1

    #parse dialignr output and merge with experimental design and library
    df=pd.read_csv(dialignR, sep='\t')
    df=df.merge(df_I, how='outer', on='run')

    df=df[['sequence', 'RT', 'charge', 'Sample', 'intensity']]
    df.columns=['PeptideSequence', 'retention_time', 'PrecursorCharge', 'Run', 'Intensity']

    df=df.merge(df_lib, how='inner', on='PeptideSequence')
    df=df.drop_duplicates()

    #reformat as PRT protein mzTab
    remaining_prh_cols='description,taxid,species,database,database_version,best_search_engine_score,ambiguity_members,modifications'.split(',')

    num_peptides=df.groupby(['ProteinName'])['PeptideSequence'].apply(list).apply(len).reset_index()['PeptideSequence']
    num_peptides_distinct=df.groupby(['ProteinName'])['PeptideSequence'].apply(list).apply(set).apply(len).reset_index()['PeptideSequence']
    df_prot=df.groupby(['ProteinName','Run'])['Intensity'].apply(sum).reset_index().pivot(index='ProteinName', columns='Run', values='Intensity').reset_index()
    col_idxs=df_prot.columns[1:].tolist()
    col_idxs_enum=[i[0]+1 for i in enumerate(col_idxs)]
    df_prot.columns=['accession']+['protein_abundance_study_variable['+str(i)+']' for i in col_idxs_enum]
    for i in col_idxs_enum:
        df_prot['protein_abundance_stdev_study_variable['+str(i)+']']='null'
        df_prot['protein_abundance_std_error_study_variable['+str(i)+']']='null'

    cols=df_prot.columns.tolist()
    df_prot['PRH']='PRT'
    df_prot=df_prot[['PRH']+cols]
    df_prot['search_engine']='OpenSwathWorkflow'
    df_prot['num_peptides']=num_peptides
    df_prot['num_peptides_distinct']=num_peptides_distinct
    for c in remaining_prh_cols:
        df_prot[c]='null'

    df_prot.to_csv('PRT_dataframe.csv',index=False, sep='\t')
    op=open('PRT_dataframe.csv','r')
    opr_prt=op.readlines()
    op.close()

    #reformat as PEP peptide mzTab
    remaining_pep_cols='retention_time_window,unique,database,database_version,search_engine_score,modifications'.split(',')

    retention_time=df.groupby(['PeptideSequence','ProteinName','PrecursorCharge','PrecursorMz'])['retention_time'].apply(np.median).reset_index()['retention_time']
    df_pep=df.groupby(['PeptideSequence','ProteinName','PrecursorCharge','PrecursorMz','Run'])['Intensity'].apply(sum).unstack('Run').reset_index()
    df_pep.columns=['sequence','accession','charge','mass_to_charge']+['peptide_abundance_study_variable['+str(i)+']' for i in col_idxs_enum]
    for i in col_idxs_enum:
        df_pep['peptide_abundance_stdev_study_variable['+str(i)+']']='null'
        df_pep['peptide_abundance_std_error_study_variable['+str(i)+']']='null'

    cols=df_pep.columns.tolist()
    df_pep['PEH']='PEP'
    df_pep=df_pep[['PEH']+cols]
    df_pep['retention_time']=retention_time
    df_pep['search_engine']='OpenSwathWorkflow'
    df_pep['charge']=df_pep['charge'].apply(int)
    for c in remaining_pep_cols:
        df_pep[c]='null'

    df_pep.to_csv('PEP_dataframe.csv',index=False, sep='\t')

    op=open('PEP_dataframe.csv','r')
    opr_pep=op.readlines()
    op.close()

    #reformat MTD header
    header=['\t'.join(['MTD','mzTab-version','1.0'])+'\n']
    header.append('\t'.join(['MTD','mzTab-mode','Summary'])+'\n')
    header.append('\t'.join(['MTD','mzTab-type','Quantification'])+'\n')
    header.append('\t'.join(['MTD','description','mztab-like output summarizing DIAproteomics search results'])+'\n')
    header.append('\t'.join(['MTD','software','nfcore/DIAproteomics V.'+diaproteomics_version])+'\n')
    header.append('\t'.join(['MTD', 'software[1]-setting[1]', 'ms1_scoring='+ms1_scoring])+'\n')
    header.append('\t'.join(['MTD', 'software[1]-setting[2]', 'rt_extraction_window='+rt_extraction_window])+'\n')
    header.append('\t'.join(['MTD', 'software[1]-setting[3]', 'mz_extraction_window='+mz_extraction_window])+'\n')
    header.append('\t'.join(['MTD', 'software[1]-setting[4]', 'mz_extraction_window_ms1='+mz_extraction_window_ms1])+'\n')
    header.append('\t'.join(['MTD', 'software[1]-setting[5]', 'mz_extraction_unit='+mz_extraction_unit])+'\n')
    header.append('\t'.join(['MTD', 'software[1]-setting[6]', 'mz_extraction_unit_ms1='+mz_extraction_unit_ms1])+'\n')
    header.append('\t'.join(['MTD', 'software[1]-setting[7]', 'dialignr_global_align_fdr='+dialignr_global_align_fdr])+'\n')
    header.append('\t'.join(['MTD', 'software[1]-setting[8]', 'dialignr_align_fdr='+dialignr_align_fdr])+'\n')
    header.append('\t'.join(['MTD', 'software[1]-setting[9]', 'dialignr_unalign_fdr='+dialignr_unalign_fdr])+'\n')
    header.append('\t'.join(['MTD', 'software[1]-setting[10]', 'dialignr_analyte_fdr='+dialignr_analyte_fdr])+'\n')
    header.append('\t'.join(['MTD', 'software[1]-setting[11]', 'dialignr_query_fdr='+dialignr_query_fdr])+'\n')
    header.append('\t'.join(['MTD','false_discovery_rate', fdr_short+':global FDR', fdr_threshold])+'\n')
    header.append('\t'.join(['MTD','protein_search_engine_score', 'pyProphet q-value'])+'\n')
    header.append('\t'.join(['MTD','peptide_search_engine_score', 'pyProphet q-value'])+'\n')
    header.append('\t'.join(['MTD','protein_quantification_unit', '[ , , Abundance, ]'])+'\n')
    header.append('\t'.join(['MTD','peptide_quantification_unit', '[ , , Abundance, ]'])+'\n')
    header.append('\t'.join(['MTD','variable_mod', 'library-based'])+'\n')
    header.append('\t'.join(['MTD','fixed_mod', 'library-based'])+'\n')

    for i in enumerate(col_idxs):
        file_org=df_I[df_I['Sample']==i[1]]['run'].tolist()[0]
        condition=df_I[df_I['Sample']==i[1]]['MSstats_Condition'].tolist()[0]
        header.append('\t'.join(['MTD','ms_run['+str(i[0]+1)+']-location',file_org])+'\n')
        header.append('\t'.join(['MTD','assay['+str(assayDict[condition])+']-ms_run_ref','ms_run['+str(i[0]+1)+']'])+'\n')
        header.append('\t'.join(['MTD','study_variable['+str(i[0]+1)+']-assay_ref','assay['+str(assayDict[condition])+']'])+'\n')
        header.append('\t'.join(['MTD','study_variable['+str(i[0]+1)+']-description',condition])+'\n')

    #output mzTab
    op=open(str(args.output),'w')
    for line in header:
        op.write(line)

    op.write('\n')
    for line in opr_prt:
        op.write(line)

    op.write('\n')
    for line in opr_pep:
        op.write(line)

    op.close() 

if __name__=="__main__":
    main()
