#!/usr/bin/env python
from __future__ import print_function
import sys
import scipy
from numpy import *
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import glob
import re
import networkx as nx
from itertools import combinations, product
from scipy.interpolate import interp1d
import argparse
import logging


"""
merge_and_align_libraries_from_easypqp.py:
This script takes multiple spectral libraries as input and can concacenate them into a single merged library.
If specified a linear RT alignment is carried out between the libraries, in order to bring them in the same RT space.
"""

__author__      = "Leon Bichmann"


#logging setup
console = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console.setFormatter(formatter)
LOG = logging.getLogger("Align RT of spectral libraries")
LOG.addHandler(console)
LOG.setLevel(logging.INFO)

# Linear retention time alignment between two libraries (reference and other)
def align_libs(reference, other, rsq_threshold):

    # read first library, groupby modified sequence and charge and store RTs
    df = reference
    df_rt = df.groupby(['ModifiedPeptideSequence', 'PrecursorCharge'])['NormalizedRetentionTime'].apply(
        median).reset_index()
    df_rt['Key'] = df_rt['ModifiedPeptideSequence'] + df_rt['PrecursorCharge'].apply(str)
    df_rt = df_rt[['Key', 'NormalizedRetentionTime']]

    # read second library, groupby modified sequence and charge and store RTs
    df_II = other
    df_rt_II = df_II.groupby(['ModifiedPeptideSequence', 'PrecursorCharge'])['NormalizedRetentionTime'].apply(
        median).reset_index()
    df_rt_II['Key'] = df_rt_II['ModifiedPeptideSequence'] + df_rt_II['PrecursorCharge'].apply(str)
    df_rt_II = df_rt_II[['Key', 'NormalizedRetentionTime']]

    # calibrate rt alignment
    df_rt_merged = pd.merge(df_rt, df_rt_II, on='Key')
    slope, intercept, r_value, p_value, std_err = stats.linregress(df_rt_merged['NormalizedRetentionTime_x'], df_rt_merged['NormalizedRetentionTime_y'])
    (a, b) = (slope, intercept)
    rsq = r_value ** 2

    # raise an error if the alignment is below a specified rsq threshold
    if rsq < rsq_threshold:
        raise Exception("Error: R-squared " + str(rsq) + " is below the threshold of " + str(rsq_threshold) + ".")

    # apply rt transformation
    for i, row in df_II.iterrows():
        rt = row["NormalizedRetentionTime"]
        new_rt = scipy.polyval([a, b], rt)
        df_II.loc[i, "NormalizedRetentionTime"] = new_rt

    return df_II


# Compute a minimum spanning tree (MST) across multiple libraries based on their shared peptide overlap
def compute_MST(libs,min_overlap):
    file_combs = combinations(libs, 2)
    G=nx.Graph()

    #compute peptide shares across all libraries
    for file_comb in file_combs:

        lib_I=pd.read_csv(file_comb[0], sep='\t') #.split('/')[-1].split('_Class')[0][3:]
        lib_II=pd.read_csv(file_comb[1], sep='\t') #.split('/')[-1].split('_Class')[0][3:]

        overlap = list(set(lib_I['ModifiedPeptideSequence'].values.tolist()) & set(lib_II['ModifiedPeptideSequence'].values.tolist()))

        # store them as connected graph, if the share is greater than a specified threshold
        if file_comb[0] not in G.nodes():
            G.add_node(file_comb[0])
        if file_comb[1] not in G.nodes():
            G.add_node(file_comb[1])

        if len(overlap)>=min_overlap:
            G.add_edges_from([(file_comb[0],file_comb[1],{'weight':float(10000)/len(overlap)})]) # the graph edge lengths are antiproportional to the peptide share

        # generate minimum spanning tree
        T = nx.minimum_spanning_tree(G)

    return T

# Carry out a pairwise RT alignment and concatenate multiple libraries along the edges of the computed minimum spanning tree of peptide share.
def combine_libs_by_edges_of_MST(T, rsq_threshold):

    #define the center of the MST as source_file
    try:
        source_file = nx.center(T)[0]

    except:
        print('It was not possible to align all libraries into the same RT space, since they don\'t share enough peptides between all runs!')
        print('All libraries must pairwise connect into a complete graph by shared peptides.')
        print('There might be one outlier sample that has no overlap with the others or the samples are too distant.')
        print('You can try lowering the parameter --overlap_for_merging, but be aware that this is the least number of peptides used for linear RT alignment!')
        sys.exit()


    #collect shortest paths from center to all nodes in MST
    short_paths=[]
    for node_combi in sorted(product([source_file], T.nodes())):

        if node_combi[0] != node_combi[1]:
            short_paths.append(nx.shortest_path(T, node_combi[0], node_combi[1]))

    short_paths.sort(key=len)

    #align all libraries to source file along shortest paths of MST
    #start with all paths of length 2
    print('>>> align libraries')
    outfiles = {}
    for path in [p for p in short_paths if len(p) == 2]:
        print(path)
        reference = pd.read_csv(path[0], sep='\t')
        other = pd.read_csv(path[1], sep='\t')

        cols=[c for c in reference.columns if not 'NormalizedRetentionTime' in c and not 'TransitionId' in c]
        aligned = align_libs(reference, other, rsq_threshold)
        concat = pd.concat([reference,aligned])
        regrouped = concat.iloc[concat[cols].drop_duplicates().index]

        outfile = './' + path[0].split('/')[-1].split('.tsv')[0] + '_' + path[1].split('/')[-1]
        regrouped.to_csv(outfile, index=False, sep='\t')
        outfiles[''.join(path)] = outfile

    #continue with all paths of length longer two
    for path in [p for p in short_paths if len(p) > 2]:
        print(path)
        reference = pd.read_csv(outfiles[''.join(path[:len(path) - 1])], sep='\t')
        other = pd.read_csv(path[-1], sep='\t')

        aligned = align_libs(reference, other, rsq_threshold)
        concat = pd.concat([reference,aligned])
        regrouped = concat.iloc[concat[cols].drop_duplicates().index]

        outfile = './' + outfiles[''.join(path[:len(path) - 1])].split('/')[-1].split('.tsv')[0] + '_' + path[-1].split('/')[-1]
        regrouped.to_csv(outfile, index=False, sep='\t')
        outfiles[''.join(path)] = outfile

    #group all aligned merged libs and drop duplicates
    outfile_list=[]
    for o in outfiles.values():
        df_out=pd.read_csv(o, sep='\t')
        outfile_list.append(df_out)

    concat = pd.concat(outfile_list)
    concat = concat.drop('TransitionId',axis='columns')
    regrouped = concat.iloc[concat[cols].drop_duplicates().index]
    combined_lib = regrouped.drop_duplicates()

    #assign new transition ids
    combined_lib['TransitionId']=range(0,combined_lib.shape[0])

    return combined_lib


# concatenate libraries without carriying out an RT alignment between them (already aligned or no overlap)
def concatenate_without_alignment(files):
    reference = pd.read_csv(files[0], sep='\t')

    for other in files[1:]:
       other_df = pd.read_csv(other, sep='\t')

       cols=[c for c in reference.columns if not 'NormalizedRetentionTime' in c and not 'TransitionId' in c]
       concat = pd.concat([reference,other_df])
       regrouped = concat.iloc[concat[cols].drop_duplicates().index]

       reference = regrouped

    combined_lib = regrouped.drop_duplicates()
    combined_lib['TransitionId']=range(0,combined_lib.shape[0])

    outfile = './' + files[0].split('/')[-1].split('.tsv')[0] + '_concatenated_lib.tsv'
    combined_lib.to_csv(outfile, index=False, sep='\t')

    return combined_lib


def main():
    model = argparse.ArgumentParser(description='Postprocess Neoepitopes predicted by MHCNuggets')

    model.add_argument(
        '-i', '--input_libraries',
        type=str,
        nargs='*',
        help='library file serving as reference'
    )

    model.add_argument(
        '-f', '--filter',
        type=str,
        nargs='*',
        help='peptide list to filter library with'
    )

    model.add_argument(
        '-u', '--min_overlap',
        type=int,
        help='min number of peptides having to overlap between libraries'
    )

    model.add_argument(
        '-a', '--align',
        type=str,
        help='Whether an alignment should be carried out or just the merging of multiple libraries'
    )

    model.add_argument(
        '-t', '--rsq_threshold',
        type=float,
        help='rsq threshold for alignment'
    )

    model.add_argument(
        '-o', '--output',
        type=str,
        help='output aligned library file'
    )

    args = model.parse_args()

    libs=args.input_libraries
    rsq_threshold=args.rsq_threshold
    min_overlap=args.min_overlap
    align=args.align

    if align=='true':
       if len(libs)>1:
          MST=compute_MST(libs, min_overlap)
          nx.draw(MST, with_labels=True, node_size=1, font_size=5, alpha=0.4)
          plt.savefig(args.output.replace('.tsv','.png'),dpi=(1000))

          combined_lib=combine_libs_by_edges_of_MST(MST, rsq_threshold)

       else:
          combined_lib=pd.read_csv(libs[0], sep='\t')

    else:
       combined_lib=concatenate_without_alignment(libs)

    #output transformed dataframe II
    if args.filter:
        combined_lib=combined_lib[combined_lib['PeptideSequence'].isin(args.filter)]

    combined_lib.to_csv(args.output, index=False, sep='\t')


if __name__=="__main__":
   main()

