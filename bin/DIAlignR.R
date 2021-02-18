#!/usr/bin/env Rscript

#load libraries
library(tools)
library(dplyr)
library(RSQLite)
library(DIAlignR)
library(BiocParallel)

#read out params from commandline execution
args <- commandArgs()
p_globalAlignFDR <- args[6]
p_analyteFDR<-args[7]
p_unalignFDR<-args[8]
p_alignFDR<-args[9]
p_queryFDR<-args[10]
p_level<-args[11]
xic_filter<-args[12]
parallel<-args[13]
workers<-args[14]

#set params for DIAlignR
params<-paramsDIAlignR()
params$globalAlignmentFdr<-p_globalAlignFDR
params$analyteFDR<-p_analyteFDR
params$unalignFDR<-p_unalignFDR
params$alignFDR<-p_alignFDR
params$maxFdrQuery<-p_queryFDR
params$maxPeptideFdr<-p_queryFDR
p_level<-toTitleCase(p_level)
params$level<-p_level
params$XICfilter<-xic_filter

if (parallel=='parallel'){
#configure mutlicore execution
BiocParallel::register(BiocParallel::MulticoreParam(workers = workers))

# run dialignR with the specified FDR settings with paralellization
alignTargetedRuns(dataPath='./', params=params, applyFun = BiocParallel::bplapply, oswMerged = TRUE)

} else {

# run dialignR with the specified FDR settings without parallelization
alignTargetedRuns(dataPath='./', params=params, oswMerged = TRUE)
}
