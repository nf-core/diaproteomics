#!/usr/bin/env Rscript

library(dplyr)
library(RSQLite)
library(DIAlignR)

args <- commandArgs()
p_globalAlignFDR <- args[6]
p_analyteFDR<-args[7]
p_unalignFDR<-args[8]
p_alignFDR<-args[9]
p_queryFDR<-args[10]

params<-paramsDIAlignR()
params$globalAlignmentFdr<-p_globalAlignFDR
params$analyteFDR<-p_analyteFDR
params$unalignFDR<-p_unalignFDR
params$alignFDR<-p_alignFDR
params$maxFdrQuery<-p_queryFDR
params$maxPeptideFdr<-p_queryFDR
params$XICfilter <- 'none'

# run dialignR with the specified FDR settings
alignTargetedRuns(dataPath='./', params=params)
