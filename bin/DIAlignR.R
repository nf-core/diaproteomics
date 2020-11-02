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

# run dialignR with the specified FDR settings
alignTargetedRuns(dataPath='./', maxFdrQuery = p_queryFDR,  globalAlignmentFdr = p_globalAlignFDR, analyteFDR=p_analyteFDR, unalignedFDR=p_unalignFDR, alignedFDR=p_alignFDR)
