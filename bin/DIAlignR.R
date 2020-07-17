#!/usr/bin/env Rscript

library(dplyr)
library(RSQLite)
library(DIAlignR)

args <- commandArgs()
p_globalAlignFDR <- args[6]
p_analyteFDR<-args[7]
p_unalignFDR<-args[8]
p_alignFDR<-args[9]

alignTargetedRuns(dataPath='./', globalAlignmentFdr = p_globalAlignFDR, analyteFDR=p_analyteFDR, unalignedFDR=p_unalignFDR, alignedFDR=p_alignFDR)
