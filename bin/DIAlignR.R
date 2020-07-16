#!/usr/bin/env Rscript

library(dplyr)
library(RSQLite)
library(DIAlignR)

alignTargetedRuns(dataPath='./', gapQuantile = 0.9, globalAlignmentFdr = 1, analyteFDR=1, unalignedFDR=1, alignedFDR=1)
