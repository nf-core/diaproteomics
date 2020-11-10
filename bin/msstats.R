#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

usage <- "Rscript msstats.R input.csv [list of contrasts or 'pairwise'] [default control condition or ''] [output prefix]"

df_list = list.files(path = "./", pattern = "*.csv")
print(length(df_list))
if (length(df_list)<1) {
  print(usage)
  stop("At least one arguments must be supplied (input csv).n", call.=FALSE)
}
if (length(df_list)<=2) {
  # contrasts
  args[3] = "pairwise"
  args[4] = ""
}
if (length(df_list)<=3) {
  # default control condition
  args[4] = ""
}
if (length(df_list)<=4) {
  # default output prefix
  args[5] = "msstats"
}

df <- read.csv(paste0("./",df_list[1]), sep=',')
print(df_list[1])
if (length(df_list)>1){
  for (df_name in df_list[2:2]) {
    print(df_name)
    df_add <- read.csv(paste0("./",df_name), sep=',')
    df <- rbind(df, df_add)
  }
}

csv_input <- df
contrast_str <- args[3]
print(contrast_str)
control_str <- args[4]
print(control_str)
out_prefix <- args[5]
folder <- dirname(args[1])
filename <- basename(args[1])
mzTab_output <- paste0(folder,'/',out_prefix,filename)

# load the MSstats library 
require(MSstats)
require(dplyr)
require(tidyr)

# read dataframe into MSstats
data <- df
print('load data')
quant <- OpenMStoMSstatsFormat(data,
                               removeProtein_with1Feature = FALSE)

# process data
print('process data')
processed.quant <- dataProcess(quant, censoredInt = 'NA', featureSubset='highQuality')

lvls <- levels(as.factor(data$Condition))
if (length(lvls) == 1)
{
  print("Only one condition found. No contrasts to be tested. If this is not the case, please check your experimental design.")
} else {
  print(length(lvls))
  if (contrast_str == "pairwise")
  {
    if (control_str == "")
    {
      print(lvls)
      l <- length(lvls)
      contrast_mat <- matrix(nrow = l * (l-1) / 2, ncol = l)
      rownames(contrast_mat) <- rep(NA, l * (l-1) / 2)
      colnames(contrast_mat) <- lvls
      c <- 1
      for (i in 1:(l-1))
      {
        for (j in (i+1):l)
        {
          comparison <- rep(0,l)
          comparison[i] <- -1
          comparison[j] <- 1
          contrast_mat[c,] <- comparison
          rownames(contrast_mat)[c] <- paste0(lvls[i],"-",lvls[j])
          c <- c+1
        }
      }
    } else {
      control <- which(as.character(lvls) == control_str)
      if (length(control) == 0)
      {
        stop("Control condition not part of found levels.n", call.=FALSE)
      }
      
      l <- length(lvls)
      contrast_mat <- matrix(nrow = l-1, ncol = l)
      rownames(contrast_mat) <- rep(NA, l-1)
      colnames(contrast_mat) <- lvls
      c <- 1
      for (j in setdiff(1:l,control))
      {
        comparison <- rep(0,l)
        comparison[i] <- -1
        comparison[j] <- 1
        contrast_mat[c,] <- comparison
        rownames(contrast_mat)[c] <- paste0(lvls[i],"-",lvls[j])
        c <- c+1
      }
    }
  } else {
    print("Specific contrasts not supported yet.")
    exit(1)
  }
  
  print ("Contrasts to be tested:")
  print (contrast_mat)
  #TODO allow for user specified contrasts
  test.MSstats <- groupComparison(contrast.matrix=contrast_mat, data=processed.quant)
  
  #TODO allow manual input (e.g. proteins of interest)
  write.csv(test.MSstats$ComparisonResult, paste0(filename, "_msstats_results.csv"))
  
  groupComparisonPlots(data=test.MSstats$ComparisonResult, type="ComparisonPlot",
                       width=12, height=12,dot.size = 2)
  
  test.MSstats$Volcano = test.MSstats$ComparisonResult[!is.na(test.MSstats$ComparisonResult$pvalue),]
  groupComparisonPlots(data=test.MSstats$Volcano, type="VolcanoPlot",
                       width=12, height=12,dot.size = 2,ProteinName=FALSE)

  # Otherwise it fails since the behaviour is undefined
  if (nrow(contrast_mat) > 1)
  {
    groupComparisonPlots(data=test.MSstats$ComparisonResult, type="Heatmap",
                         width=12, height=12,dot.size = 2)
  }
  
  #for (comp in rownames(contrast_mat))
  #{
  #  groupComparisonPlots(data=test.MSstats$ComparisonResult, type="ComparisonPlot",
  #                       width=12, height=12,dot.size = 2, sig=1),
  #                       which.Comparison = comp,
  #                       address=F)
  #  # try to plot all comparisons
  #}

  
  # annotate how often the protein was quantified in each condition (NA values introduced by merge of completely missing are set to 1.0)
  ############ also calculate missingness on condition level

  # input: ProcessedData matrix of MSstats
  # output: 
  #   calculate fraction of na in condition (per protein)
  # Groups:   PROTEIN [762]
  #   PROTEIN                 `1`   `2`
  #   <fct>                 <dbl> <dbl>
  # 1 sp|A1ANS1|HTPG_PELPD   0    0.5  
  # 2 sp|A2I7N3|SPA37_BOVIN  0    0.5  
  # 3 sp|A2VDF0|FUCM_HUMAN   0    0.5  
  # 4 sp|A6ND91|ASPD_HUMAN   0.5  0.5  
  # 5 sp|A7E3W2|LG3BP_BOVIN  0.5  0.5  
  # 6 sp|B8FGT4|ATPB_DESAA   0    0.5

  getMissingInCondition <- function(processedData)
  {
    p <- processedData

    print('start missingness')

    # count number of samples per condition
    n_samples = p %>% group_by(GROUP) %>% summarize(n_samples = length(unique((as.numeric(SUBJECT))))) 
    print('summary done')

    p <- p %>% 
     filter(!is.na(INTENSITY)) %>% # remove rows with INTENSITY=NA
     select(PROTEIN, GROUP, SUBJECT) %>%
     distinct() %>% 
     group_by(PROTEIN, GROUP) %>% 
     summarize(non_na = n())  # count non-NA values for this protein and condition
    print('prepare for join')

    p <- left_join(p, n_samples) %>% 
         mutate(missingInCondition = 1 - non_na/n_samples) # calculate fraction of missing values in condition
    print('join done')

    # create one column for every condition containing the missingness
    p <- spread(data = p[,c("PROTEIN", "GROUP", "missingInCondition")], key = GROUP, value = missingInCondition)
    return(p)
  }

  mic <- getMissingInCondition(processed.quant$ProcessedData)

  print('missing done')

  test.MSstats$ComparisonResult <- merge(x=test.MSstats$ComparisonResult, y=mic, by.x="Protein", by.y="PROTEIN")

  print('comparison done')

  commoncols <- intersect(colnames(mic), colnames(test.MSstats$ComparisonResult))

  test.MSstats$ComparisonResult[, commoncols]<-test.MSstats$ComparisonResult %>% select(all_of(commoncols)) %>% mutate_all(list(replace = function(x){replace(x, is.na(x), 1)})) 

  print('aggregate results done')

  #write comparison to CSV (one CSV per contrast)              
  writeComparisonToCSV <- function(DF) 
  {
    write.table(DF, file=paste0(filename,"_comparison_",unique(DF$Label),".csv"), quote=FALSE, sep='\t', row.names = FALSE)
    return(DF)
  }

  print('write done')
  test.MSstats$ComparisonResult %>% group_by(Label) %>% do(writeComparisonToCSV(as.data.frame(.)))
}
