#!/usr/bin/env Rscript

library(ggplot2)
library(gplots)
library(reshape2)

#parse commands
args <- commandArgs()
Sample <- args[6]
#Condition<- args[2]




# plot library RT deviation
df_list = list.files(path = "./", pattern = "*.mzML.tsv")
df <- read.csv(paste0("./",df_list[1]), sep='\t')

for (df_name in df_list[2:2])
  print(df_name)
  df_add <- read.csv(paste0("./",df_name), sep='\t')
  df <- rbind(df, df_add)
  

df <- df[which(df$decoy==0),]

pdf(paste0(Sample, '_rt_deviation.pdf'))
ggplot(df, aes(RT, delta_rt)) + geom_bin2d(bins = 500) + scale_fill_viridis_c() + ylim( -500, 500) + theme_bw()
dev.off()

# plot peptide and protein counts
data <- data.frame(
  Evidence=c('Peptides','Proteins'),
  Count=c(length(unique(df$Sequence)),length(unique(df$ProteinName)))
)

pdf(paste0(Sample,'_protein_peptide_stats.pdf'))
ggplot(data, aes(x="Evidence", y=Count, fill=Evidence)) + geom_bar(position="dodge", stat="identity", width=1) + theme_bw()
dev.off()


# plot peptide and protein counts
data <- data.frame(
  Charge=names(table(df$Charge)),
  Value=as.numeric(table(df$Charge))
)

pdf(paste0(Sample, '_charge_stats.pdf'))
ggplot(data, aes(x="", y=Value, fill=Charge)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + theme_bw()
dev.off()


# plot peptide quantities across runs
df_list = list.files(path = "./", pattern = "*quantities.csv")

for (df_name in df_list)
  df <- read.csv(paste0("./",df_name), sep=',')

  df$intensity <- log10(df$intensity)

  pivot <- as.data.frame(dcast(df, sequence ~ run, value.var = "intensity", fun.aggregate = sum), index='sequence', stringsAsFactors = FALSE)
  rownames(pivot) <- pivot$sequence
  pivot$sequence <- NULL
  pivot[mapply(is.infinite, pivot)] <- NA
  pivot[is.na(pivot)] <- 0

  pdf(gsub(df_name,pattern = '.csv', replacement = '.pdf'))
  heatmap.2(data.matrix(pivot), xlab = "MS runs", ylab = "Peptides", labRow = FALSE, cexCol=0.15, srtCol=45, trace="none", key.title = 'Intensity range', col = colorRampPalette(c("darkblue","white","darkred"))(100), margins=c(5,8))
  dev.off()
