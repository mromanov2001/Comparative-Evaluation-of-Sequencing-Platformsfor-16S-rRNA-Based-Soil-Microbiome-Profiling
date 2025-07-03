args <- commandArgs(TRUE)

# load parameters
work_dir <- args[[1]]
metadata <- read.csv(args[[2]], sep = ",", header = T, colClasses = "character")
statistic <- c(args[[3]])

# set work directory
setwd(work_dir)

# import libraries
library(dplyr)
library(reshape2)
library(tidyr)

raw <- list.files('stat/raw/')
filtered <- list.files('stat/raw__trimmed__filtered/')

df.raw <- NULL
for (i in raw){
    df <- read.csv(paste0('stat/raw/', i), sep = "\t")[c(1,4,5),]
    df$barcode <- sub(".tsv", "", i)
    df.raw <- rbind(spread(df, Metrics, dataset)[c(1,3,2,4)], df.raw)
}

colnames(df.raw)[-1] <- paste0("raw_", colnames(df.raw)[-1])

df.filtered <- NULL
for (i in filtered){
    df <- read.csv(paste0('stat/raw__trimmed__filtered/', i), sep = "\t")[c(1,4,5),]
    df$barcode <- sub(".tsv", "", i)
    df.filtered <- rbind(spread(df, Metrics, dataset)[c(1,3,2,4)], df.filtered)
}

colnames(df.filtered)[-1] <- paste0("filtered_", colnames(df.filtered)[-1])

df.stat <- merge(df.raw, df.filtered, by = 1)
df.stat$perc_filtered <- 100*(1-(as.numeric(df.stat$filtered_number_of_reads)/as.numeric(df.stat$raw_number_of_reads)))
df.stat <- merge(metadata, df.stat, by = 1)

write.table(df.stat, statistic, sep = "\t", quote = F, row.names = F)