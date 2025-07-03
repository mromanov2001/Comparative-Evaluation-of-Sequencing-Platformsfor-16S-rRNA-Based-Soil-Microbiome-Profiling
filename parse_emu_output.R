args <- commandArgs(TRUE)

# load args
work_dir <- args[[1]]
metadata <- read.csv(args[[2]], sep = ",", header = T, colClasses = "character")
emu_table <- c(args[[3]])
taxonomy_table <- c(args[[4]])
emu_report_table <- c(args[[5]])

# set work directory
setwd(work_dir)

# import libraries
library(dplyr)
library(reshape2)
library(tidyr)

# make taxonomy and abundance tables
df.all <- NULL
df.tax.all <- NULL

for (f in metadata$barcode){
    df <- read.csv(paste0("emu/", f, "/", f, "_rel-abundance.tsv"), sep = "\t")

    # make taxonomy table
    df.tax <- df[c(10:3)]
    df.tax.all <- rbind(df.tax, df.tax.all)

    df <- df[c(3,14)]
    df$estimated.counts <- round(df$estimated.counts)
    colnames(df)[2] <- "abundance"
    df <- df[-nrow(df),]
    df$barcode <- sub("_rel-abundance.tsv", "", basename(f))

    # # make abundance table
    df.all <- rbind(df.all, df)
}

# reformat tables
## relab table
df.all <- spread(df.all, species, abundance, fill = 0)
rownames(df.all) <- df.all$barcode
df.all <- df.all[-1]
df.all <- df.all[order(colSums(df.all), decreasing = T)]
colnames(df.all) <- gsub(" ", "_", colnames(df.all))
colnames(df.all) <- gsub("\\[|\\]", "", colnames(df.all))
colnames(df.all) <- gsub("-", "_", colnames(df.all))

## tax table
df.tax.all <- unique(df.tax.all)
df.tax.all$species <- gsub(" ", "_", df.tax.all$species)
df.tax.all$species <- gsub("\\[|\\]", "", df.tax.all$species)
df.tax.all$species <- gsub("-", "_", df.tax.all$species)
df.tax.all <- df.tax.all[df.tax.all$species %in% colnames(df.all),]
df.tax.all <- df.tax.all[-2]
df.tax.all <- df.tax.all[c(7:1)]

# make emu report table
emu_report <- as.data.frame(t(df.all))
emu_report <- merge(df.tax.all, cbind(sampleid = row.names(emu_report), emu_report), by = 1)
emu_report <- emu_report[c(7:1,8:ncol(emu_report))]
emu_report <- rbind(c(rep(NA, 7), metadata$sampleid), emu_report)

# write tables
write.table(df.all, emu_table, sep = "\t", quote = F)
write.table(df.tax.all, taxonomy_table, sep = "\t", quote = F, row.names = F)
write.table(emu_report, emu_report_table, sep = "\t", quote = F, row.names = F)