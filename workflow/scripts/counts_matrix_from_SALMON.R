# R script takes outputs from salmon and uses tximport package to create a matrix
# that is genes x samples

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
salmon_output_folder <- args[1] #folder where salmon files are written
counts_dir <- args[2] #where to save to collated count files
colDataFile <- args[3]

writeCounts <- function(colDataFile, salmon_output_folder, counts_dir, type) {
  colData <- read.table(colDataFile)
  if(type == 'transcripts') {
    files <- file.path(salmon_output_folder, rownames(colData), "quant.sf")
  } else if ( type == 'genes') {
    files <- file.path(salmon_output_folder, rownames(colData), "quant.genes.sf")
  }
  names(files) <- rownames(colData)
  txi <- tximport::tximport(files, type = "salmon", txOut = TRUE)

  dds <- DESeq2::DESeqDataSetFromTximport(txi, colData, ~group)

  #save raw counts
  write.table(x = DESeq2::counts(dds),
              file = file.path(counts_dir, paste0("counts_from_SALMON.", type,".tsv")),
              quote = F, sep = '\t')
}

writeTPMcounts <- function(colDataFile, salmon_output_folder, counts_dir, type) {
  colData <- read.table(colDataFile)
  if(type == 'transcripts') {
    files <- file.path(salmon_output_folder, rownames(colData), "quant.sf")
  } else if ( type == 'genes') {
    files <- file.path(salmon_output_folder, rownames(colData), "quant.genes.sf")
  }
  names(files) <- rownames(colData)

  tpmTables <- lapply(names(files), function(n) {
    f <- files[[n]]
    df <- subset(fread(f), select = c('Name', 'TPM'))
    colnames(df) <- c('Name', n)
    return(df)
  })
  names(tpmTables) <- names(files)

  #tpmMatrix = data.frame(Reduce(function(...) merge(..., all = TRUE, by = c('Name')), tpmTables))
  tpmMatrix = Reduce(function(...) merge(..., all = TRUE, by = c('Name')), tpmTables)

  rownames(tpmMatrix) <- tpmMatrix$Name
  tpmMatrix <- tpmMatrix[-1]
  write.table(x = tpmMatrix,
              file = file.path(counts_dir, paste0("TPM_counts_from_SALMON.", type,".tsv")),
              quote = F, sep = '\t')
}

raw_counts_dir <- file.path(counts_dir, "raw_counts", "salmon")
writeCounts(colDataFile, salmon_output_folder, raw_counts_dir, type = 'transcripts')
writeCounts(colDataFile, salmon_output_folder, raw_counts_dir, type = 'genes')

normalized_counts_dir <- file.path(counts_dir, "normalized", "salmon")
writeTPMcounts(colDataFile, salmon_output_folder, normalized_counts_dir, type = 'transcripts')
writeTPMcounts(colDataFile, salmon_output_folder, normalized_counts_dir, type = 'genes')
