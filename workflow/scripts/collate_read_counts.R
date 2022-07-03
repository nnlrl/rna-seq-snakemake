# Collate read counts into one big matrix

args <- commandArgs(trailingOnly = TRUE)

input_dir <- args[1]
out_file <- args[2]

count_files <- dir(input_dir, pattern = ".read_counts.csv$", full.names = TRUE)

# get read counts for each sample into a list
counts <- lapply(count_files, function(f) {
  x <- data.table::fread(f)[,c(1,7)]
  names(x)[2] <- sub(".*/(.*?)_Aligned\\.sortedByCoord\\.out\\.bam", "\\1", names(x)[2])
  x
})

# merge list of data frames
counts_all <- as.data.frame(Reduce(function(dtf1, dtf2)
  merge(dtf1, dtf2, by = "Geneid", all.x = TRUE),
       counts))
rownames(counts_all) <- counts_all$Geneid
counts_all$Geneid <- NULL

# save results to out file
write.table(counts_all, out_file, quote = FALSE,
            sep = '\t')
