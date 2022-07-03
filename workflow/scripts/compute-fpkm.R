# Copyright © 2022 Nnlrl. All rights reserved.
# @Created by : Nnlrl
# @Created on : 2022-04-15
# @Title      : TODO
# @Objective  : TODO

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(GenomicFeatures)

txdb <- makeTxDbFromGFF(snakemake@input[["gtf"]], format="gtf")
ebg <- exonsBy(txdb, by="gene")

adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}

Counts <- read.table(file = snakemake@input[["counts"]],
                     header = TRUE, row.names = 1)

ebgList <- sum(width(reduce(ebg)))
genes <- intersect(rownames(Counts), names(ebgList))
Length <- as.vector(ebgList[genes])

TPM <- t(t(Counts / t(Length)) * 1e6 / colSums(Counts / t(Length)))
FPKM <- t(t(Counts / t(Length)) * 1e9 / colSums(Counts))

write.table(adjustdata(as.data.frame(TPM)), file=snakemake@output[["tpm"]],
            append=FALSE, quote=FALSE, sep="\t", row.names=F, col.names=T)

write.table(adjustdata(as.data.frame(FPKM)), file=snakemake@output[["fpkm"]],
            append=FALSE, quote=FALSE, sep="\t", row.names=F, col.names=T)
