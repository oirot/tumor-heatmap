#setting defaults for reference genome and annotations
readTumoralGenomes <- function(mpfFilePath, referenceGenome = 
                                 BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, 
                               annotations = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene){
  genomes <- readGenomesFromMPF(mpfFilePath, 
                     numBases = 3, type = "Alexandrov", trDir = FALSE, 
                     refGenome = referenceGenome, 
                     transcriptAnno = annotations, 
                     verbose = FALSE)
}