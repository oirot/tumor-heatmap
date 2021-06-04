
genomes <- NULL


heatmap.tumor <- function(mpfFilePath, refGenome = 
                            BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, 
                          transcriptAnno = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
                          trDir = FALSE, 
                          type = "AlexV3.2",
                          numBases = 3,
                          enforceUniqueTrDir = TRUE,
                          signatures){
  
  genomes <- readGenomesFromMPF(mpfFilePath, numBases = numBases, 
                                type = "Alexandrov", trDir = trDir,
                                refGenome = refGenome, 
                                transcriptAnno = transcriptAnno, verbose = FALSE)
  #TODO:remove
  assign("genomes", genomes, envir = .GlobalEnv)
  print("Genome loaded")
  signatures <- read.signature(signaturePath = signatures, signatureType = type)
  print("Signature loaded")
  exposure <- decomposeTumorGenomes(genomes = genomes, signatures = signatures)
  print("Exposure vecotrs estimated")
  exp <- do.call(cbind, exposure)
  #exp <- apply(exp, MARGIN = 1, function(sig) {
  #  sig - mean(sig) * sd(sig)
  #})
  heatmap(exp, Colv = NA, scale = "none")
  
  
}