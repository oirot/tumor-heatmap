
genomes <- NULL


tumorHeatmap <- function(mpfFilePath, refGenome = 
                            BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, 
                          transcriptAnno = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
                          trDir = FALSE, 
                          type = "AlexV3.2",
                          numBases = 3,
                          enforceUniqueTrDir = TRUE,
                          signatures){
  #read the tumoral genomes
  genomes <- readGenomesFromMPF(mpfFilePath, numBases = numBases, 
                                type = "Alexandrov", trDir = trDir,
                                refGenome = refGenome, 
                                transcriptAnno = transcriptAnno, verbose = FALSE)
  
  if(is(signatures, "charachter")){
    #if the signatures are a string treat them as a path to the actual signatures
    signatures <- read.signature(signaturePath = signatures, signatureType = type)
  }
  #otherwise assume that the signature are in the object itself
    
  #TODO:remove
  assign("genomes", genomes, envir = .GlobalEnv)
  print("Genome loaded")
  exposure <- decomposeTumorGenomes(genomes = genomes, signatures = signatures)
  print("Exposure vecotrs estimated")
  exp <- do.call(cbind, exposure)
  heatmap(exp, Colv = NA, scale = "none")
  
  
}