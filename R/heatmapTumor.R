
tumorHeatmap <- function(mpfFilePath, 
                         signatures,
                         refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, 
                         transcriptAnno = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
                         trDir = FALSE, 
                         type = signatureTypes$alexandrov32,
                         numBases = 3,
                         enforceUniqueTrDir = TRUE,
                         verbose = FALSE){
  
  if(is(signatures, "character")){
    #if the signatures are a string treat them as a path to the actual signatures
    signatures <- readSignatures(signaturePath = signatures, signatureType = type)
  }
  #otherwise assume that the signature are in the object itself
  if(!isSignatureSet(signatures)){
    stop("Invalid signatures")
  }
  
  
  if(verbose) print("Loading genomes... ")
  #read the tumoral genomes
  genomes <- readGenomesFromMPF(mpfFilePath, numBases = numBases, 
                                type = getGenomeType(type), trDir = trDir,
                                refGenome = refGenome, 
                                transcriptAnno = transcriptAnno, verbose = FALSE)
  if(verbose) print("DONE")
  
  if(verbose) print("Estimating exposure vecotrs...")
  exposure <- decomposeTumorGenomes(genomes = genomes, signatures = signatures)
  if(verbose) print("DONE")
  
  exp <- do.call(cbind, exposure)
  heatmap(t(exp), Colv = NA, scale = "row", revC = TRUE, 
          main = "Heatmap of the exposure vectors", 
          xlab = "Signature", ylab = "Sample")
}