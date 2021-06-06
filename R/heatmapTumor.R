#' Computes the exposure vectors from an MPF file and plots 
#' an heatmap clustering the different tumor genomes by their 
#' similarities in the calculated exposure different tumarl genomes by their 
#'
#' @param mpfFilePath (Mandatory) The path of the MPF file from which to load the tumor
#' genomes
#'  
#' @param signatures (Mandatory) Either a path to the signature file or a vector of 
#' path to signatures file (in the case of Shiraishi signatures). Also, it can be
#' dircetly an object containing the signatures. The object can be A list of 
#' vectors, data frames or matrices. Vectors are used for Alexandrov signatures, 
#' data frames or matrices for Shiraishi signatures.
#'
#' @param signaturesTypes (Optional) A string representing the signature type
#' a list of signature types is available as the object `signatureTypes` 
#' @seealso [signatureTypes]. Use alexandrov2 for Alexandrov V2 signatures, 
#' alandrov32 for Alexandrov V3.2 ones and shiraishi for Shiraishi signatures.
#'
#' @param refGenome (Mandatory) Total number of bases (mutated base and flanking
#' bases) to be used for sequence patterns. Must be odd. Default: 5
#'
#' @param trDir (Mandatory) Specifies whether the transcription direction is 
#' taken into account in the signature model. If so, only mutations within
#' genomic regions with a defined transcription direction can be 
#' considered. Default: TRUE
#' 
#' @param enforceUniqueTrDir (Optional) Used only if trDir is TRUE. If
#' enforceUniqueTrDir is TRUE (default), then mutations which map to a region
#' with multiple overlapping genes with opposing transcription directions will
#' be excluded from the analysis. If FALSE, the transcript direction 
#' encountered first in the transcript database (see transcriptAnno) is
#' assigned to the mutation. 
#' 
#' @param refGenome refGenome	(Mandatory) The reference genome (BSgenome)
#' needed to extract sequence patterns. Default: BSgenome object for hg19.
#' 
#' @param transcriptAnno (Optional) Transcript annotation (TxDb object) used to
#' determine the transcription direction. This is required only if trDir is 
#' TRUE. Default: TxDb object for hg19.
#' 
#' @param verbose (Optional) Print information about reading and processing the
#' mutation data. Default: TRUE
#' 
#' @details
#' Many of the parameters of this function are also described in the function
#' `readGenomesFromMPF{decompTumor2Sig}` @seealso [readGenomesFromMPF()] 
#'  @seealso [decomposeTumorGenomes()] as that function and package are at the
#' base of this function and many parameteres are passed directly to those
#' functions. A more detailed explanation of the parameters and the function can 
#' be found there.
#' 
#' @return The exposure vector is returned and a heatmap of their similarity is
#' plotted
#' 
#' @importFrom decompTumor2Sig isSignatureSet
#' @importFrom decompTumor2Sig readGenomesFromMPF
#' @importFrom decompTumor2Sig decomposeTumorGenomes
#' @export


tumorHeatmap <- function(mpfFilePath, 
                         signatures,
                         signaturesType = signatureTypes$alexandrov32,
                         numBases = 5,
                         trDir = TRUE, 
                         enforceUniqueTrDir = TRUE,
                         refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, 
                         transcriptAnno = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
                         verbose = FALSE){
  
  if(is(signatures, "character")){
    #if the signatures are a string treat them as a path to the actual signatures
    signatures <- readSignatures(signaturePath = signatures, type = signaturesType)
  }
  #otherwise assume that the signature are in the object itself
  if(!isSignatureSet(signatures)){
    stop("Invalid signatures")
  }
  
  
  if(verbose) print("Loading genomes... ")
  #read the tumor genomes
  genomes <- decompTumor2Sig::readGenomesFromMPF(mpfFilePath, numBases = numBases, 
                                type = getGenomeType(signaturesType), trDir = trDir,
                                refGenome = refGenome, 
                                transcriptAnno = transcriptAnno,
                                verbose = verbose)
  if(verbose) print("DONE")
  if(!sameSignatureFormat(signatures, genomes)){
    stop(paste0("The parameter used to load the genomes need to be the same",
    " as the ones used to generate the signatures. Please verify that this is ",
    "the case"))
  }
  if(verbose) print("Estimating exposure vecotrs...")
  exposure <- decompTumor2Sig::decomposeTumorGenomes(genomes = genomes, signatures = signatures)
  if(verbose) print("DONE")
  
  exposures <- t(do.call(cbind, exposure))
  heatmap(exposures, Colv = NA, scale = "row", revC = TRUE, 
          main = "Heatmap of the exposure vectors", 
          xlab = "Signature", ylab = "Sample")
  
  return(exposures)
}