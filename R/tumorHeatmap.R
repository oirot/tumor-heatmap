##read the signature from a file or if its a proper signature obj do nothing
.readSignaturesFromFileOrObject <-
  function(signatures, signaturesType) {
    if (is(signatures, "character")) {
      ##if the signatures are a string treat them as a path to
      ##the actual signatures
      signatures <- .readSignatures(signaturePath=signatures,
                                   type=signaturesType)
    }
    ##otherwise assume that the signature are in the object itself
    if (!isSignatureSet(signatures)) {
      stop("Invalid signatures")
    }
    
    return(signatures)
  }


##read the genome from the MPF file
.readGenome <- function(mpfFilePath,
                       numBases,
                       signaturesType,
                       trDir,
                       refGenome,
                       transcriptAnno,
                       verbose) {
  if (verbose)
    message("Loading genomes... ")
  ##read the tumor genomes
  genomes <- readGenomesFromMPF(
    mpfFilePath,
    numBases=numBases,
    type=.getGenomeType(signaturesType),
    trDir=trDir,
    refGenome=refGenome,
    transcriptAnno=transcriptAnno,
    verbose=verbose
  )
  if (verbose)
    message("DONE")
  
  return(genomes)
}


##compute the exposure vectors
.computeExposures <- function(genomes, signatures, verbose) {
  if (verbose)
    message("Estimating exposure vecotrs...")
  exposure <-
    decomposeTumorGenomes(genomes=genomes, signatures=signatures)
  exposure <- t(do.call(cbind, exposure))
  if (verbose)
    message("DONE")
  
  return(exposure)
}


.plotHeatmapAndDendogram <- function(exposures) {
  heatmap(
    exposures,
    Colv=NA,
    scale="row",
    revC=TRUE,
    main="Heatmap of the exposure vectors",
    xlab="Signature",
    ylab="Sample",
    cexRow=0.8,
    cexCol=0.8
  )
}


#' Heatmap of exposure vectors from tumor genomes
#'
#' `tumorHeatmap` function computes the exposure vectors from an MPF file and
#' plots an heatmap, clustering the different tumor genomes by their
#' similarities in the exposure vectors.
#'
#'
#' @usage tumorHeatmap(mpfFilePath,
#' signatures,
#' signaturesType=signatureTypes$alexandrov32,
#' numBases=5,
#' trDir=TRUE,
#' enforceUniqueTrDir=TRUE,
#' refGenome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
#' transcriptAnno=TxDb.Hsapiens.UCSC.hg19.knownGene::
#' TxDb.Hsapiens.UCSC.hg19.knownGene,
#' verbose=FALSE,
#' plot=TRUE)
#'
#'
#'
#' @param mpfFilePath (Mandatory) The path of the MPF file from
#'  which to load the tumor genomes
#'
#' @param signatures (Mandatory) Either a path to the signature
#' file or a vector
#` of path to signatures file (in the case of Shiraishi signatures). Also, it
#' can be dircetly an object containing the signatures. The object can be a
#' list of vectors, data frames or matrices. Vectors are
#' used for Alexandrov
#' signatures, data frames or matrices for Shiraishi signatures.
#'
#' @param signaturesType (Optional) A string representing the signature type.
#' A list of the allowd signature types is available as the object
#'  `signatureTypes`. Use alexandrov2 for Alexandrov V2 signatures,
#' alandrov32 for Alexandrov V3.2 ones and shiraishi for Shiraishi signatures.
#' DEFAULT: signatureTypes$alexandrov32
#'
#' @param numBases (Optional) Total number of bases (mutated base and flanking
#'  bases) to be used for sequence patterns. Must be odd. Default: 5
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
#' @param refGenome (Mandatory) The reference genome (BSgenome)
#' needed to extract sequence patterns. Default: BSgenome object for hg19.
#'
#' @param transcriptAnno (Optional) Transcript annotation (TxDb object) used to
#' determine the transcription direction. This is required only if trDir is
#' TRUE. Default: TxDb object for hg19.
#'
#' @param verbose (Optional) Print information about reading and processing the
#' mutation data. Default: TRUE
#'
#' @param plot (Optional) Output the heatmap of the computed exposures
#' vectors. Default: TRUE
#'
#' @details
#' Many of the parameters of this function are also described in the function
#' `readGenomesFromMPF()` and `decomposeTumorGenomes()` from the package
#' `decompTumor2Sig` as those functions and package are at the
#' base of this function and many parameteres are passed directly to those
#' functions. A more throughout explanation of the parameters and the functions
#' can be found in their documentation, liked in the see also section.
#' @seealso [decompTumor2Sig]
#' The links to the documentation of the mentioned functions from the package
#' `decompTumor2Sig`
#' \link[decompTumor2Sig]{readGenomesFromMPF}
#' \link[decompTumor2Sig]{decomposeTumorGenomes}
#'
#' @return The exposure vector is returned and a heatmap of their similarity is
#' plotted
#'
#' @examples
#'##read breast cancer genomes from Nik-Zainal et al (PMID: 22608084)
#'gfile <- system.file("extdata", "Nik-Zainal_PMID_22608084-MPF.txt.gz",
#'                   package="decompTumor2Sig")
#'
#'##get the filenames with the shiraishi signatures
#'sigfiles <- system.file("extdata",
#'                       paste0("Nik-Zainal_PMID_22608084-pmsignature-sig",
#'                              1:4,".tsv"),
#'                       package="decompTumor2Sig")
#'##compute the exposure vectors and plot the heatmap
#'exposures <- tumorHeatmap(gfile, sigfiles,
#'                         signaturesType=signatureTypes$shiraishi)
#'
#'
#' @importFrom decompTumor2Sig isSignatureSet
#' @importFrom decompTumor2Sig readGenomesFromMPF
#' @importFrom decompTumor2Sig decomposeTumorGenomes
#' @importFrom decompTumor2Sig sameSignatureFormat
#' @importFrom stats heatmap
#' @importFrom utils read.table
#' @importFrom methods is
#' @export
tumorHeatmap <- function(mpfFilePath,
                         signatures,
                         signaturesType=signatureTypes$alexandrov32,
                         numBases=5,
                         trDir=TRUE,
                         enforceUniqueTrDir=TRUE,
                         refGenome=BSgenome.Hsapiens.UCSC.hg19::
                           BSgenome.Hsapiens.UCSC.hg19,
                         transcriptAnno=TxDb.Hsapiens.UCSC.hg19.knownGene::
                           TxDb.Hsapiens.UCSC.hg19.knownGene,
                         verbose=FALSE,
                         plot=TRUE) {
  signatures <-
    .readSignaturesFromFileOrObject(signatures, signaturesType)
  genomes <- .readGenome(mpfFilePath,
                        numBases,
                        signaturesType,
                        trDir,
                        refGenome,
                        transcriptAnno,
                        verbose)
  ##check that the genome and the signatures are compatible
  if (!sameSignatureFormat(signatures, genomes)) {
    stop(
      paste0(
        "The parameter used to load the genomes need to be the same",
        " as the ones used to generate the signatures. Please verify",
        " that this is the case"
      )
    )
  }
  
  exposures <- .computeExposures(genomes, signatures, verbose)
  .plotHeatmapAndDendogram(exposures)
  
  return(exposures)
}