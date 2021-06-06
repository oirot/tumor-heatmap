#list of allowed signature types
#' A list of the available signatures
#' 
#' @export
signatureTypes <- list(alexandrov2 = "Alexandrov.2", alexandrov32 = "Alexandrov.32",
                       shiraishi = "Shiraishi")

getGenomeType <- function(signatureType){
  return(strsplit(signatureType, ".", fixed = TRUE)[[1]][1])
}

#TODO: add check for the format
#' @importFrom utils read.table
readAlexandrovV32Signatures <- function(filenames){
  signatures <- as.matrix(read.table(filenames, sep="\t", header = TRUE,
                                     row.names = 1))
  mutationIDs <- rownames(signatures)
  signatures <- split(signatures, col(signatures, as.factor = TRUE))
  signatures <- lapply(signatures, function(l){
    names(l) <- mutationIDs
    return(l)
  })
  return(signatures)
}

#'@importFrom decompTumor2Sig readAlexandrovSignatures
#'@importFrom decompTumor2Sig readShiraishiSignatures
#vector of functions to read the allowd signature types
readSignaturesFunctions <- c(decompTumor2Sig::readAlexandrovSignatures,
                             readAlexandrovV32Signatures,
                             decompTumor2Sig::readShiraishiSignatures)

#returns the correct function to read the signatures from the type
getReadSignatureFunction <- function(type){
  type <- match(type, signatureTypes) 
  return(readSignaturesFunctions[[type]])
}