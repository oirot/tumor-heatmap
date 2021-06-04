#TODO: add check for the format
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

readAlexandrovV2Signatures <- function(filename){
  
  signatures <- readAlexandrovSignatures(filename)
  return(signatures)
}

readShiraishiSignatures <- function(filenames){
  signatures <- readShiraishiSignatures(files = filenames) 
  return(signatures)  
}

#list of allowed signature types
signatureTypes <- list(alex2 = "AlexV2", alex3 = "AlexV3.2", shi = "shi")
#vector of functions to read the allowd signature types
readSignaturesFunctions <- c(readAlexandrovV2Signatures,
                             readAlexandrovV32Signatures,
                             readShiraishiSignatures)

read.signature <- function(signaturePath, signatureType){
  type <- match(signatureType, signatureTypes) 
  signatures <- readSignaturesFunctions[[type]](signaturePath)
  return(signatures)
}
