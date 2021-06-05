readSignatures <- function(signaturePath, signatureType){
  signatures <- getReadSignatureFunction(signatureType)(signaturePath)
  return(signatures)
}
