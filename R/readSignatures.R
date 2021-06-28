## read the signatures according to their signature type
.readSignatures <- function(signaturePath, type) {
    signatures <- .getReadSignatureFunction(type)(signaturePath)
    return(signatures)
}
