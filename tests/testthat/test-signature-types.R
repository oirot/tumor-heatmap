test_that("Get genome type from signature type", {
    expect_equal(.getGenomeType(signatureTypes$alexandrov2), "Alexandrov")
    expect_equal(.getGenomeType(signatureTypes$alexandrov32), "Alexandrov")
    expect_equal(.getGenomeType(signatureTypes$shiraishi), "Shiraishi")
})

test_that("Returning correct signature reading function", {
    expect_equal(
        .getReadSignatureFunction(signatureTypes$alexandrov2),
        decompTumor2Sig::readAlexandrovSignatures
    )
    expect_equal(
        .getReadSignatureFunction(signatureTypes$shiraishi),
        decompTumor2Sig::readShiraishiSignatures
    )
    expect_equal(
        .getReadSignatureFunction(signatureTypes$alexandrov32),
        .readAlexandrovV32Signatures
    )
})
