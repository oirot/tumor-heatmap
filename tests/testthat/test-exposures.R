test_that("Heatmap with Shiarishi signatures", {
  #read breast cancer genomes from Nik-Zainal et al (PMID: 22608084) 
  gfile <- system.file("extdata", "Nik-Zainal_PMID_22608084-MPF.txt.gz", 
                       package="decompTumor2Sig")
                       
  # get the filenames with the shiraishi signatures 
  sigfiles <- system.file("extdata",
                          paste0("Nik-Zainal_PMID_22608084-pmsignature-sig",
                                 1:4,".tsv"),
                          package="decompTumor2Sig")
  
  # compute the exposure vectors and plot the heatmap
  exposures <- tumorHeatmap(gfile, sigfiles, 
                            signaturesType = signatureTypes$shiraishi)
  #switch off the graphical device
  dev.off()
  expect_is(exposures, "matrix")
})

test_that("Heatmap with Alexandrov V3.2 signatures", {
  #read breast cancer genomes from Nik-Zainal et al (PMID: 22608084) 
  gfile <- system.file("extdata", "Nik-Zainal_PMID_22608084-MPF.txt.gz", 
                       package="decompTumor2Sig")
  
  sigfile <- system.file("extdata", "COSMIC_v3.2_SBS_GRCh37.txt", 
                       package="tumorHeatmap")
  exposures <- tumorHeatmap(gfile, sigfile, 
                            signaturesType = signatureTypes$alexandrov32,
                            numBases = 3,
                            trDir = FALSE)
  #switch off the graphical device
  dev.off()
  expect_is(exposures, "matrix")
})

test_that("Wrong signatures", {
  #read breast cancer genomes from Nik-Zainal et al (PMID: 22608084) 
  gfile <- system.file("extdata", "Nik-Zainal_PMID_22608084-MPF.txt.gz", 
                       package="decompTumor2Sig")
  
  sigfiles <- system.file("extdata",
                          paste0("Nik-Zainal_PMID_22608084-pmsignature-sig",
                                 1:4,".tsv"),
                          package="decompTumor2Sig")
  expect_error(tumorHeatmap(gfile, sigfiles, 
                            signaturesType = signatureTypes$alexandrov32,
                            numBases = 3,
                            trDir = FALSE))
})

test_that("Wrong parameters", {
  #read breast cancer genomes from Nik-Zainal et al (PMID: 22608084) 
  gfile <- system.file("extdata", "Nik-Zainal_PMID_22608084-MPF.txt.gz", 
                       package="decompTumor2Sig")
  
  sigfile <- system.file("extdata", "COSMIC_v3.2_SBS_GRCh37.txt", 
                         package="tumorHeatmap")
  expect_error(tumorHeatmap(gfile, sigfile, 
                            signaturesType = signatureTypes$alexandrov32,
                            numBases = 3,
                            trDir = TRUE))
})