context("Testing functions that deal with transmission.")

d<-small_meta$collect[4]+1
small_meta$iSNV<-small_meta$iSNV<-round(runif(nrow(small_meta),1,10),0)
small_meta$iSNV[7]<-0
test_that("Correct SPECID is returned", {
  expect_equal(get_close(small_meta,d,enrollid = "50001",case = "Donor",iSNV = T),"HS1563")
  expect_equal(get_close(small_meta,d-6,enrollid = "50001",case = "Donor",iSNV = T),"HS1563")
  expect_equal(get_close(small_meta,d+6,enrollid = "50001",case = "Donor",iSNV = T),"MH8688")
  expect_warning(get_close(small_meta,d,enrollid = "50003",case = "Donor",iSNV = T))
})
