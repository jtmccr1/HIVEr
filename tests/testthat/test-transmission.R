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


test_that("sampling pairings correctly",{
  x<-com_sample_trans(small_community.comp,2,c("HS1595","MH0000"),test=T)
  sample_data<-x[[1]]
  model_out<-x[[2]]

  expect_equal(nrow(sample_data),4)
  expect_equal(length(unique(sample_data$pair_id)),2)

  expect_equal(unique(model_out$trial),c(1,2))
  expect_equal(nrow(model_out),100)




})
