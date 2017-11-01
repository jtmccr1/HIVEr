context("Testing functions that deal with transmission.")
set.seed(11)
# This ensures we don't get a warning from the logit functions due
# to not having a diverse enough freq1 set up in the attatched data

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
  sample_data<-x[[1]] # this is the data used in the last run
  model_out<-x[[2]]

  expect_equal(nrow(sample_data),4)
  expect_equal(length(unique(sample_data$pair_id)),2)

  expect_equal(unique(model_out$trial),c(1,2))
  expect_equal(nrow(model_out),8) # 2 variants, 2 samples, 2 trials


})
test_that("sampling pairings with resued donor",{
  expect_warning(com_sample_trans(small_community.comp,2,c("HS1595","MH0001")),regex = "Removed")

  w <- capture_warnings(com_sample_trans(small_community.comp,2,c("HS1595","HS1595","MH0000")))
  if(length(w)>0){# There are warnings are they all from glm
    expect_match(w, "glm.fit: fitted probabilities numerically 0 or 1 occurred", all = TRUE)
  }
  x<-com_sample_trans(small_community.comp,2,c("HS1595","HS1595","MH0000"),test=T)
  sample_data<-x[[1]]
  model_out<-x[[2]]
  expect_equal(length(which(sample_data$SPECID1=="HS1595")),4) # two variants for 2 pairs
  expect_equal(length(which(sample_data$SPECID1=="MH0000")),2) # two variants for 1 pairs

})

test_that("Sampling without donor restriction",{

  x<-sample_trans(small_community.comp,2,2,test=T)
  sample_data<-x[[1]]
  model_out<-x[[2]]

  expect_equal(nrow(sample_data),4)
  expect_equal(length(unique(sample_data$pair_id)),2)

  expect_equal(unique(model_out$trial),c(1,2))
  expect_equal(nrow(model_out),100)
})
