context("Testing processing functions")

# variants 7 positions 7 major variants 3 minor ones that are real the outcome should be 10
#

test_that("Testing quality function",{
  q<-quality(variants)
  expect_equal(nrow(q),10)
  expect_equal(nrow(subset(q,freq.var<0.5)),3)
  variants$gc_ul<-1.1e5
  w<-capture_warnings(quality(variants))
  expect_match(w,"were sequenced twice")
  qsmall<-quality(subset(variants,Id=="HS1381_A"))
  expect_equal(nrow(qsmall),12)
  variants$gc_ul<-1e2
  w<-capture_warnings(quality(variants))
  expect_match(w,"It will be removed")

})

test_that("Testing correct id",{
  cov_sample<-dplyr::tibble(Id=c("129.0","MH00000"))
  x<-correct_id(cov_sample,Id)
  expect_equal(x$Id,c("129","MH00000"))
})

test_that("Testing diversity cut",{
  diverse_sites(small_isnv,1,season,pcr_result,pos,chr)->x
  expect_equal(1099,unique(x$pos))
})

test_that("Testing monomorphic sites",{
  wacky_isnv<-small_isnv
  wacky_isnv$freq.var[wacky_isnv$mutation=="PB2_G1602A"]<-0.8

  monomorphic(wacky_isnv,SPECID,season,pcr_result,pos,chr)->x
  out<-unique(x$freq.var[x$mutation=="PB2_G1602A"])
  expect_equal(out,1)
  })
