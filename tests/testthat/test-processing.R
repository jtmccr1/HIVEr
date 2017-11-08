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
  expect_error(quality(variants))

})
