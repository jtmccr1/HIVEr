context("Testing functions that compare isnv calls between samples.")

position.1 <- data.frame(ENROLLID1 = c(300294,300294),
                       ENROLLID2 = c(300293,300293),
                       mutation = c("M_C806A","M_C806C"),
                       freq1 = c(0,0),
                       freq2 = c(0.9,0.1),
                       chr = c("M","M"),
                       pos = c(806,806),
                       ref = c("C","C"),
                       var = c("A","C"),stringsAsFactors = F)

result.1 <- data.frame(ENROLLID1 = c(300294,300294),
                         ENROLLID2 = c(300293,300293),
                         mutation = c("M_C806A","M_C806C"),
                         freq1 = c(0,1),
                         freq2 = c(0.9,0.1),
                         chr = c("M","M"),
                         pos = c(806,806),
                         ref = c("C","C"),
                         var = c("A","C"),stringsAsFactors = F)

test_that("Reference allele is added in place", {
  expect_equal(equal_compare(position.1), result.1)
})

position.2 <- data.frame(ENROLLID1 = c(300294),
                       ENROLLID2 = c(300293),
                       mutation = c("M_C806A"),
                       freq1 = c(0),
                       freq2 = c(1),
                       chr = c("M"),
                       pos = c(806),
                       ref = c("C"),
                       var = c("A"),stringsAsFactors = F)

result.2 <- data.frame(ENROLLID1 = c(300294,300294),
                       ENROLLID2 = c(300293,300293),
                       mutation = c("M_C806A","M_C806C"),
                       freq1 = c(0,1),
                       freq2 = c(1,0),
                       chr = c("M","M"),
                       pos = c(806,806),
                       ref = c("C","C"),
                       var = c("A","C"),stringsAsFactors = F)

test_that("Reference allele is in a new row", {
  expect_equal(equal_compare(position.2), result.2)
})

position.3 <- data.frame(ENROLLID1 = c(300294),
                       ENROLLID2 = c(300293),
                       mutation = c("M_C806A"),
                       freq1 = c(1),
                       freq2 = c(1),
                       chr = c("M"),
                       pos = c(806),
                       ref = c("C"),
                       var = c("A"),stringsAsFactors = F)

test_that("Remove sites that are fixed and equal", {
  expect_equal(nrow(equal_compare(position.3)), 0)
})

position.4 <- data.frame(ENROLLID1 = c(300294),
                         ENROLLID2 = c(300293),
                         mutation = c("M_C806A"),
                         freq1 = c(0),
                         freq2 = c(0.5),
                         chr = c("M"),
                         pos = c(806),
                         ref = c("C"),
                         var = c("A"),stringsAsFactors = F)
position.5 <- data.frame(ENROLLID1 = c(300294),
                         ENROLLID2 = c(300293),
                         mutation = c("M_C806A"),
                         freq1 = c(0.5),
                         freq2 = c(0),
                         chr = c("M"),
                         pos = c(806),
                         ref = c("C"),
                         var = c("A"),stringsAsFactors = F)
position.6 <- data.frame(ENROLLID1 = c(300294),
                         ENROLLID2 = c(300293),
                         mutation = c("M_C806A"),
                         freq1 = c(0),
                         freq2 = c(0),
                         chr = c("M"),
                         pos = c(806),
                         ref = c("C"),
                         var = c("A"),stringsAsFactors = F)
test_that("Error if something is wrong with the data", {
  expect_warning(equal_compare(position.4))
  expect_warning(equal_compare(position.5))
  expect_error(equal_compare(position.6))


})

#=============================================
# Getting frequency test
#
freq_out = data.frame(mutation =c("PB2_G58A","PB2_G58G","PB2_G1099A","PB2_G1099G"),
                      chr = rep("PB2",times=4),
                      pos = rep(c(58,1099),each=2),
                      ref = rep("G",times=4),
                      var = rep(c("A","G"),times=2),
                      season = rep("2014-2015",times=4),
                      pcr_result = rep("A/H3N2",times=4),
                      freq1 = c(1,0,0.04,0.96),
                      freq2 = c(0,1,0.1,0.9),
                      SPECID1 = rep("HS1595",times=4),
                      SPECID2 = rep("HS1563",times = 4),stringsAsFactors = F)

test_that("Getting frequencies",{
  expect_equal(get_freqs(c("HS1595","HS1563"),small_isnv),freq_out)
})

test_that("Getting distance",{
  expect_equal(dist_tp(c("HS1595","HS1563"),small_isnv),2.12)
})


small_freqs<-get_freqs(c("HS1595","HS1563"),small_isnv)
test_that("polishes frequency data to polymorphic sites",{
  expect_equal(nrow(setdiff(polish_freq(small_freqs,0.02,freq1),polish_freq(small_freqs,0.02,freq2))),0)
  expect_equal(nrow(setdiff(freq_out[3:4,],polish_freq(small_freqs,min_freq = 0.02,freq2))),0)
})
