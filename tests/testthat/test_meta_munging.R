context("Testing functions that wrangle the data frames")

meta_one <- only_one(small_meta)
test_that("Correctly handling mulitple SPECID", {
  expect_equal(meta_one$SPECID[meta_one$ENROLLID=="50001"],"HS1563")
  expect_equal(meta_one$SPECID[meta_one$ENROLLID=="50068"],"MH8692")
  expect_equal(meta_one$SPECID[meta_one$ENROLLID=="50069"],"HS1569")
  expect_equal(meta_one$SPECID[meta_one$ENROLLID=="50993"],"HS1595")

})
