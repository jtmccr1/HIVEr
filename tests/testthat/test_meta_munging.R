context("Testing functions that wrangle the data frames")

meta_one <- only_one(small_meta)
test_that("Correctly handling mulitple SPECID", {
  expect_equal(meta_one$SPECID[meta_one$ENROLLID=="50001"],"HS1563") # Both good (gc test)
  expect_equal(meta_one$SPECID[meta_one$ENROLLID=="50068"],"MH8692") # One good
  expect_equal(meta_one$SPECID[meta_one$ENROLLID=="50069"],"HS1569") # Both bad
  expect_equal(meta_one$SPECID[meta_one$ENROLLID=="50993"],"HS1595") # Only one
})

tp<-getting_tp(meta_one)
just_one<-subset(tp,HOUSE_ID==5017)
# There are 3 valid on index day and 3 removed because of muliple donors
mult<-subset(tp,HOUSE_ID==5001)

jump<-subset(tp,HOUSE_ID==5228)
test_that("Correctly handl transmission pairs",{
  expect_equal(c(just_one$ENROLLID1,just_one$ENROLLID2),c("50069","50068"))
  expect_true(just_one$valid)

  expect_equal(mult$valid[mult$ENROLLID2=="50002"],c(F,F,F))
  expect_equal(sum(mult$valid[mult$ENROLLID2!="50002"]),3)
  expect_true(all(tp$ENROLLID1!=tp$ENROLLID2))

  expect_true(tp$valid[tp$ENROLLID1=="50969" & tp$ENROLLID2=="50967"])
  expect_true(tp$valid[tp$ENROLLID1=="50967" & tp$ENROLLID2=="50968"])
  expect_false(tp$valid[tp$ENROLLID1=="50969" & tp$ENROLLID2=="50968"])
  expect_true(all(tp$onset1<=tp$onset2))

  expect_true(1068 %in% tp$HOUSE_ID) # these were lost at one point

})
