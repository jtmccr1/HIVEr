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

intra <- get_double(small_meta,2)
test_that("Identifies intahost pairs",{
  expect_equal(unique(intra$ENROLLID), c("50001","50003","50004","50068","50069","50968","50969"))
  expect_equal(nrow(intra),7*2)
})

longform_pairs(tp)->long_tp

test_that("longform test",{
  expect_equal(nrow(tp)*2,nrow(long_tp))
  expect_true(all(tp$ENROLLID1 %in% long_tp$ENROLLID))
  expect_true(all(tp$ENROLLID2 %in% long_tp$ENROLLID))

})


#==========================================================
# Short form
#
get_double(small_meta)->small_doub
small_doub<-dplyr::select(small_doub,ENROLLID,pcr_result,HOUSE_ID,SPECID,onset,collect,vaccination_status,DPI,season,gc_ul,sequenced,home_collected,snv_qualified)
columns = c("SPECID","onset","collect","vaccination_status",
            "DPI","gc_ul","sequenced","home_collected","snv_qualified")
columns_eror = c("SPECID","onset","collect","vaccination_status",
                 "DPI","sequenced","home_collected","snv_qualified")
short_doub <-short_pairs(small_doub,columns,ENROLLID,HOUSE_ID,pcr_result,season)
test_that("short form test",{
  expect_error(short_pairs(small_doub,columns_eror,ENROLLID,HOUSE_ID,pcr_result,season))
  expect_error(short_pairs(small_doub,columns_eror,pcr_result,season))
  expect_true(all(short_doub$onset1<=short_doub$onset2))
  expect_true(all(short_doub$SPECID1 <= short_doub$SPECID2))


})

