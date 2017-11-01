qual<-read.csv("../../Analysis/HIVE/data/processed/qual.snv.csv",stringsAsFactors = F)
subset(qual,SPECID %in% c("HS1595","HS1563"))->small_isnv
subset(small_isnv,chr=="PB2"& ref=="G")->small_isnv
small_isnv%>% select(HOUSE_ID,ENROLLID,SPECID,mutation,chr,pos,ref,var,freq.var,season,pcr_result)->small_isnv
# add fixed difference
extra_fixed<-small_isnv[2,]
extra_fixed$pos <- 58
extra_fixed$mutation <- "PB2_G58A"


small_isnv<-rbind(small_isnv,extra_fixed)

small_isnv$freq.var[small_isnv$SPECID=="HS1595" & small_isnv$mutation=="PB2_G1099A"] <- 0.04
small_isnv$freq.var[small_isnv$SPECID=="HS1595" & small_isnv$mutation=="PB2_G1099G"] <- 0.96


small_isnv$freq.var[small_isnv$SPECID=="HS1563" & small_isnv$mutation=="PB2_G1099A"] <- 0.1
small_isnv$freq.var[small_isnv$SPECID=="HS1563" & small_isnv$mutation=="PB2_G1099G"] <- 0.9

devtools::use_data(small_isnv,overwrite = T)
