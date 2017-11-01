meta<-read.csv("../../Analysis/HIVE/data/reference/all_meta.sequence_success.csv",colClasses = c('onset'='Date','collect'='Date'),stringsAsFactors = F) # Read in the meta data
meta_isnv<-read.csv("../../Analysis/HIVE/data/meta_snv_qual.csv",colClasses = c('onset'='Date','collect'='Date'),stringsAsFactors = F) # Read in the meta data
subset(meta,SPECID %in% c("HS1595","HS1563") | ENROLLID %in% c("50001","50068","50069") | HOUSE_ID %in% c(5017,5001,5228,3195,1068))->small_meta

 #merge and set NAs to 0 and note it iwll be diff than data.
# house 3195 doesn't have a onset date. It use to through an error. it is here to ensure it doesn't
devtools::use_data(small_meta,overwrite = T)


# These are removed in meta_snv_qual "MH3040" "MH8689" "MH8691" "HS1568" "HS1569" "MH8693"
