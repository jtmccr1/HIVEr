trans_freq<-read.csv("../../Analysis/HIVE/results/transmission_pairs_freq.poly.donor.csv",
                     stringsAsFactors = F)
subset(trans_freq,pair_id %in% c(1,4,22))->small_trans

devtools::use_data(small_trans,overwrite = T)
