

small_isnv$HOUSE_ID <- 000000
get_freqs(c("HS1595","HS1563"),small_isnv)->small_dups
small_dups.comp<-polish_freq(small_dups,freq1,0.02)

extra<-small_dups.comp
extra$SPECID1<-"MH0000"
small_dups.comp$freq2[1]<-0

small_dups.comp<-rbind(small_dups.comp,extra)
small_dups.comp<-rbind(small_dups.comp,small_dups.comp)
small_dups.comp<-rbind(small_dups.comp,small_dups.comp)

small_dups.comp$pair_id<-rep(1:(nrow(small_dups.comp)/2),each=2)

small_dups.comp$freq1<-runif(n = nrow(small_dups.comp),0,1)

small_dups.comp$found<-small_dups.comp$freq2>0.02
small_community.comp<-small_dups.comp
small_dups.comp$SPECID2<-"MH0001"
small_community.comp<-rbind(small_community.comp,small_dups.comp)

small_community.comp$pair_id<-rep(1:(nrow(small_community.comp)/2),each=2)

devtools::use_data(small_community.comp,overwrite = T)
