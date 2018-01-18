require(magrittr)
require(tidyverse)

A<-read_csv("./data-raw/HS1381_A.removed.ref.sum.filtered.AA.csv")
B<-read_csv("./data-raw/HS1381_B.removed.ref.sum.filtered.AA.csv")
covs <-read_rbind(c("./data-raw/HS1381_A.removed.cov.csv","./data-raw/HS1381_B.removed.cov.csv"))

A<-subset(A,chr=="PB1" &
            (pos==1705|pos==530 |pos==389 |
               pos==1798 |pos==2032|pos==1617| pos==931))

B<-subset(B,chr=="PB1" &
         (pos==1705|pos==530 |pos==389 |
            pos==1798 |pos==2032|pos==1617| pos==931))


#The real variants are PB1_G530A,PB1_A1705G, & PB1_G1798G
variants<-rbind(A,B)

LAURING_ID_LOOKUP<-data_frame(Id=c(unique(variants$Id),"HS0000"))

LAURING_ID_LOOKUP<- LAURING_ID_LOOKUP %>% tidyr::separate(Id,into=c("LAURING_ID","dup"),
                                                          sep="_",remove=F,fill="left")


variants<-dplyr::left_join(variants,LAURING_ID_LOOKUP,by="Id")
variants<-dplyr::left_join(variants,covs,by=c("Id","chr","pos"="chr.pos"))

variants$gc_ul<-1.4e4
variants<-dplyr::select(variants,-dplyr::starts_with("X"))

variants<-dplyr::select(variants,-dplyr::starts_with("Unn"))
devtools::use_data(variants,overwrite = T)
