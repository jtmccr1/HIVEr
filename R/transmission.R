#' Finding SPECID for transmission pairs
#'
#' Identifies the SPECID collected nearest a given date, and on the correct side of
#' the date when specified. The criteria is
#'
#' 1) Sample closest to transmission that are on the "right side" of tranmission
#' (before transmission for donor if possible)
#'
#' 2) Titer is the tie breaker when applicable.
#' 3) If no iSNV were found in a donor sample
#'
#' @param meta data frame with meta data con
#' @param date The estimated trasmission date
#' @param enrollid The enrollid
#' @param case c("donor","recipient")
#' @param iSNV boolan if TRUE then we take the sample that has iSNV if the
#' "best fit" doesn't. This is the historic use of the function. Recall that we
#' do not have samples before transmission in many cases anyway.
#'
#' @return The SPECID that best fits the critea
#' @examples
#'
#' d<-small_meta$collect[4]+1
#' get_close(small_meta,d,enrollid = "50001",case = "Donor",iSNV = T)
#' @export

get_close<-function(meta,date,enrollid,case,iSNV=T){
  id_snv=subset(meta,ENROLLID==enrollid & snv_qualified==T)
  if(nrow(id_snv)==1){ # it there is only one sample from this individual
    return(id_snv$SPECID)
  }
  else if(nrow(id_snv)==2){ # there are two samples for this individual
    date<-enquo(date)
    id_snv=dplyr::mutate(id_snv,dtrans=.data$collect-!!date)
    # For donors we want the largest dtrans possible less than 0 or the smallest
    # possitive dtrans. For the recipient we want the smallest possitive dtrans.
    # It should be impossible for the recipient to have a negative dtrans but ENROLLID
    # 50538 is a recipient and does have a sample taken prior to onset.
    #We'll treat this like any other case. For now.

    if(all(id_snv$dtrans<0)){
      # all are less than 0 so this puts the sample with least negative dtrans on top.
      # titer is the tie braker
      id_snv=id_snv[order(id_snv$dtrans,id_snv$gc_ul,decreasing=T),]
    }else{
      # There is 1 or no dtrans less than 0 so this puts the sample with the
      # lowest dtrans (either the negative 1 or the smallest positive 1) on top.
      # titer is the tie braker
      id_snv=id_snv[order(id_snv$dtrans,-id_snv$gc_ul,decreasing=F),]
    }
    if(id_snv$iSNV[1]==0 & id_snv$iSNV[2]>0 & iSNV==T){
      # finally if there were no polymorphisms in the first sample and
      # there were a few in the second then we take the second
      warning(paste0(case," SPECID ", id_snv$SPECID[2], "  being used for ENROLLID ",
                     id_snv$ENROLLID[1], " because no iSNV were found in the prefered
                     SPECID" , id_snv$SPECID[2]))
      return(id_snv$SPECID[2])
    }else{
      if(id_snv$gc_ul[1]<id_snv$gc_ul[2]){
        message(paste0("SPECID ", id_snv$SPECID[1],
                       "  being used for ENROLLID ", id_snv$ENROLLID[1],
                       " based on the time to the transmission date",
                       "even though this sample has a lower titer than " , id_snv$SPECID[2]))
      }
      return(id_snv$SPECID[1])
    }
  }else{
    stop("stopping neither 1 nor 2 specid found for this sample")
  }
}

