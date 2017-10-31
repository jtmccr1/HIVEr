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
                     SPECID" , id_snv$SPECID[2]),"\n")
      return(id_snv$SPECID[2])
    }else{
      if(id_snv$gc_ul[1]<id_snv$gc_ul[2]){
        message(paste0("SPECID ", id_snv$SPECID[1],
                       "  being used for ENROLLID ", id_snv$ENROLLID[1],
                       " based on the time to the transmission date",
                       "even though this sample has a lower titer than " , id_snv$SPECID[2],"\n"))
      }
      return(id_snv$SPECID[1])
    }
  }else{
    stop("stopping neither 1 nor 2 specid found for this sample")
  }
}

#' Randomly sample a transmission data.frame based on a list
#' of donors
#'
#' Each run we will ramdomly choose one of the possible transmission pairs
#' from data with a donor ENROLLID provided in the list and evalute the probability
#' of transmission. Smoothing every 2%.
#' @param data data frame on which to sample. contains freq1 and freq2
#' @param runs The number of times to sample
#' @param SPECIDs A list of SPECIDs that can be used
#' @return a data frame with columns freq1, trial, and probability of transmission.
#' @examples
#'
#' d<-small_meta$collect[4]+1
#' get_close(small_meta,d,enrollid = "50001",case = "Donor",iSNV = T)
#' @export
#'
com_sample_trans<-function(data,runs,SPECIDs,test=F){

  # randomly sample n rows from a data frame
  sample_n <- function(df,n){
    return(df[sample(nrow(df),n),])
  }

  start.df <- data.frame(freq1=seq(0.02,1,0.02))
  # We will smooth every 2%
  model.df<-data.frame(freq1=rep(seq(0.02,1,0.02),runs),
                      trial=rep(1:runs,each=length(seq(0.02,1,0.02))),
                      prob=NA)
  #print(model.df)
  for (i in 1:runs){
    # We only want pairings with an eligable donor
    possible_pairs <- subset(data,SPECID1 %in% SPECIDs)
    #pick one for each ENROLLID

    pairings <- possible_pairs %>%
                dplyr::group_by(SPECID1) %>%
                dplyr::summarize(pair_id = sample(.data$pair_id,1))

    #print(pairings)
    sampled_data <- subset(data,pair_id %in% pairings$pair_id)
    #print(sampled_data)

    #smooths the data.
    logit<-glm(formula =found~freq1,family=binomial(logit),data=sampled_data) # Fit a logit model to the data
    #print(logit)
    # Get the predictions on using the start frequencies

    final.df<-dplyr::mutate(start.df,
                     prob=predict(logit,data.frame(freq1=freq1),
                                  type="response"),trial=i)
    #print(final.df)
    #print(model.df)
    #data.sampled$prob = logit$fitted.values
    #final.df=data.frame(freq1=data.sampled$freq1,prob=data.sampled$prob,trial=i) # Get the predictions on using the start frequencies

    ind<-which(model.df$trial==i)
    model.df$prob[ind]<-final.df$prob# add to the final output
  }
  if(!test){
  return(model.df)
  }
  if(test){
    return(list(sampled_data,model.df))
  }
}
