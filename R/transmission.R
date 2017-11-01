#' Finding SPECID for transmission pairs
#'
#' Identifies the SPECID collected nearest a given date, and on the correct side of
#' the date when specified. The criteria is
#'
#' 1) Sample closest to transmission that are on the "right side" of tranmission
#' (before transmission for donor if possible)
#'
#' 2) Titer is the tie breaker when applicable.
#'
#' 3) If no iSNV were found in a donor sample and iSNV==T we take the sample with iSNV present.
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
#' get_close(small_meta,d,enrollid = "50001",case = "Donor",iSNV = TRUE)
#' @export


get_close<-function(meta,date,enrollid,case,iSNV=TRUE){
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

#' Randomly sample a transmission tibble based on a list
#' of donors
#'
#' Each run we will ramdomly choose one of the possible transmission pairs
#' from data with a donor ENROLLID provided in the list and evalute the probability
#' of transmission. Smoothing every 2%.
#' @param data data frame on which to sample. contains freq1 and freq2
#' @param runs The number of times to sample
#' @param SPECIDs A list of SPECIDs that can be used as donors
#' @param test A boolean switch. If activatve returns a list of 2 tibbles. The first is
#' the last sample_data tibble the second is the normal output
#' @return a data frame with columns freq1, trial, and probability of transmission.
#' @examples
#' print(small_community.comp)
#' com_sample_trans(small_community.comp,2,c("HS1595","MH0000"))
#' @export
#'
com_sample_trans<-function(data,runs,SPECIDs,test=F){

  # It may be the case that there isn't another community pair for a given
  # donor. We identify it here. remove it and warn the user.
  orig_SPECID<-SPECIDs

  SPECIDs<-SPECIDs[SPECIDs %in% data$SPECID1]

  if(length(orig_SPECID)>length(SPECIDs)){
    warning(paste0("Removed donor SPECID: ",orig_SPECID[!(orig_SPECID%in%SPECIDs)]," Not in data\n"))
  }

  # How many comparisons are made.
  comparisons <- tibble(SPECID1=rep(SPECIDs,times=runs),
                             trial=rep(1:runs,times=length(SPECIDs)),
                             pair_id=NA)


  # we want to set up the tibble to hold the results here.
  row_df <- data %>% dplyr::group_by(SPECID1)%>%
   dplyr::summarise(mutations = length(unique(mutation)))
  row_df<- row_df %>% dplyr::group_by(SPECID1)%>%
    dplyr::mutate(mult = length(which(SPECIDs ==SPECID1)),
                  rows=mult*mutations)

  rows_needed<-sum(row_df$rows)

  model.df<-tibble(freq1=rep(NA,rows_needed*runs),
                      trial=rep(1:runs,each=rows_needed),
                      prob=NA)

  for (i in 1:runs){
    # We only want pairings with an eligable donor
    possible_pairs <- subset(data,SPECID1 %in% SPECIDs)
    #pick one for each ENROLLID

    # This ensures we pick different recipients for each donor within a run
    pass_unique_test<-F
    while(!pass_unique_test){
    pairings<-tibble(SPECID1=SPECIDs)
    pairings <- pairings %>%
                dplyr::rowwise() %>%
                dplyr::mutate(pair_id = base::sample(possible_pairs$pair_id[possible_pairs$SPECID1==.data$SPECID1],1,replace=F))

    unique_test<-dplyr::ungroup(pairings) %>% dplyr::group_by(SPECID1) %>%
      dplyr::summarise(unique= length(pair_id)==length(unique(pair_id)))

    pass_unique_test<-all(unique_test$unique==T)

    }

    sampled_data <- subset(data,pair_id %in% pairings$pair_id)

        logit<-glm(formula =found~freq1,family=binomial(logit),data=sampled_data) # Fit a logit model to the data
    #print(logit)
    # Get the fitted values on using the start frequencies

    sampled_data$prob<-logit$fitted.values

        # final.df<-dplyr::mutate(start.df,
    #                  prob=predict(logit,tibble(freq1=freq1),
    #                               type="response"),trial=i)
    #print(final.df)
    #print(model.df)
    #data.sampled$prob = logit$fitted.values
    #final.df=tibble(freq1=data.sampled$freq1,prob=data.sampled$prob,trial=i) # Get the predictions on using the start frequencies

    ind<-which(model.df$trial==i)
    model.df$prob[ind]<-sampled_data$prob# add to the final output
    model.df$freq1[ind]<-sampled_data$freq1


    #record how many unique pairings tested.
    #
    comp_ind<-which(comparisons$trial==i)
    # We want the pair data to be in the same order as the initial coparison data frame
    pairings<-pairings[order(match(pairings$SPECID1,SPECIDs)),]
    comparisons$pair_id[comp_ind]<-pairings$pair_id


  }

  comparison_meta<-comparisons %>% dplyr::group_by(SPECID1)%>%
    dplyr::summarise(comparisons = length(unique(pair_id)))


  message(paste0("Number donor SPECID used : ",nrow(comparison_meta),"\n",
                 "Number unique pairs made : ",sum(comparison_meta$comparisons),"\n",
                 "Mean unique pairs / donor : ",mean(comparison_meta$comparisons),"\n",
                 "Mininum pairs / donor : ", min(comparison_meta$comparisons)," for SPECID ",
                 comparison_meta$SPECID1[comparison_meta$comparisons==min(comparison_meta$comparisons)],"\n"
                 ))

  if(!test){
  return(model.df)
  }
  if(test){
    return(list(sampled_data,model.df))
  }
}

#' Randomly sample a transmission data.frame without any
#' restriction on donor SPECID
#'
#' Each run we will ramdomly choose 'pairs' of the possible transmission pairs
#' from data  smoothing every 2%.
#' @param data data frame on which to sample. contains freq1 and freq2
#' @param runs The number of times to sample
#' @param pairs The number of pairs to include in each run.
#' @param test A boolean switch. If activatve returns a list of 2 tibbles. The first is
#' the last sample_data tibble the second is the normal output
#' @return a data frame with columns freq1, trial, and probability of transmission.
#' @examples
#'
#' sample_trans(small_community.comp,2,3)
#' @export
#'
sample_trans<-function(data,runs,pairs,test=F){
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
  possible_pairs <- unique(data$pair_id)

  for (i in 1:runs){
    #print(pairings)

    run_pairs<-sample(possible_pairs,pairs)
    sampled_data <- data[data$pair_id %in% run_pairs,]
    #smooths the data.
    logit<-glm(formula =found~freq1,family=binomial(logit),data=sampled_data) # Fit a logit model to the data
    #print(logit)
    # Get the predictions on using the start frequencies
    final.df<-dplyr::mutate(start.df,
                            prob=predict(logit,data.frame(freq1=freq1),
                                         type="response"),trial=i)
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
