
#' Picking one sample per person
#'
#' In the analysis it is sometimes useful to just have one sample per ENROLLID.
#' This is particularily true when we are simply interest in data specific to
#' the indivudal and not the sample. This function reduces a meta file that
#' may contain many SPECID/ENROLLID season pcr_result and reduces it so that just one SPECID
#' exists for each ENROLLID. We choose the best SPECID defined by
#'  the sequenced sample
#'  the sample that qualified for snv calling
#'  the sample with highest titer.
#'  If no samples were seqeunced we choose randomly.
#'
#' @param meta_df A data frame containing one row for each SPECID in the data set.
#' @return A data frame with one row for each ENROLLID in the data frame.
#' @examples
#'print(small_meta)
#'only_one(small_meta)
#'
#'
only_one_helper<-function(x){
  # if there is only one sample then we are done
  if(nrow(x)==1){
    return(x)
  }else if(nrow(x)==2){
    # There are two samples here with this person and pcr result
    # pick the ones that we sequenced
    pick<-which(x$sequenced==T)
    if(length(pick)==1){
      #if there is only one then we're done
      return(x[pick,])
    }
    else if(length(pick)==0){
      # we didn't sequence any of these
      pick<-1
      return(x[pick,]) # take the first one
    }
    else if(length(pick)==2){
      # we sequenced both
      # See how many qualified for snv calling
      pick<-which(x$snv_qualified==T)
      # Either 0 or 2 did
      if(length(pick)==0 | length(pick)==2){
        # take the one with a higher titer.
        pick<-which(x$gc_ul==max(x$gc_ul))
      }
      return(x[pick,])
    }
  }else{
    stop("more than 2 samples for this person?")
  }
}

only_one<-function(meta_df){ #  one per person
  out <- meta_df %>% dplyr::group_by(ENROLLID,pcr_result,season) %>%
    dplyr::do(only_one_helper(.))
  return(out)
}


# ==============================================================================
# Transmission pairs

#' Validating transmission pairs
#'
#' Here we take a data frame where each row is a putative transmission pair
#' and determine which ones are valid. Valid pairs may only have the same
#' onset date if they are the index cases of the house. Such cases get a double
#' tag as well. Otherwise we only except pairs that do not jump a middle case,
#' and there can not be 2 possible donors (same onset date) for a recipient.
#' That is thrown out.
#'
#' @param tp A data frame containing one row for each possible transmission pair
#' @return A data frame with one row for each possible transmission pair
#'
#' @describeIn getting_tp

finding_valid_helper<-function(x){
  # for each possible recipient in the house we will find which donor is valid
  # Get the difference in onset date
  x=dplyr::mutate(x,diff_onset=onset2-onset1)
  # valid pairs at this point must have a difference on onset date >0.
  valid_onsets=x$diff_onset[x$diff_onset>0]

  # the recipient is sick on the house index case and onset1==onset2 and there is no other possible donor date - these are all sick on the first day and need to be handeled in both directions

  if(identical(x$onset1,x$onset2) & all(x$onset2 == min_onset)){
    x$valid<-T
    x$double<-T
  }
  else if(length(valid_onsets)>0){
    # There are some valid donors
    # If there is only 1 valid donor here with onset date closest to that of recipient (i.e donors not sick on the same day)
    if(length(which(x$diff_onset==min(valid_onsets)))==1){

      x$valid[x$diff_onset==min(valid_onsets)]=T
    }
  }
  # otherwise we leave them as false

  x<-subset(x,select=-c(diff_onset)) # get rid of this column
  return(x)
}

finding_valid<-function(tp){
  # This function will take in a household and determine whether or not each pair
  # is it is a  valid transmission event - Note we have not subset on those with useful sequence data yet.
  # This the index date for the house it is used later
  min_onset<-min(tp$onset1)

  if(nrow(tp)==1){ # there is only one pair here so it's valid
    tp$valid=T

    # there is only one pair and they are sick on the same day
    # so this is valid but we will need to look at transmission both ways

    if(tp$onset1==tp$onset2){
      tp$double=T
    }
    return(tp)
  }
  else{
    # There are multiple possible pairs here
    # All that is left is to remove the pairs that skip an individual - arc over the top.
    # For this we will look at each recipeint and validate the pair that has the minimum diff in onset date
    # but whose difference is >0
    tp %>% dplyr::group_by(ENROLLID2) %>%
      do(finding_valid_helper(.))->tp_multi
    return(tp_multi)
  }
}






#' Finding transmission pairs
#'
#' Here we take a data frame and look for all cases where two individuals are
#' sick within a week of each other. This does not compute wether or not the
#' pairs match our criteria yet
#'
#' @param meta_one A data frame containing one row for each ENROLLID in the data set.
#' @return A data frame with one row for each possible transmission pair
#' @examples
#' one_meta<-only_one(small_meta)
#' getting_tp(one_meta)
#'
#'

getting_tp_helper<-function(df){
  house<-unique(df$HOUSE_ID)
  stopifnot(length(house)==1,length(unique(df$pcr_result))==1) # Verify only 1 house here and only 1 strain
  if(length(unique(df$ENROLLID))!=length(unique(df$SPECID))) stop("Please supply a data frame with only one entry/ENROLLID")
  pairs<-data.frame(ENROLLID1=NA,ENROLLID2=NA,onset1=NA,onset2=NA,transmission=NA,
                    sequenced1=NA,sequenced2=NA,gc_ul1=NA,gc_ul2=NA,
                    snv_qualified1=NA,snv_qualified2=NA,
                    valid=NA,double=NA)
  # verify more than 1 person in the household
  if(length(unique(df$ENROLLID))>1){

    # Order the data frame based on onset date. Earliest is first
    df<-df[order(df$onset,df$ENROLLID,decreasing = F),]
    # loop through and get all possible transmission pairs and make 1
    # row for each possible transmission pair I know this is not tidy but
    # that's ok for now
    for(i in 1:(nrow(df)-1)){
      ENROLLID1=df$ENROLLID[i]
      onset1=df$onset[i]
      sequenced1=df$sequenced[i]
      gc_ul1=df$gc_ul[i]
      snv_qualified1=df$snv_qualified[i]
      for(j in (i+1):nrow(df)){
        ENROLLID2=df$ENROLLID[j]
        onset2=df$onset[j]
        sequenced2=df$sequenced[j]
        gc_ul2=df$gc_ul[j]
        snv_qualified2=df$snv_qualified[j]
        transmission=onset2-1
        if((onset2-onset1)<=7){ #if onset is withing a week of each other

          out<-data.frame(ENROLLID1=ENROLLID1,ENROLLID2=ENROLLID2,
                          onset1=onset1,onset2=onset2,
                          transmission=transmission,
                          sequenced1=sequenced1,sequenced2=sequenced2,
                          gc_ul1=gc_ul1,gc_ul2=gc_ul2,
                          snv_qualified1=snv_qualified1,snv_qualified2=snv_qualified2,
                          valid=F,double=F) # all start as not valid
          pairs<-rbind(subset(pairs,!(is.na(ENROLLID1))),out)# We don't want that first line with the NA this just adds the out data frame with the transmission pair if it meets the onset date
        }
      }
    }
    valid_pairs<-finding_valid(pairs)
    valid_pairs<- dplyr::mutate(valid_pairs,# Just adding columns to summarize the pairs.
                         sequenced_pair=sequenced1==sequenced2 & sequenced2==T,
                         titer_pair= sequenced1==sequenced2 & sequenced2==T & gc_ul1>1e3 & gc_ul2>1e3,
                         snv_qualified_pair = snv_qualified1==T & snv_qualified2==T)
    return(valid_pairs)
    #return(pairs)
  } else {
    return(pairs[F,]) # We need to return a data.frame for each entry. This is empty
  }
}
getting_tp<-function(meta_one){ # The data frame is all for 1 house. The only criteria here is that the samples are from the same house, and onset is within a week of eachother. The data.frames are ordered such that onset2 is on the same day or later than onset1
  result <- meta_one %>% dplyr::group_by(HOUSE_ID,pcr_result,season) %>%
    dplyr::do(getting_tp_helper(.))
  return(result)
}



