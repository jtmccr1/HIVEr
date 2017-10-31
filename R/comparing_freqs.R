
#' Fair base call comparisons.
#'
#' DeepSNV calls minority variants according to a plasmid control which means that often the reference bases is the
#' major allele in the sample (but not always). The DeepSNV script provides all alleles present at a variable site.
#' In this data set there is only ever 2 but there could be more. In comparing two samples it is likely that sample 1
#' contains a minority variant at some position that in sample 2 is characterized by only one allele. If the allele in
#'  sample 2 is the plasmid reference allele there will be no record of it in the data frame. This function adds that
#'  data so that it is clear the samples contain the same major allele (or different major alleles if that is the case).
#'  If both samples contain only the major allele and it is the same allele an empty data frame is returned as we are not
#'  interested in cases were the alleles are the same and there is no minor allele.
#'
#'
#'  If there is only one allele but we might expect two given the frequency then we warn the user. This can happen
#'  when one allele is removed because of poor quality. We do not adjust the frequency of the
#'  remaining alleles in this case. If the present allele matches the only allele in the other sample
#'  then the position is removed as the data suggests both samples have only one allele and
#'  the same allele (at that).
#'
#'
#'
#'
#' @param position A data frame containing one row for each base called at a loci for the two samples in question.
#' @return A data frame with either the missing allele updated (in the appropriate sample), the missing allele added as
#' a new entry in the data frame, or an empty data frame (in the case where the is only one allele and it is fixed in both samples)
#' @examples
#' # Here we add the reference base in the correct row.
#' position <- data.frame(ENROLLID1 = c(300294,300294),
#'                        ENROLLID2 = c(300293,300293),
#'                        mutation = c("M_C806A","M_C806C"),
#'                        freq1 = c(0,0),
#'                        freq2 = c(0.9,0.1),
#'                        chr = c("M","M"),
#'                        pos = c(806,806),
#'                        ref = c("C","C"),
#'                        var = c("A","C"),stringsAsFactors = F)
#' print(position)
#' equal_compare(position)
#'
#' # Here we add a row when the major alleles differ.
#' position <- data.frame(ENROLLID1 = c(300294),
#'                        ENROLLID2 = c(300293),
#'                        mutation = c("M_C806A"),
#'                        freq1 = c(0),
#'                        freq2 = c(0.99),
#'                        chr = c("M"),
#'                        pos = c(806),
#'                        ref = c("C"),
#'                        var = c("A"),stringsAsFactors = F)
#' print(position)
#' equal_compare(position)
#'
#' # Here we remove the entry as the samples are the same
#' position <- data.frame(ENROLLID1 = c(300294),
#'                        ENROLLID2 = c(300293),
#'                        mutation = c("M_C806A"),
#'                        freq1 = c(1),
#'                        freq2 = c(0.99),
#'                        chr = c("M"),
#'                        pos = c(806),
#'                        ref = c("C"),
#'                        var = c("A"),stringsAsFactors = F)
#' print(position)
#' equal_compare(position)
#' @export

equal_compare <- function(position){ # take in all the variants found in a sample pair at a given posistion and correct for differences in infering.Meaning the reference base was called here but only in one sample because in the other sample there wasn't a minor allele. Each position with a minor variant should have at least 2 variants,in both samples. after this is done. A minor varaint and the reference. If a position has novariants in either sample then it is excluded from this analysis. In this case only the reference base was present in both samples.  If only major variants are present then the major allele is fixed and there are no other minor variants to infer.

  pos_sum_1 <- sum(position$freq1) # What is the sum of all variants at this site in the first sample
  pos_sum_2 <-  sum(position$freq2) # What is the sum of all variants at this site in the second sample
  stopifnot(pos_sum_1+pos_sum_2>0)#"No variant at this position - something terrible happend somewhere"
  if((pos_sum_1> 0 & min(position$freq1)<0.95 & pos_sum_1==min(position$freq1)))
    warning(paste("There is only 1 allele in sample",position$SPECID1,
    "Its frequency is :",position$freq1), " it may be removed.\n")
  if((pos_sum_2> 0 & min(position$freq2)<0.95 & pos_sum_2==min(position$freq2)))
    warning(paste("There appears to only be one allele in sample",position$SPECID2,
     "Its frequency is :", position$freq2)," it may be removed.\n") # only a minor variant called

  if(pos_sum_1==0){ # there isn't a variant in the first sample.The reference is fixed there. the sum of nothing is 0. Try it.  sum() =0.  sum(c()) =0

    x<-which(position$ref==position$var) # Where are the infered calls
    stopifnot(length(x)<2) # should be 1 or 0.
    if (length(x)==1){ # The reference has been infered in the second sample but there are no varaints in the first so it was not infered
      position$freq1[x]<-1 # The reference base is fixed here
    } else if(length(x)==0){ # the reference was not infered in the second sample and no variants here. Add the reference now. There must have been another base fixed in the second sample.
      extra<-dplyr::mutate(position,var=ref,mutation=paste0(chr,"_",ref,pos,var),freq1=1,freq2=0) # adding a line with the reference base fixed in sample 1 but not sample 2
      position<-rbind(position,extra)
      #print("need to add a few here")
    }
    # Do it again if there are no varaints in the second sample.
  } else if(pos_sum_2==0){ # there isn't a variant in the second sample. The reference is fixed. also I told you so. (line 196)
    x<-which(position$ref==position$var)
    stopifnot(length(x)<2) # should be 1 or 0
    if (length(x)==1){ # The reference has been infered in the first sample but there are no varaints in the second so it was not infered
      position$freq2[x]<-1
    } else if(length(x)==0){ # the reference was not infered in the first sample and no variants here. Add the reference now
      extra<-dplyr::mutate(position,var=ref,mutation=paste0(chr,"_",ref,pos,var),freq1=0,freq2=1) # adding a line with the reference base fixed in sample 2 but not sample 1
      position<-rbind(position,extra)
    }
  }

  if(nrow(position)>1){ # if after all this there is only 1 row in the position then the variant base is fixed in both samples and not interesting to us. We don't look at all the sites that have the reference base only in both.
    return(position)
  }else{
    return(position[F,])
  }
}


#' Get and compare the frequencies of mutations in two samples.
#'
#' This function compares the frequencies of alleles found in two samples. It relies
#' on equal_compare to make fair comaparisons at each loci.
#'
#' @param pairs A vector of SPECID. These strings must be found in the SPECID column of the iSNV data frame
#' @param snv A data frame of isnv calls from deepSNV with a SPECID column
#' @return a data frame with a row for each allele that differs between the two samples. Each loci that differs
#' should have two rows (one for each allele). If there are no iSNV or none that differ then an empty data frame
#' is returned.
#' @examples
#' # There  one shared iSNV here and one fixed difference the rest are the same.
#' print(small_isnv)
#' get_freqs(c("HS1595","HS1563"),small_isnv)
#' @export
#' @seealso \link{equal_compare}
#' @importFrom magrittr %>%


get_freqs<-function(pairs,snv){
  # take in a data frame of pairs of Ids. Only 1 pair,
  # and a list of snv calls. and output the comparison between the 2.
  # each iSNV and its frequency in both samples
  snv<-subset(snv,SPECID %in% pairs,select=
                c(SPECID,mutation,chr,pos,ref,var,freq.var,season,pcr_result))

  # just need the snv in the samples we're looking at.
  if(nrow(snv)>0){ # There are mutations.
    # We only compare samples from the same season and strain so in each case the reference base in the same
    stopifnot(length(unique(snv$season))==1,length(unique(snv$pcr_result))==1)

    mut_table <- tidyr::spread(snv,SPECID,freq.var,fill=0) # a data frame with mutation down the first row and then frequency in either sample in the next 2.

    mut_table$SPECID1<-pairs[1] # add column with first sample ID
    mut_table$SPECID2<-pairs[2] # add column with second sample ID
    names(mut_table)[which(names(mut_table)==as.character(pairs[1]))] <-'freq1' # rename this column as the frequency in the first sample
    names(mut_table)[which(names(mut_table)==as.character(pairs[2]))] <-'freq2' # dido

    # This function can only be run on samples that qualified for variant identification.
    # If no variants were found in the sample then the SPECID will be missing from mut_table column and so
    # freq1 or freq2 will be missing since nothing was found we set that to 0 here and add the column.
    # equal compare will replace these cases with the reference at 1.
    if(!('freq1' %in% names(mut_table))){
      mut_table$freq1 <- 0
    }
    if(!('freq2' %in% names(mut_table))){
      mut_table$freq2 <- 0
    }
    mut_table <- dplyr::select(mut_table,mutation, chr, pos, ref, var, season, pcr_result, freq1,freq2, SPECID1, SPECID2)
    mut_table <- mut_table[order(mut_table$chr,mut_table$pos),]

    # fill in differences based on infered. In this case we are left with only sites that are polymorphic in one or between the 2 samples
    all.freq <- mut_table %>%
      dplyr::group_by(chr,pos) %>%
      dplyr::do(equal_compare(.))

    if(nrow(all.freq)>0){ # Some differences exists
      return(all.freq)
    }
    else{ # some snv were initially present (before frequency filter) but after taking into account differences in inference the samples are the same.
      x<-data.frame(mutation=NA,freq1=NA,freq2=NA,SPECID1=NA,
                    SPECID2=NA,chr=NA,ref=NA,"pos"=NA,var=NA)
      return(x[F,])
    }
  }
  else{ # No variants found in either sample
    x<-data.frame(mutation=NA,freq1=NA,freq2=NA,SPECID1=NA,
                  SPECID2=NA,chr=NA,ref=NA,"pos"=NA,var=NA)
  return(x[F,])
  }
}

#' Get the genetic distance between two samples .
#'
#' This function compares the frequencies of alleles found in two samples and
#' calculates the L1 norm or manhattan distance
#'
#' @param pairs A vector of SPECID. These strings must be found in the SPECID column of the iSNV data frame
#' @param snv A data frame of isnv calls from deepSNV with a SPECID column
#' @return the L1 norm
#' @examples
#' # There  one shared iSNV here and one fixed difference the rest are the same.
#' print(small_isnv)
#' dist_tp(c("HS1595","HS1596"),small_isnv)
#' @export


dist_tp<-function(pairs,snv){
  data.df<-get_freqs(pairs = pairs,snv = snv)

  if(nrow(data.df)==0){
    # This is the case when there are no variants found in either sample
      d=0
    }else{
      y<-as.matrix(cbind(data.df$freq1,data.df$freq2))
      d=dist(t(y),method = "manhattan")
    }
  #print(paste(pairs,d,sep=":"))
  return(as.numeric(d))
  # This is the same as L1-norm I checked it by hand (1/18/17)
} # data frame of 1 pair



#' Polishing compared frequencies.
#'
#' This function polishes a frequency comparision data frame. The output
#' is a subset  of the initial data frame only containing loci that are
#' polymorphic in the specified sample.This function is needed because there are cases where a major allele
#' exists but there is not minor allele because of sequencing errors. In much of the analysis
#' we are only interested in polymophisms. This function ensures that is all we observe.
#'
#' @param df data frame with columns freq1 and freq2
#' @param relative either freq1 or freq2 . unquoted. only look at sites with frequency above 0 at freq1,freq2
#' @param min_freq the minumun frequency allowed. frequencies below will be set to 0. those above to 1

#' @return a data frame with polymophoric loci
#' @examples
#' get_freqs(c("HS1595","HS1563"),small_isnv)->small_dups
#' polish_freq(small_dups,0.02,freq1)
#'
#' @export

polish_freq<-function(df,relative,min_freq=0){
  rel_quo<-enquo(relative)
  min_freq_quo<-enquo(min_freq)
  mess = paste0("Requiring polymophic at least in ", quo_name(rel_quo))
  message(mess)

  # We apply a freqeuncy cut off. it will be at least 0
  frequency_cut<-rlang::expr(UQ(rel_quo)>UQ(min_freq_quo))
  df<-dplyr::filter(df,UQ(frequency_cut))

  # We want polymorphic sites in the given sample. Due to frequency oddities the major allels may not
  # be exactly 1-min_freq i.e. we don't set any frequencies by hand. minor frequencies are set by deepsnv
  # and this includes an estimated error rate. Major frequencies are set by allele counts.

  frequency_cut_minor<-rlang::expr(UQ(rel_quo)<1)



  # How many alleles have frequencies in the choosen sample above the min and below 1 at each loci
  counts<- df %>% dplyr::group_by(chr,pos,SPECID1,SPECID2) %>%
  dplyr::summarize(alleles=length(which(UQ(frequency_cut) & UQ(frequency_cut_minor ))))

  if(max(counts$alleles>2)) stop(paste0("More that two alleles found at a position in a sample"))
  df<-merge(df,counts)

  return(subset(df,alleles==2,select=-c(alleles)))
  }
