
#' Mutation rate and Ne from spectrum
#'
#' This function determines the Ne and mutation rate assuming all infections begin as a clone
#' and then mutations are built upon this base. The probability of a site going from 0 to
#' x given mu, Ne, and t generations given in Rouzine 2001 JVI as
#'
#' \deqn{g(f,t) = \frac{2 \mu N_e}{f} e^{- \frac{2N_ef}{t}}}
#'
#' This model makes some assumptions.
#' 1) Every infection began with every potential loci with a mutation at frequency 0 - we only model minor alleles
#' 2) The population develops under neutral processes
#' 3) Because we intrate numerically we assume anything present at less than 0.001 (0.1%)
#' is essentially at 0.
#' 4) There are 13133 sites that can contain mutations - the number of nucleotides in our concatenated 2014-2015 sample
#'
#' @param data The data frame containing iSNV calls. It must have the following columns:
#' freq.var,DPI-days post infection, and only contain 1 allele/ site. Sites not present will be assumed
#' to not have a variant present.
#' @param mu The mutation rate per nucleotide per cellular infectious cycle
#' @param Ne The effective population size
#' @param gen_time The generation time in hours. Default to 6
#' @param neg Boolean return the negative log likelihood not. Default =FALSE
#' @param acc An accuracy data frame
#' @return The log likihood for the data given the parameters and the model
#'
#' @export
mut_model_LL<-function(data,mu,Ne,gen_time=6,neg=FALSE,acc=accuracy_stringent){

   if(max(data$freq.var>0.5)){
    data<-data[data$freq.var<0.5,]
    warning("Removing major alleles")
  }
  #First we get the contribution of each iSNV that was found
  found<- data %>% dplyr::rowwise()%>%
    dplyr::mutate(prob = g_ft(mu=mu,Ne=Ne,t=as.numeric(DPI,units="days")*(24/gen_time),
                              x=freq.var,gc_ul = gc_ul,acc=acc))

    found_LL = sum(log(found$prob)) # These are all the found snv

  # Now we'll take into account all the points that aren't there
  data <- dplyr::mutate(data,gc_ul.acc = 10^(floor(log10(gc_ul)))) # Get the gc_ul in terms of accruacy- ie round down

  data$gc_ul.acc[data$gc_ul.acc>1e5]<-1e5

  min.count = data %>% dplyr::group_by(SPECID,DPI,gc_ul.acc) %>%
    dplyr::summarize(count = length(mutation)) # How many iSNV in this samples

  min.count = dplyr::mutate(min.count,fixed = 13133-count)
  # flu genome size (concatenated coding region for NY - 2014-2015) How many non mutated points
  # Just the sites without iSNV
# We are interested in DPI and gc_ul.acc as that determines the Likelihood of no variants. after that we
# just need to know the number of sites that match these characteristics
  sum_table.missed = min.count %>% dplyr::group_by(DPI,gc_ul.acc) %>%
    dplyr::summarize(total_fixed = sum(fixed)) # how many fixed points in these samples

  sum_table.missed$DPI= as.numeric(sum_table.missed$DPI,unit="days")
  #print(sum_table.missed)

  sum_table.missed = sum_table.missed %>%
                        dplyr::rowwise()%>%
                         dplyr::mutate(
                            # likelihood for each DPI and acc of not observing a mutation
                            like = not_detected(mu = mu,Ne= Ne,t = DPI*(24/gen_time)
                                                ,gc_ul = gc_ul.acc),
                            # Add up for each DPI and acc
                            LL = total_fixed*log(like))

  missed_LL<-sum(sum_table.missed$LL) # Get total likelihood for missed

  if(neg==T){
    return(-(found_LL+missed_LL))
  }else{
    return((found_LL+missed_LL)) # return total likelihood for found and missed samples.
  }
}

#' Likelihood of frequency starting at 0
#'
#' This function determines the Ne and mutation rate assuming all infections begin as a clone
#' and then mutations are built upon this base. The probability of a site going from 0 to
#' x given mu, Ne, and t generations given in Rouzine 2001 JVI as
#'
#' \deqn{g(f,t) = \frac{2 \mu N_e}{f} e^{- \frac{2N_ef}{t}}}
#'
#' The output is adjusted to account for the sensitivity. This isn't important in fitting,
#' but it is nice in testing to know the prob of being found+ prob of not being found =1
#' @param x The measured frequency - must be between 0 and 1
#' @param mu The mutation rate per nucleotide per cellular infectious cycle
#' @param Ne The effective population size
#' @param t Generations
#' @param gc_ul titer default to NULL use if interested in adjusting for sensitivity
#' @param acc sensitivity data frame default is NULL use if interested in adjusting for sensitivity
#' @export



g_ft<-function(x,mu,Ne,t,gc_ul=NULL,acc=NULL){ # x is the frequency
  prob <- (2*mu*Ne/x) * exp(-1*2*Ne*x/t)
  if(!is.null(gc_ul)){
    acc_gc=10^(floor(log10(gc_ul)))
    if(acc_gc>1e5){
      acc_gc <- 1e5
    }
    if(x>=0.02 & x<0.05){
      sense= acc$sensitivity[which(acc$gc_ul==acc_gc & acc$freq==0.02)]
    }else if(x>=0.05 & x<0.1){
      sense= acc$sensitivity[which(acc$gc_ul==acc_gc & acc$freq==0.05)]
    }
    else {
      sense=1
    }
    return( prob*sense )
  }

  else{
    return(prob)
  }
}



#' Mutational model- non polmorphic
#'
#' What's the probability the site still is at frequency 0?
#' This assumes anything below 0.001 (0.1%) is at 0. It yields 1- integral g_ft from 0.001 to 1
#' @param mu The mutation rate per nucleotide per cellular infectious cycle
#' @param Ne The effective population size
#' @param t Generations

p0 <- function(mu,Ne,t){ # 1 - polymorphic - Polymorphic being defined as between 0.001 and 1
  poly <- integrate(g_ft,lower=0.001,upper=1,mu=mu,Ne=Ne,t=t)
  1-poly$value
}

#' Probability of not detecting a polymorphism
#' This is the probability the site is less than 0.001 (assumed 0 in this model) or
#' not found by our sequencing method.
#'
#' @param mu The mutation rate per nucleotide per cellular infectious cycle
#' @param Ne The effective population size
#' @param t The number of generations
#' @param gc_ul The sample titer
#' @param acc  dataframe holding the sensitivity
#' @return The probability for the data given the parameters and the model

not_detected<-function(mu,Ne,t,gc_ul,acc=accuracy_stringent){
  # p0
  naught = p0(mu,Ne,t) # probability it is at 0
  # below cut
  below = integrate(g_ft,lower=0.001,upper=0.02,
                    mu=mu,Ne=Ne,t=t) # below detection
  # missed
  acc_gc=10^(floor(log10(gc_ul)))
  if(acc_gc>1e5){
    acc_gc <- 1e5
  } # This will round the gc down to the nearest log as we have always done to be conservative
  uncert_term=c()
  f=c(0.02,0.05,0.10)
  for(i in 1:(length(f)-1)){
    uncert=1-acc$sensitivity[which(acc$gc_ul==acc_gc & acc$freq==f[i])]
    # The prob the variant is missed because itis between f[i] and f[i+1] given the sample size
    integrand = integrate(g_ft,lower=f[i],upper=f[i+1],
                          mu=mu,Ne=Ne,t=t)
    uncert_term[i]=integrand$value*uncert # the probability it is present * the probability of not seeing it
  }
  missed = sum(uncert_term)

  return(naught+below$value+missed) # the probability of being present and missed or not present
}





