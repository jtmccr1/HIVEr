#' Is this a whole number
#'
#' Function that tets whether input is a whole number and
#' returns T/F. This is used creating the pmf of discrete
#' distributions. This is taken from the example in as.integer
#' documentation included in base.
#'
#'
#' @param x A number
#' @return Logical is the number an integer
#'
#' @examples
#' is_wholenumber(2)
#'

is_wholenumber <-function(x, tol = .Machine$double.eps^0.5){
   abs(x - round(x)) < tol
}

#' PMF of a zero truncated Poisson distribution
#'
#' PMF of zero truncated Poisson distribution at a point x.
#'
#'
#' @param x A number or vector of numbers
#' @param lambda The mean of the Poisson distribution that is truncated
#' not the mean of the truncated distribution.
#' @return Logical is the number an integer
#'
#' @examples
#' dzpois(1,3)
#'
#' # If Nb and l are vectors the output is a matrix of the form
#' # |---------Nb_j---------|
#' # |                      |
#' # |                      |
#' # l_i                   l_i
#' # |                      |
#' # |                      |
#' # |                      |
#' dzpois(1:2,1:3)
#' @export
dzpois<-function(x,lambda){ # the pmf of the zero truncated poisson distribution
  if(x<=0 | is_wholenumber(x)==F){
    return(0)
  }else{
    (exp(-1*lambda)*lambda^x)/(factorial(x)*(1-exp(-1*lambda)))

  }
}
dzpois<-Vectorize(dzpois,vectorize.args = c("x"))

#' Probability of drawing an allele n times in n draws
#'
#' Probability of drawing an allele n times in n draws. This is a helper
#' function used in transmission modeling.
#'
#' @param p The probability (frequency of allele)
#' @param n The number of draws. or a vector of draw values
#' @return The probability of only drawing an allele at the given frequency.
#' Or a vector of such probabilities
#'
#' @examples
#' p_all(0.5,1)
#' p_all(0.5,c(1,2,3,4,5,6))
p_all<-function(p,n){ # probability all success - only finding the variant at frequency p in n draws
  p^n
}
p_all<-Vectorize(p_all,vectorize.args="n")


#' The probability of lambda for a loci
#'
#' The probability of observing the data given lambda and the presence absence
#' model.In this fit we take the minority frequency to be  correct
# and set the major frequency to 1-minority. This accounts for
# the fact that frequencies are related by the expression :
# minor allele + major allele + errror =1.
#Here we make the major allele frequency = major allele + error.
#The error is always small. if it exceeds 1% then we through an error here.

#'
#' @param data a data frame of one donor polymorphic loci. It must have
#' columns chr,pos,freq1, and found
#' @param l a vector of lambda values to test
#' @param Nb_max The maximum bottleneck size to test
#' @param model The model we are using must be either "PA" or "BetaBin"
#' @param threshold limit of variant calling detection
#' @param acc a data frame with accuracy metrics
#'
#' @return a tibble with columns for each lambda and the probability of
#' observing the data for that site.
#'
#' @examples
#' polish_freq(small_dups,freq1,0.02)->x
#' x$found=c(T,T)
#' math_fit(x,1:2)
#'
#' x$found=c(F,T)
#' math_fit(x,1:2,100,"PA")
#'
#' @export
math_fit=function(data,l,Nb_max,model,threshold,acc){
  # this  calculation is for each position in the genome.
  stopifnot(length(unique(data$chr))==1, length(unique(data$pos))==1)
  stopifnot(model %in% c("PA","BetaBin"))
  Nb <- 1:Nb_max
  Nb_given_l<-matrix(dzpois(Nb,l),ncol = length(l),byrow=T)
  # This gives a 100 by 1000 matrix.
  # |-------lambda_j-------|
  # |                      |
  # |                      |
  #  Nb_i                   Nb_i
  # |                      |
  # |                      |
  # |                      |

  # we now need to determine the probability of observing the data for each bottleneck
  # and each model.
  if(model=="PA"){
    #  In this fit we take the minority frequency to be  correct
    # and set the major frequency to 1-minority. This accounts for
    # the fact that frequencies are related by the expression :
    # minor allele + major allele + errror =1.
    #Here we make the major allele frequency = major allele + error.
    #The error is always small. if it exceeds 1% then we through an error here.

    if(1-sum(data$freq1)>0.01){
      stop("The sum of the frequencies is less than 99%")
    }
    data$freq1[data$freq1==max(data$freq1)]<-1-min(data$freq1)


    found<-data[data$found==T,] # only alleles that were transmitted
    Nb<-1:Nb_max # Here are the bottlenecks
    # we need the probability of getting a bottleneck of size Nb give lambda


    if(nrow(found)==0 | nrow(found)>2){
      stop(paste0("No variant transmitted for this site or",
                  "there are more than 2 variants here"))
    }else if(nrow(found)==1){ # one variant found here. All successes
      prob<-p_all(p=found$freq1,n=Nb)
      # this is a vector of probabilities for each
      # n prob[i]= the probability of only getting that
      # variant in Nb[i] (i.e. draws)
    }else if(nrow(found)==2){
    # if at least on of each allele was transmitted

      first_var<-p_all(p=found$freq1[1],n=Nb) # all this one
      second_var<-p_all(p=found$freq1[2],n=Nb) # all the other one
      one_each<-1-(first_var+second_var)
      # at least one of each -
      # This is a vector as above since R adds and subtracts the
      # elements of the vectors as expected
      prob<-one_each
    }
  }else if(model=="BetaBin"){
    if(nrow(data[data$freq1>0.5,])>0){
      data<-data[df$freq1<0.5,]
      warning("The beta binomials model only uses minor alleles. Subsetting the data now.")
    }
    # 2 is the recipient
    v_r = data$freq2
    v_d = data$freq1
    gc_ul = data$gc_ul2
    threshold = 0.02

    prob = L.Nb.beta(v_r,v_d,Nb,gc_ul,threshold,acc)
  }


  conditional_prob<-matrix(prob,nrow=1) %*%  Nb_given_l
  # This is matrix muliplication - it results in P(l)=P(Data|Nb)P(Nb|l)
  # summed over all Nb for each l
  #
  #                        |-------lambda_i-------|
  #                        |                      |
  #                        |                      |
  # [...P(Data|Nb_i)...] * Nb_i                   Nb_i  = [.....P(lambda_i).....]
  #                        |                      |
  #                        |                      |
  #                        |                      |

  conditional_prob<-as.vector(conditional_prob)
  return(tibble(lambda=l,prob=conditional_prob))
}

#' Fit the transmission model
#'
#' This is wrapper around math_fit that fits the presence absence model to
#' each pair present in the data set.
#' @param data a data frame or tibble it will be split by chr pos and pair_id
#' @param l a vector of lambda values to test
#' @param Nb_max The maximum bottleneck size to test
#' @param model PA or BetaBin
#' @param threshold limit of variant calling detection
#' @param acc a data frame with accuracy metrics
#' @param ... other columns to group by in final output.
#' @return a tibble with columns pair_id,lambda,LL (log likelihood), and pair_id if desired.
#' @export

trans_fit<-function(data,l,Nb_max,model,threshold,acc,...){
  group <- rlang::quos(...,lambda)
  probs<-data %>% dplyr::group_by(chr,pos,pair_id) %>%
    dplyr::do(math_fit(.data,l,Nb_max,model,threshold,acc))
  # For each genomic position in question
  LL.df<-probs %>% dplyr::group_by(!!!group) %>%
    dplyr::summarize(LL=sum(log(prob)))
  # Get the  log likelihood of the each lambda for this pair - we can sum accros pairs later.
  return(LL.df)
}

#' Summarize model likelihoods
#'
#' This function summarizes the likelihood from the model_fit functions.
#' returning the max lambda and 95% confidence interval. The confidence interval
#' assumes a smooth curve.
#'
#' @param data data frame or tibble of model fit data.
#' @return a tibble with columns
#' lambda : best fit
#' lower_95 : lower 95% bound (lambda)
#' upper_95 : upper 95% bound (lambda)
#' Nb : mean Nb of a zero truncated poisson given the lambda
#'
#' @examples
#' small_trans %>% pa_fit(1:2,100) %>% model_summary()
#' @export

model_summary<-function(data){
  Nb<-data$lambda[which(data$LL==max(data$LL))] # Get the  max lambda
  good_range<-subset(data,LL> (max(LL)-1.92)) # get the bottlenecks that fall in this region the 95% confidence intereval
  lower<-good_range$lambda[1]
  upper <- good_range$lambda[nrow(good_range)]
  return(tibble(lambda=Nb,lower_95=lower,
                upper_95=upper,mean_Nb = lambda/(1-exp(-1*lambda))))
}


#' Simulate transmission Presence absence model
#'
#' Simulate the transmission of alleles using the
#' presence absence model. Essentially this takes
#' the frequency of the variant allele, and bottleneck
#' size and simulates the found column. It also selects a bottleneck size
#' based on lambda and a zero truncated Poisson. If simulatin the whole data
#' set each pair should be run separately as to not use the same bottleneck
#'
#' @param data a data frame with with chr,pos,freq1, and  pair_id columns
#' @param lambda The lambda of the zero truncated Poisson
#'
#' @return data frame the same as data but with a simulated found column
#'
#' @examples
#' pa_sim(small_trans,1.3)
#' @export
pa_sim<-function(data,lambda){
  pair<-unique(data$pair_id)
  if(length(pair)>1){
    warning(paste0("Running on ",length(pair)," pairs. All will have the same bottleneck."))
  }
  Nb<-rzpois(1,lambda)
  out<-data %>% dplyr::group_by(chr,pos,pair_id) %>%
    dplyr::do(pa_sim_helper(.,Nb))
  return(out)
}


#' PA sim helper
#'
#' Helper function to simulate presence absence data. This takes in one loci.
#' the data frame should have 2 rows (one for each allele.)
#' Simulate the transmission of alleles using the
#' presence absence model. Essentially this takes
#' the frequency of the variant allele, and bottleneck
#' size and simulates the found column.
#'
#' @param data a data frame with with freq1 and pair_id columns and 2 rows
#' @param Nb the bottleneck size
#'
#' @return data frame the same as data but with a simulated found column
#'

pa_sim_helper<-function(data,Nb){
  # data refers to the 2 mutations at this position, bottlenecks-
  #  a data frame with the bottle_neck,
  #model - the name of the column in the bottlenecks that
  # contains the bottleneck size we want to use.
  if(nrow(data)!=2){
    stop(print("There are not 2 mutations at this point."))
  }
  pair<-unique(data$pair_id)

  if(length(pair)!=1){
    stop(print("There should only be one pair id at this point."))
  }
  freq_success<-1-min(data$freq1)
  # This is the major allele frequency - calculated the same as in the model

  success<-rbinom(1,Nb,freq_success)
  # the number of success in n trials.
  #  How many of the major alleles in the draw
  #
  data$found<-F
  if(success==Nb){ # only found major variants
    data$found[data$freq1==max(data$freq1)]=T
  } else if(success==0){ # only the minor
    data$found[data$freq1==min(data$freq1)]=T
  } else if(success>0 & success<Nb){   # both were found
    data$found=T
  } else{
    stop("Error!")
  }
  return(data)
}

#' Randomly sample a zero truncated Poisson
#'
#' Randomly sample a zero truncated Poisson distribution given
#' the lambda of the underlying Poisson (that is of course truncated).
#' @param n number samples
#' @param lambda The mean of the untruncated Poisson
#' @return a random sample
#' @examples
#' rzpois(10,2.3)
#' @export
rzpois<-function(n,lambda){
  #from http://giocc.com/zero_truncated_poisson_sampling_algorithm.html
  out = vector(mode="double",length=n)
  for(i in 1:n){
    k = 1
    t = exp(-lambda) / (1 - exp(-lambda)) * lambda
    s = t
    u = runif(1)
    while(s < u){
      k = k+1
      t = t*(lambda / k)
      s = s+t
    }
    out[i] <- k
  }
  return(out)
}

#' Simulate probability of transmission  based on frequency
#'
#' @param data data frame of transmission data
#' @param runs how many simulations to run
#' @param lambda What is the lamda of the zero truncated Poisson to use
#' @param FUN what function to use to simulate the found column pa_sim or betabin_sim unquoted
#' @param threshold - optional to pass to betabinomial simulator
#' @param acc optional to pass to betabinomial simulation
#' @param ... columns to group for simulations. These will get the same
#' bottleneck size within each replication. It should probably only ever be
#' pair_id.(ie each person gets the same bottleneck within a run)
#'
#' @return a data frame of simulated logit fits to frequency in donor vs.
#' probability of transmission
#' @examples
#' simulations(small_trans,10,3.2,pa_sim)
#' simulations(small_trans,10,3.2,betabin_sim,0.02,accuracy_stringent)
#'
#'
#' @export
simulations<-function(data,runs,lambda,FUN,threshold=NULL,acc=NULL,...){
  # iSNV data, how many iterations,
  # what is lambda value and what function is used to simulate the data.

  group_by<-rlang::quos(...)
  rows_needed<-nrow(data)
  model.df<-tibble(freq1=rep(NA,rows_needed*runs),
                   trial=rep(1:runs,each=rows_needed),
                   prob=NA)
  pairs<-length(unique(data$pair_id))
    for (i in 1:runs){

    if(is.null(threshold)){
    trial.df<-data %>% dplyr::group_by(!!!group_by)%>%
      dplyr::do(FUN(.,lambda))
    }else{
      trial.df<-data %>% dplyr::group_by(!!!group_by)%>%
        dplyr::do(FUN(.,lambda,threshold,acc))
    }
    logit<-glm(formula =found~freq1,family=binomial(logit),data=trial.df) # Fit a logit model to the data
    trial.df$prob<-logit$fitted.values
    trial.df$trial<-i
    ind<-which(model.df$trial==i)
    model.df$prob[ind]<-trial.df$prob# add to the final output
    model.df$freq1[ind]<-trial.df$freq1
  }
  return(model.df)
}
#=============================================================================
#                           BETA BINOMIAL FUNCTIONS
#=============================================================================
#

#' @describeIn L.Nb.beta likelihood the variant is  found in recipient


like_found.beta<-function(v_r,v_d,Nb,gc_ul,threshold,acc=accuracy_stringent){
  acc_gc=10^(floor(log10(gc_ul)))
  if(acc_gc>1e5){
    acc_gc <- 1e5
  }
  # If the frequency is below our threshold then the probability of being found is 0
  if(v_r<threshold){
    return(0)
  }
  # This will round the gc down to the nearest log as discussed below.
  if(v_r>=0.02 & v_r<0.05){
    sense= acc$sensitivity[which(acc$gc_ul==acc_gc & acc$freq==0.02)]
  }else if(v_r>=0.05 & v_r<0.1){
    sense= acc$sensitivity[which(acc$gc_ul==acc_gc & acc$freq==0.05)]
  }
  else {
    sense=1
  }
  prob = c()
  if(v_r<=(1-threshold)){
    for(k in 0:Nb){
      prob[k+1]=dbeta(x=v_r,shape1 = k,shape2=Nb-k)*dbinom(x=k,size=Nb,prob = v_d)*sense
    }
    prob=sum(prob)
  }else if(v_r>(1-threshold)){
    # it is fixed - whats the probability the other allele was not found
    lost_allele_freq = 1-v_d
    # The for loop over k and sum is in the like_lost.beta.uncert function
    prob = like_lost.beta.uncert(v_r=0,v_d=lost_allele_freq,Nb,gc_ul,threshold,acc=acc)
    }
    return(prob)
}

#' @describeIn L.Nb.beta likelihood the variant is not found in recipient
like_lost.beta.uncert<-function(v_r,v_d,Nb,gc_ul,threshold,acc=accuracy_stringent){ # sum over k
  stopifnot(v_r<threshold) # ensure not found in the recipient
  acc_gc=10^(floor(log10(gc_ul)))
  if(acc_gc>1e5){
    acc_gc <- 1e5
  }
  # This will round the gc down to the nearest log as discussed below.
  prob=c()
  for(k in 0:Nb){
    uncert_term=c()
    f=c(0.02,0.05,0.10)
    for(i in 1:(length(f)-1)){
      uncert=1-acc$sensitivity[which(acc$gc_ul==acc_gc & acc$freq==f[i])]
      # The prob the variant is missed because itis between f[i] and f[i+1] given
      # the sample size
      uncert_term[i]=(pbeta(q = f[i+1],shape1 = k,shape2 = Nb-k)-
                        pbeta(q = f[i],shape1 = k,shape2 = Nb-k))*
        dbinom(x=k,size=Nb,prob = v_d)*uncert
    }
    #probability the variant is below the cut off or present but missed
    prob[k+1]=pbeta(q = threshold,shape1 = k,shape2 = Nb-k)*dbinom(x=k,size=Nb,prob = v_d)+sum(uncert_term)
  }
  sum(prob)
}

#' Betabinomial Likelihood functions
#'
#' What is the probability of observing the data given a bottleneck
#' size under the beta binomial model. This is a wrapper around the two likelihood
#' functions
#' @param v_r allele frequency in the recipient
#' @param v_d allele frequency in the donor
#' @param Nb Bottleneck size (may be a vector)
#' @param gc_ul titer in recipeint sample
#' @param threshold detection level threshold
#' @param acc a data frame with accrucacy data must contain, columns
#'  freq,sensitivity, and gc_ul
#'
#' @return The likelihood of observing the data given the bottleneck size.
#' @export
L.Nb.beta<-function(v_r,v_d,Nb,gc_ul,threshold,acc=accuracy_stringent){
  if(v_r>=threshold){
        like=like_found.beta(v_r=v_r,v_d=v_d,Nb=Nb,gc_ul=gc_ul,threshold=threshold,acc=acc)
  }else if(v_r<threshold){ # Not found
      like=like_lost.beta.uncert(v_r=v_r,v_d=v_d,Nb=Nb,gc_ul=gc_ul,threshold=threshold,acc=acc)
    }
  return(like)
}
L.Nb.beta<-Vectorize(L.Nb.beta,vectorize.args = "Nb")


#' Simulate transmission Beta binomial model
#'
#' Simulate the transmission of alleles using the
#' betabinomial model. Essentially this takes
#' the frequency of the variant allele, and bottleneck
#' size and simulates the found column. It also selects a bottleneck size
#' based on lambda and a zero truncated Poisson. If simulation the whole data
#' set each pair should be run separately as to not use the same bottleneck
#'
#' @param data a data frame with with chr,pos,freq1, and  pair_id columns
#' @param lambda The lambda of the zero truncated Poisson
#' @param threshold limit of variant calling detection
#' @param acc a data frame with accuracy metrics
#'
#' @return data frame the same as data but with a simulated found column
#'
#' @examples
#' betabin_sim(small_trans,1.2,0.02,accuracy_stringent)
#' @export
betabin_sim<-function(data,lambda,threshold,acc=accuracy_stringent){

  pair<-unique(data$pair_id)
  if(length(pair)>1){
    warning(paste0("Running on ",length(pair)," pairs. All will have the same bottleneck."))
  }

  Nb<-rzpois(1,lambda)
  out<-data %>% dplyr::group_by(chr,pos,pair_id) %>%
    dplyr::do(betabin_sim_helper(.,Nb,threshold,acc))
  return(out)



}


#' beta bin sim helper
#'
#' Helper function to simulate beta binomial data. This takes in one loci.
#' the data frame should have 2 rows (one for each allele.)
#'  Essentially this takes
#' the frequency of the variant allele, and bottleneck
#' size and simulates the found column.
#'
#' @param data a data frame with with freq1 and pair_id columns and 2 rows
#' @param Nb the bottleneck size
#' @param threshold limit of variant calling detection
#' @param acc a data frame with accuracy metrics
#'
#' @return data frame the same as data but with a simulated found column
#'
betabin_sim_helper<-function(data,Nb,threshold,acc){
  # data refers to the 2 mutations at this position, bottlenecks-
  #  a data frame with the bottle_neck,
  #model - the name of the column in the bottlenecks that
  # contains the bottleneck size we want to use.
  if(nrow(data)!=2){
    stop(print("There are not 2 mutations at this point."))
  }
  pair<-unique(data$pair_id)

  if(length(pair)!=1){
    stop(print("There should only be one pair id at this point."))
  }

  minor<-min(data$freq1) # This is the minor allele frequency
  major<-1-min(data$freq1)

    # includes uncertainty -
    # with v_r = 0 or 1 these are probabilities otherwise they are densities.
    only_minor= L.Nb.beta(v_r=1,v_d=minor,Nb=Nb,gc_ul=unique(data$gc_ul2),threshold=threshold,acc)
      # includes uncertainty - The likelihood the allele was lost
    only_major = L.Nb.beta(v_r=1,v_d=major,Nb=Nb,gc_ul=unique(data$gc_ul2),threshold=threshold,acc) # includes uncertainty - The likelihood the allele was lost

  both = 1-(only_minor+only_major)

  # flip the coin
  pick = runif(1,0,1)
  data$found=F

  # 0  only_minor         only_major         both                        1
  # |------------------|---------------|---------------------------------|

  if(pick<=only_minor){
    # onlny the minor was found
    data$found[data$freq1==min(data$freq1)]=T
  }else if(pick>only_minor & pick<(only_minor+only_major)){
    # only the major was found
    data$found[data$freq1==max(data$freq1)]=T
  }else if(pick>=(only_minor+only_major)){
    # both were found
    data$found=T
  }
  return(data)
}

