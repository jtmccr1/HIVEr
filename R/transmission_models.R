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
#' math_fit(x,1:2)
#'
#' @export
math_fit=function(data,l=seq(0.01,10,0.01),Nb_max=100){
  # this  calculation is for each position in the genome.
  stopifnot(length(unique(data$chr))==1, length(unique(data$pos))==1)
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
  Nb_given_l<-matrix(dzpois(Nb,l),ncol = length(l),byrow=T)
  # This gives a 100 by 1000 matrix.
  # |-------lambda_j-------|
  # |                      |
  # |                      |
  #  Nb_i                   Nb_i
  # |                      |
  # |                      |
  # |                      |

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

#' Fit presence absence model
#'
#' This is wrapper around math_fit that fits the presence absence model to
#' each pair present in the data set.
#' @param data a data frame or tibble it will be split by chr pos and pair_id
#' @param l a vector of lambda values to test
#' @param Nb_max The maximum bottleneck size to test
#' @param ... other columns to group by in final output.
#' @return a tibble with columns pair_id,lambda,LL (log likelihood), and pair_id if desired.
#' @export

pa_fit<-function(data,l,Nb_max,...){
  group <- rlang::quos(...,lambda)
  probs<-data %>% dplyr::group_by(chr,pos,pair_id) %>%
    dplyr::do(math_fit(.data,l,Nb_max))
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
    warning(paste0("Runing on ",length(pair)," pairs. All will have the same bottleneck."))
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
#' @examples
#' pa_sim_helper(small_trans[1:2,],3)

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
#' @param FUN what function to use to simulate the found column
#' @parma ... columns to group for simulations. These will get the same
#' bottleneck size within each replication. It should probably only ever be
#' pair_id.(ie each person gets the same bottleneck within a run)
#'
#' @return
#'
#'
#' @export
simulations<-function(data,runs,lambda,FUN,...){
  # iSNV data, how many iterations,
  # what is lambda value and what function is used to simulate the data.

  group_by<-rlang::quos(...)
  rows_needed<-nrow(data)
  model.df<-tibble(freq1=rep(NA,rows_needed*runs),
                   trial=rep(1:runs,each=rows_needed),
                   prob=NA)
  pairs<-length(unique(data$pair_id))
    for (i in 1:runs){


    trial.df<-data %>% dplyr::group_by(!!!group_by)%>%
      dplyr::do(FUN(.,lambda))
    logit<-glm(formula =found~freq1,family=binomial(logit),data=trial.df) # Fit a logit model to the data
    trial.df$prob<-logit$fitted.values
    trial.df$trial<-i
    ind<-which(model.df$trial==i)
    model.df$prob[ind]<-trial.df$prob# add to the final output
    model.df$freq1[ind]<-trial.df$freq1
  }
  return(model.df)
}
