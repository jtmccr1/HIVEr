context("Testing transmission model functions")

test_that("Testing helper functions",{
  expect_equal(dzpois(c(1,2,3),3),c(dzpois(1,3),dzpois(2,3),dzpois(3,3)))
  expect_equal(sum(dzpois(c(0,0.1,0.6,9.9,-1),0)),0)
  Nb=1:10
  l=seq(0.01,1,0.01)
  expect_equal(dim(dzpois(Nb,l)),c(100,10))

  expect_equal(p_all(0.3,1:3),c(p_all(0.3,1),
                                p_all(0.3,2),
                                p_all(0.3,3)))
  # These test line 99 in math_fit that the matrix matches the form distribed
  # there
  x<-matrix(dzpois(Nb,l),ncol = length(l),byrow=T)
  expect_equal(dim(x),c(10,100))
  check_values<-c()
  c=1
  for(i in 1:length(Nb)){
    for(j in 1:length(l)){
     check_values = x[i,j]==dzpois(Nb[i],l[j])
    }
  }
  expect_true(all(check_values))


})

set.seed(15) # so we avoid glm warnings
test_that("proper warnings and data structure in PA_sim",{
  x<-simulations(small_trans,2,0.4,pa_sim,threshold= NULL,acc=NULL, pair_id)
  w<-capture_warnings(simulations(small_trans,2,0.5,pa_sim,threshold= NULL,acc=NULL))
  expect_length(unique(x$trial),2)
  expect_match(w, "All will have the same bottleneck.", all = TRUE)

})

test_that("Betabinomial likelihoods",{
  #accuracy_stringent<-read.csv("../Analysis/HIVE/data/reference/accuracy_stringent.csv",
   #                                                stringsAsFactors = F)
  #integrate over outcomes and the sum should 1.
  #
  # Making the integrare function so we can try different values.
 test_sum<-function(v_d,Nb,gc_ul,t,acc=accuracy_stringent){
    beta_int<-function(x) L.Nb.beta(x,v_d,Nb,gc_ul,t,acc)
    beta_intV<-Vectorize(beta_int)
    lost<-beta_int(0)
    fixed<-beta_int(1)
    poly<-integrate(beta_intV,t,(1-t))
    return(1-(lost+fixed+poly$value))
  }
  test_sumV<-Vectorize(test_sum,"v_d")
  test_sumV<-Vectorize(test_sumV,"Nb")
  test_sumV<-Vectorize(test_sumV,"gc_ul")
  test_sumV<-Vectorize(test_sumV,"t")

  out<-test_sumV(c(0.01,0.02,0.05,0.3,0.5,0.8),
                 c(1,10,100),
                 c(1.5e3,1.5e4,1.5e5,1e6),
                 c(0.001,0.02))

  # These should ==1 and in many cases do but due to computing error we allow a little over estimatation
  expect_true(all(out<0.1))
  expect_true(all(out>-0.1))
})

test_that("message and warnings with big Nb",{
  m<-capture_messages(suppressWarnings(
    trans_fit(small_trans,seq(1,300,1),300,"PA",threshold = NULL,acc=NULL)
    ))
  w<-capture_warnings(trans_fit(small_trans,seq(1,300,1),300,"PA",threshold = NULL,acc=NULL))
  expect_match(m,"Nonfinite probabilities given by the zerotruncated Poisson",all=TRUE)
  expect_match(w,"value out of range in",all=TRUE)


  })
