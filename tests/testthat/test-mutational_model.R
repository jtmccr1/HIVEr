context("Testing mutational model functions")

test_that("Probabilities make sense",{
  g_ft.V<-Vectorize(g_ft,vectorize.args = "x")

  test_prob <- function(mu,Ne,t,gc_ul){
    found<-integrate(g_ft.V,lower=0.02,upper=1,
             mu=mu,Ne=Ne,t=t,gc_ul,accuracy_stringent)
    lost<-not_detected(mu,Ne,t,gc_ul,accuracy_stringent)
    lost+found$value
  }
  trial_1<-test_prob(1e-8,10,30,1e4)
  trial_2<-test_prob(1e-8,70,30,1e3)
  trial_3<-test_prob(1e-2,30,10,1e6)
expect_true(abs(1-trial_1)<0.0002)
expect_true(abs(1-trial_2)<0.0002)
expect_true(abs(1-trial_3)<0.0002)

}
)


