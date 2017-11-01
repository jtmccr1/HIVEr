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
