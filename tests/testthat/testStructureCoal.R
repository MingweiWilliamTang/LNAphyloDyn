test_that("Test SEIR model Matches the package",{
  births <- matrix(c('0','parms$beta * I','0','0'),ncol=2)
  rownames(births)=c('E','I')
  migrates <- matrix(c('0','0','parms$mu * E','0'),ncol=2)
  rownames(migrates) = c('E','I')
  deaths <- c(E = '0', I = 'parms$gamma * I')

  dm.det <- build.demographic.process(births=births
                                      ,migrations = migrates,
                                      , deaths = deaths
                                      , parameterNames=c('beta','mu', 'gamma')
                                      , rcpp=FALSE
                                      , sde = FALSE)
  tre.tygf = dm.det(theta = list( beta = 10,mu=3 ,gamma = 6 ),x0  = c(E=5, I = 5 ),t0 = 0,t1 = 10,res = 101)

  path0 =cbind(seq(0,10,length.out = 101),matrix(unlist(rev(tre.tygf[[4]])), ncol=2,byrow = T))

  myODE = ODE_rk45(c(5,5), seq(0,10,length.out = 101),param = c(10/6,3,6),x_r=c(1),x_i = c(0,3,0,2),model = "SEIR2")
  sum(abs(path0[1:10,] - myODE[1:10,]))<5

  expect_equal(sum(abs(path0[1:10,] - myODE[1:10,]))<5,TRUE)
})



test_that("Test Chol decomposition 2by2", {
  M1 = matrix(rnorm(4,0,4),ncol=2)
  M1 = M1 %*% t(M1)
  L1 = chol(M1)
  expect_true(sum(abs(inv2(M1) - solve(M1)) < 0.00000001))
  expect_equal(L1[1,2],0)
  expect_true(sum(abs(L1 %*% t(L1) - M1)) < 0.00000001)
})



