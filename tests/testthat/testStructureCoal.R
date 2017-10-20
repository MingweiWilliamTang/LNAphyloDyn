test_that("Test SEIR model Matches the package",{
  library(phydynR)
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
  L1 = chols(M1)
  expect_true(sum(abs(inv2(M1) - solve(M1))) < 0.01)
  expect_equal(L1[1,2],0)
  expect_true(sum(abs(L1 %*% t(L1) - M1)) < 0.0001)
})


test_that("Test SIR model", {
  a = ODE_rk45(c(1000000,5),seq(0,1.4,length.out = 2001),param = c(1.5,20,0.3),x_r = c(1000000,0.6),x_i=c(1,2,0,1),model = "SIR")
  b = ODE_rk45(c(1000000,5),seq(0,1.4,length.out = 2001),param = c(1.5,20,0.3,1000),x_r = c(1000000,0.6),x_i=c(1,2,0,1),model = "SIR")
  expect_equal( sum(abs(a - b)),0)
  KF1 = KF_param(a, c(1.5,33,1000,0.5),gridsize = 50, x_r = c(1000000,0.6),x_i=c(1,3,0,1),model = "SIR")
  KF2 = KF_param_chol(a, c(1.5,33,1000,0.5),gridsize = 50, x_r = c(1000000,0.6),x_i=c(1,3,0,1),model = "SIR")
  expect_equal( sum(abs(KF1$A - KF2$A)),0)
  expect_equal(sum(abs(KF1$Sigma[,,1] - KF2$Sigma[,,1] %*% t(KF2$Sigma[,,1]))) < 0.000001,TRUE)
  aaa = Traj_sim_general_noncentral(a[seq(1,2001,50),],KF2,1.5 )
  expect_equal(aaa$logMultiNorm + aaa$logOrigin, log_like_traj_general2(aaa$SimuTraj, a[seq(1,2001,50),], KF1,50,1.5))
})

