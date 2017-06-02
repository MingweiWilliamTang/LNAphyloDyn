# SIRS LNA_Period

updateAlphas_SIRS = function(MCMC_obj,MCMC_setting,i){
  alpha1 = MCMC_obj$par[1] / (MCMC_setting$N -  MCMC_obj$par[1] -  MCMC_obj$par[2])
  alpha2 = MCMC_obj$par[2] / (MCMC_setting$N -  MCMC_obj$par[1] -  MCMC_obj$par[2])
  #alpha1 = MCMC_obj$par[1] / MCMC_obj$par[2]
  alpha1_new = alpha1 * exp(runif(1,-MCMC_setting$pa, MCMC_setting$pa))
  alpha2_new = alpha2 * exp(runif(1,-MCMC_setting$pa, MCMC_setting$pa))

  state_new = c(X = MCMC_setting$N * alpha1_new / (alpha1_new + alpha2_new + 1),
   Y = MCMC_setting$N * alpha2_new / (alpha1_new + alpha2_new + 1),Z = MCMC_setting$N / (alpha1_new + alpha2_new + 1))

  #state_new = c(X = MCMC_setting$N * alpha1_new/(alpha1_new + 1),
   #             Y = MCMC_setting$N/(alpha1_new + 1))

  Ode_Traj_thin_new <- ODE2(state_new, MCMC_setting$times,MCMC_obj$par[4:7],MCMC_setting$period)


  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = SIRS_KOM_Filter(Ode_Traj_thin_new,MCMC_obj$par[4:7],MCMC_setting$gridsize,period = MCMC_setting$period)


  LatentTraj_new =cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:4] - MCMC_obj$Ode_Traj_coarse[,2:4] +
                          Ode_Traj_coarse_new[,2:4])
  logMultiNorm_new = log_like_trajSIRS(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),betaDyn(MCMC_obj$par[4],MCMC_obj$par[7],Ode_Traj_coarse_new[,1],period = MCMC_setting$period) * MCMC_setting$N,
                                 MCMC_setting$t_correct,
                                 MCMC_setting$gridsize)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[8],MCMC_setting$gridsize)
  }
  if (is.nan(logMultiNorm_new)) {
    logMultiNorm_new = -Inf
    # countInf = countInf + 1
  }
  a = dlnorm(alpha1_new,MCMC_setting$b1,MCMC_setting$a1,log = T) + dlnorm(alpha2_new,MCMC_setting$b2,MCMC_setting$a2,log = T) + coalLog_new + #dgamma(log(alpha2_new),MCMC_setting$b2,MCMC_setting$a2,log = T) +
    logMultiNorm_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog -
    ( MCMC_obj$LogAlpha1)

  #print(logMultiNorm_new-log_like_traj2(MCMC_obj$LatentTraj,MCMC_setting$times,log(state_new),MCMC_obj$par[3],MCMC_obj$par[4],MCMC_setting$gridsize,MCMC_setting$t_correct ))

  #print(logMultiNorm_new - MCMC_obj$logMultiNorm)
  if(is.na(a)){a = - Inf}
  # print(c(logMultiNorm_new,MCMC_obj$logMultiNorm,dgamma(log(alpha1_new),60,MCMC_setting$a1,log = T), dgamma(log(alpha2_new),60,MCMC_setting$a2,log = T)))
  AR = 0
  if (log(runif(1,0,1)) < a) {
    AR = 1
    state = state_new
    MCMC_obj$par[1] = state_new[1]
    MCMC_obj$par[2] = state_new[2]
    MCMC_obj$par[3] = state_new[3]
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$FT = FT_new
    MCMC_obj$LogAlpha1 = dlnorm(alpha1_new,MCMC_setting$b1,MCMC_setting$a1,log = T) +  dlnorm(alpha2_new,MCMC_setting$b2,MCMC_setting$a2,log = T)
    #    MCMC_obj$LogAlpha2 = dgamma(log(alpha2_new),MCMC_setting$b2,MCMC_setting$a2,log = T)
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}




updateS1_SIRS = function(MCMC_obj, MCMC_setting, i){
  s1 = MCMC_obj$par[4] / MCMC_obj$par[5] * MCMC_setting$N
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  s1_new = s1 + runif(1,-MCMC_setting$ps1, MCMC_setting$ps1)
  if(s1_new <1 || s1_new > 10){
    # theta1_new = s1_new * MCMC_obj$par[4] / MCMC_setting$N
    # MCMC_obj$par[3] = theta1_new
    return(list(MCMC_obj = MCMC_obj, AR = 0))
  }
  theta1_new = s1_new * MCMC_obj$par[5] / MCMC_setting$N
  param_new = c(theta1_new, MCMC_obj$par[5:7])
  #    Ode_Traj_thin_new <- ode(y = log(state), times = times,
  #                        func = SIR.log.ode2, parms = param_new)
  Ode_Traj_thin_new <- ODE2(MCMC_obj$par[1:3], MCMC_setting$times,param_new,period = MCMC_setting$period)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = SIRS_KOM_Filter(Ode_Traj_thin_new, param_new, MCMC_setting$gridsize,
                           period = MCMC_setting$period)

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:4] - MCMC_obj$Ode_Traj_coarse[,2:4] +
                             Ode_Traj_coarse_new[,2:4])

  logMultiNorm_new = log_like_trajSIRS(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)
  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),betaDyn(theta1_new,MCMC_obj$par[7],Ode_Traj_coarse_new[,1],period = MCMC_setting$period) * MCMC_setting$N,
                                 MCMC_setting$t_correct,
                                 MCMC_setting$gridsize)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[8],MCMC_setting$gridsize)
  }

  # coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,MCMC_obj$par[8],MCMC_setting$gridsize)

  if(is.nan(logMultiNorm_new)){
    a = -1
    #print("NA")
  }else{
    a = min(c(exp((logMultiNorm_new - MCMC_obj$logMultiNorm + coalLog_new - MCMC_obj$coalLog)
    ),1))
  }
  AR = 0
  if(is.na(a)){
    AR = 0
    print("NA appears when update s1")
    print(s1_new)
    print(MCMC_obj$par[4] / MCMC_obj$par[5] * MCMC_setting$N)
  }else if (runif(1,0,1) < a) {
    AR = 1
    MCMC_obj$par[4] = theta1_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    # s1=s1_new
    MCMC_obj$FT = FT_new
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  #  MCMC_obj$LatentTraj = ESlice(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,log(MCMC_obj$par[1:2]),
  #                           MCMC_setting$Init,MCMC_setting$t_correct,MCMC_obj$par[5],reps = MCMC_setting$reps,MCMC_setting$gridsize)
  #q =  MCMC_obj$logMultiNorm
  # MCMC_obj$logMultiNorm = log_like_trajSIRS(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
  #print(MCMC_obj$logMultiNorm - q)
  # MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,MCMC_obj$LatentTraj,MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)

  return(list(MCMC_obj = MCMC_obj, AR = AR))

}
###########################

updateS2_SIRS = function(MCMC_obj, MCMC_setting, i){
  s2_new = MCMC_obj$par[5] * exp(runif(1,-MCMC_setting$ps2,MCMC_setting$ps2))
  theta1_new = MCMC_obj$par[4] / MCMC_obj$par[5] * s2_new
  theta2_new = s2_new

  param_new = c(theta1 = theta1_new, theta2 = theta2_new, MCMC_obj$par[6:7])

  Ode_Traj_thin_new <- ODE2(MCMC_obj$par[1:3], MCMC_setting$times,
                           param_new,period = MCMC_setting$period)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = SIRS_KOM_Filter(Ode_Traj_thin_new,param_new,MCMC_setting$gridsize,
                           period = MCMC_setting$period)

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:4] - MCMC_obj$Ode_Traj_coarse[,2:4] +
                           Ode_Traj_coarse_new[,2:4])
  logMultiNorm_new = log_like_trajSIRS(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),betaDyn(theta1_new,MCMC_obj$par[7],Ode_Traj_coarse_new[,1],period = MCMC_setting$period) * MCMC_setting$N,
                                 MCMC_setting$t_correct,
                                 MCMC_setting$gridsize)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[8],MCMC_setting$gridsize)
  }

  # coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,MCMC_obj$par[8],MCMC_setting$gridsize)


  if (is.nan(logMultiNorm_new)) {
    a = -1
  }else{
    a = min(c(exp(dlnorm(s2_new,MCMC_setting$c1,MCMC_setting$c2,log = T) + logMultiNorm_new + coalLog_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog -
                    MCMC_obj$LogS2 ),1))
  }
  # print(theta2_new)
  AR = 0
  if(is.na(a)){
    AR = 0
    print("NA appears")
  }else if(runif(1,0,1) < a) {
    AR = 1
    MCMC_obj$par[4] = theta1_new
    MCMC_obj$par[5] = theta2_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$LogS2 = dlnorm(s2_new,MCMC_setting$c1,MCMC_setting$c2,log = T)
    MCMC_obj$FT = FT_new
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new

  }

  return(list(MCMC_obj = MCMC_obj, AR = AR))
}



#####################



updateReSus_SIRS = function(MCMC_obj, MCMC_setting, i){
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  # prior for theta3 log(theta3) ~ N(-29.8,0.5)
  theta3_new = MCMC_obj$par[6] * exp(runif(1,-0.1,0.1))

  param_new = c(MCMC_obj$par[4:5], theta3_new, MCMC_obj$par[7])
  #    Ode_Traj_thin_new <- ode(y = log(state), times = times,
  #                        func = SIR.log.ode2, parms = param_new)
  Ode_Traj_thin_new <- ODE2(MCMC_obj$par[1:3], MCMC_setting$times,param_new,period = MCMC_setting$period)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = SIRS_KOM_Filter(Ode_Traj_thin_new,param_new,MCMC_setting$gridsize,
                           period = MCMC_setting$period)

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:4] - MCMC_obj$Ode_Traj_coarse[,2:4] +
                           Ode_Traj_coarse_new[,2:4])
  logMultiNorm_new = log_like_trajSIRS(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)
  # coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,MCMC_obj$par[8],MCMC_setting$gridsize)

  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),betaDyn(MCMC_obj$par[4],MCMC_obj$par[7],Ode_Traj_coarse_new[,1],period = MCMC_setting$period) * MCMC_setting$N,
                                 MCMC_setting$t_correct,
                                 MCMC_setting$gridsize)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[8],MCMC_setting$gridsize)
  }

  if(is.nan(logMultiNorm_new)){
    a = -1
    #print("NA")
  }else{
    a = min(c(exp((logMultiNorm_new - MCMC_obj$logMultiNorm + coalLog_new - MCMC_obj$coalLog + dlnorm(theta3_new,MCMC_setting$e1,MCMC_setting$e2,T) -
                     dlnorm(MCMC_obj$par[6],MCMC_setting$e1,MCMC_setting$e2,T))
    ),1))
  }
  AR = 0
  if(is.na(a)){
    AR = 0
    print("NA appears when update rsus")
    print(theta3_new)
  }else if (runif(1,0,1) < a) {
    AR = 1
    MCMC_obj$par[6] = theta3_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    # s1=s1_new
    MCMC_obj$FT = FT_new
    # print(MCMC_obj$par)
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }

  return(list(MCMC_obj = MCMC_obj, AR = AR))
}



update_Scale_SIRS = function(MCMC_obj, MCMC_setting, i){
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  # prior for theta4 theta ~ beta(2,2)
  A_new = MCMC_obj$par[7] + (runif(1,-0.18,0.18))
  if(A_new <=0 || A_new >= 1){
    # theta1_new = s1_new * MCMC_obj$par[4] / MCMC_setting$N
    # MCMC_obj$par[3] = theta1_new
    return(list(MCMC_obj = MCMC_obj, AR = 0))
  }
  param_new = c(MCMC_obj$par[4:6],A_new)
  #    Ode_Traj_thin_new <- ode(y = log(state), times = times,
  #                        func = SIR.log.ode2, parms = param_new)
  Ode_Traj_thin_new <- ODE2(MCMC_obj$par[1:3], MCMC_setting$times,param_new,period = MCMC_setting$period)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = SIRS_KOM_Filter(Ode_Traj_thin_new,param_new,MCMC_setting$gridsize,
                           period = MCMC_setting$period)

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:4] - MCMC_obj$Ode_Traj_coarse[,2:4] +
                           Ode_Traj_coarse_new[,2:4])
  logMultiNorm_new = log_like_trajSIRS(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)
  # coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,MCMC_obj$par[8],MCMC_setting$gridsize)

  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),betaDyn(MCMC_obj$par[4],A_new,Ode_Traj_coarse_new[,1],period = MCMC_setting$period) * MCMC_setting$N,
                                 MCMC_setting$t_correct,
                                 MCMC_setting$gridsize)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[8],MCMC_setting$gridsize)
  }

  if(is.nan(logMultiNorm_new)){
    a = -1
    #print("NA")
  }else{
    a = min(c(exp((logMultiNorm_new - MCMC_obj$logMultiNorm + coalLog_new - MCMC_obj$coalLog +  dbeta(A_new,2,2,log=T) - dbeta(MCMC_obj$par[7],2,2,log = T))
    ),1))
  }
  AR = 0
  if(is.na(a)){
    AR = 0
    print("NA appears when update s1")
    print(s1)
  }else if (runif(1,0,1) < a) {
    AR = 1
    MCMC_obj$par[7] = A_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    # s1=s1_new
    MCMC_obj$FT = FT_new
    # print(MCMC_obj$par)
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  #  MCMC_obj$LatentTraj = ESlice(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,log(MCMC_obj$par[1:2]),
  #                           MCMC_setting$Init,MCMC_setting$t_correct,MCMC_obj$par[5],reps = MCMC_setting$reps,MCMC_setting$gridsize)
  #q =  MCMC_obj$logMultiNorm
  # MCMC_obj$logMultiNorm = log_like_trajSIRS(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
  #print(MCMC_obj$logMultiNorm - q)
  # MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,MCMC_obj$LatentTraj,MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)

  return(list(MCMC_obj = MCMC_obj, AR = AR))

}




##################

updateLambda_SIRS = function(MCMC_obj,MCMC_setting, i){
  lambda_new = MCMC_obj$par[8] * exp(runif(1,-0.3,0.3))
  coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(MCMC_obj$LatentTraj),MCMC_setting$t_correct,lambda_new,MCMC_setting$gridsize)
  a = min(c(exp(coalLog_new - MCMC_obj$coalLog + dnorm(lambda_new,MCMC_setting$d1,MCMC_setting$d2,log = T) -
                  MCMC_obj$LogLambda), 1))
  AR = 0

  #print(coalLog_new - MCMC_obj$coalLog + dgamma(log(lambda_new),MCMC_setting$d1,MCMC_setting$d2,log = T) -
  #       MCMC_obj$LogLambda)

  if(runif(1,0,1) < a){
    # rec[i,5] = 1
    AR = 1
    MCMC_obj$LogLambda = dnorm(lambda_new,MCMC_setting$d1,MCMC_setting$d2,log = T)
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$par[8] = lambda_new
  }

  return(list(MCMC_obj = MCMC_obj, AR = AR))
}


###########################


updateTraj = function(MCMC_obj,MCMC_setting,i){

  #print(c(MCMC_obj$par,MCMC_obj$coalLog + MCMC_obj$logMultiNorm))

  MCMC_obj$LatentTraj = ESlice_SIRS(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,MCMC_obj$par[1:3],
                               MCMC_setting$Init,params4 = c(MCMC_obj$par[c(4,7)],MCMC_setting$period, MCMC_setting$N),
                               t_correct = MCMC_setting$t_correct,lambda = MCMC_obj$par[8],reps = MCMC_setting$reps,MCMC_setting$gridsize,
                               volz = MCMC_setting$likelihood == "volz")

  #q =  MCMC_obj$logMultiNorm
  MCMC_obj$logMultiNorm = log_like_trajSIRS(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
  #print(MCMC_obj$logMultiNorm - q)

  if(MCMC_setting$likelihood == "volz"){
    MCMC_obj$coalLog = volz_loglik_nh(MCMC_setting$Init, LogTraj(MCMC_obj$LatentTraj),
                                      betaDyn(MCMC_obj$par[4],MCMC_obj$par[7],MCMC_obj$LatentTraj[,1],MCMC_setting$period) * MCMC_setting$N,
                                      MCMC_setting$t_correct,
                                      MCMC_setting$gridsize)
  }else{
    MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,LogTraj(MCMC_obj$LatentTraj ),MCMC_setting$t_correct,
                                   MCMC_obj$par[8],MCMC_setting$gridsize)
  }


 # MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,LogTraj(MCMC_obj$LatentTraj),MCMC_setting$t_correct,MCMC_obj$par[8],MCMC_setting$gridsize)
  # print(MCMC_obj$coalLog)
  return(list(MCMC_obj=MCMC_obj))
}


MCMC_initialize_SIRS = function(MCMC_setting){ #, prior_par = c(10,20,-2.3,200,40)){

  logMultiNorm = NaN
  coalLog = NaN
  ########
  while(is.nan(logMultiNorm)||is.nan(coalLog)){
    # print(MCMC_setting$control)
    if(is.null(MCMC_setting$control$alpha)){
      alpha1 = exp(rnorm(1,MCMC_setting$b1,MCMC_setting$a1))
      alpha2 = exp(rnorm(1,MCMC_setting$b2,MCMC_setting$a2))
      }else{
      alpha1 = MCMC_setting$control$alpha[1]
      alpha2 = MCMC_setting$control$alpha[2]
    }

    S = MCMC_setting$N * alpha1 / (alpha1 + alpha2 + 1)
    I = MCMC_setting$N * alpha2 / (alpha1 + alpha2 + 1)
    R = MCMC_setting$N / (alpha1 + alpha2 + 1)
    state = c(X = S, Y = I, R = R)

    if(is.null(MCMC_setting$control$s1)){
      s1 = runif(1,1,7)
    }else{
      s1 = MCMC_setting$control$s1
    }

    if(is.null(MCMC_setting$control$s2)){
      s2 = exp(rnorm(1,MCMC_setting$c1, MCMC_setting$c2))
    }else{
      s2 = MCMC_setting$control$s2
    }

    theta2 = s2
    theta1 = s1 * s2 / MCMC_setting$N


    if(is.null(MCMC_setting$control$mu)){
      theta3 = exp(rnorm(1,MCMC_setting$e1,MCMC_setting$e2))
    }else{
      theta3 = MCMC_setting$control$mu
    }

    if(is.null(MCMC_setting$control$A)){
      theta4 = runif(1,0,1)
    }else{
      theta4 = MCMC_setting$control$A
    }

    param = c(theta1 = theta1, theta2 = theta2, theta3 = theta3, theta3 = theta4)
    # print(param)
    #print(state)
    Ode_Traj_thin = ODE2(state, MCMC_setting$times, param,
                           period = MCMC_setting$period)

    Ode_Traj_coarse = Ode_Traj_thin[MCMC_setting$gridset,]


    FT = SIRS_KOM_Filter(Ode_Traj_thin,param,MCMC_setting$gridsize,period = MCMC_setting$period)



    if(is.null(MCMC_setting$control$traj)){
      Latent = Traj_sim_SIRS(state,Ode_Traj_coarse,FT,MCMC_setting$t_correct)
      LatentTraj = Latent$SimuTraj
      logMultiNorm = Latent$loglike
    }else{
      LatentTraj = MCMC_setting$control$traj
      if( sum(abs(LatentTraj[1,c(2,3)]) - c(S,I)) > 1){
        print("not consistent")
      }
      logMultiNorm = log_like_trajSIRS(LatentTraj,Ode_Traj_coarse,
                                         FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
    }

    if(is.null(MCMC_setting$control$lambda)){
      lambda = rnorm(1,MCMC_setting$d1,MCMC_setting$d2)
    }else{
      lambda = MCMC_setting$control$lambda
    }

    if(MCMC_setting$likelihood == "volz"){
      coalLog = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj),betaDyn(theta1,theta4,LatentTraj[,1],MCMC_setting$period) * MCMC_setting$N,
                               MCMC_setting$t_correct,
                               MCMC_setting$gridsize)

      #print(MCMC_setting$t_correct)
    }else{
      coalLog= coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj ),MCMC_setting$t_correct,
                           lambda,MCMC_setting$gridsize)
    }

    print(coalLog)

    if(!is.nan((coalLog))){
      coalLog = ifelse(coalLog<= - 10000000,NaN, coalLog)
    }
    plot(LatentTraj[,1],LatentTraj[,3],type="l")
  }
  LogAlpha1 = dlnorm(alpha1,MCMC_setting$b1,MCMC_setting$a1,log = T)
  LogMu = dlnorm(theta3, MCMC_setting$e1, MCMC_setting$e2, log = T)
  #LogAlpha2 = dgamma(log(alpha2),MCMC_setting$b1,MCMC_setting$a2,log = T)
  LogS2 = dlnorm(s2,MCMC_setting$c1,MCMC_setting$c2,log = T)
  LogLambda = dnorm(lambda,MCMC_setting$d1,MCMC_setting$d2,log = T)

  #print(log_like_trajSIR_BD(LatentTraj,Ode_Traj_coarse,FT,MCMC_setting$gridsize,90))
  #print()
  #plot(Ode_Traj_coarse[,3])
  plot(LatentTraj[,1],LatentTraj[,3],type="l")
  if(MCMC_setting$likelihood == "volz"){
    paras =  c(S,I,R,theta1,theta2,theta3,theta4)
  }else{
    paras =  c(S,I,R,theta1,theta2,theta3,theta4,lambda)
  }
  MCMC_obj = list(par = paras,LatentTraj = LatentTraj, logMultiNorm = logMultiNorm,
                  Ode_Traj_coarse = Ode_Traj_coarse, FT = FT, coalLog = coalLog,
                  LogAlpha1 = LogAlpha1, LogS2 = LogS2, LogMu,LogLambda = LogLambda)
  ##########
  # MCMC_para = matrix(nrow = niter,ncol = 2)
  cat("Initialize MCMC \n")
  print(paste("size = ", MCMC_setting$N))
  print(paste("S0 = ",S," I0 = ", I, "R0 = ", R))
  print(paste("R0 = ",s1," gamma = ", s2, " beta = ", theta1))
  print(paste("mu = ", theta3, " A = ", theta4))
  return(MCMC_obj)
}




SIRS_LNA_MCMC_standard = function(coal_obs,times,t_correct,N,gridsize=1000, niter = 1000, burn = 0, thin = 5,period = 40,
                                           a1 = 10, a2 = 20, b1 = 60 , b2 = 60, c1=-2.3,c2 = 0.4,d1 = 250, d2 =40, e1 = -2.8, e2 = 0.5,
                                           pa = 0.1, ps1 = 0.25, ps2 = 0.5, pga = 0, pA = 0.18,
                                  control = list(), updateVec = c(1,1,1,1,1), likelihood = "volz"){

  MCMC_setting = MCMC_setup(coal_obs,times,t_correct,N,gridsize,niter,burn,thin,period = period,
                            a1, a2,b1,b2,c1,c2,d1, d2,e1,e2,
                            pa,ps1,ps2,pga,pA,control = control,likelihood = likelihood)

  MCMC_obj = MCMC_initialize_SIRS(MCMC_setting)
  d = 0
  if(MCMC_setting$likelihood == "volz"){
    params = matrix(nrow = niter, ncol = 7)
    ARMS = matrix(nrow = niter, ncol = 7)
    d = 7
  }else{
    params = matrix(nrow = niter, ncol = 8)
    ARMS = matrix(nrow = niter, ncol = 8)
    d = 8
  }

  #print(c(MCMC_obj$coalLog,MCMC_obj$logMultiNorm))
  # MCMC_obj$LogS2 =
  # MCMC_obj$LogLambda

  #MCMC_obj$par[1] = 1500
  #MCMC_obj$par[2] = 1000
  #MCMC_obj$par[3] = 9500
  #MCMC_obj$par[5] = 0.1
  #MCMC_obj$par[6] = 0.005
  #MCMC_obj$par[7] = 0.5
 # MCMC_obj$par[8] = 400
  l = numeric(niter)
  l1 = l
  l2 = l
  l3 = l
  tjs = NULL
  for (i in 1:MCMC_setting$niter) {
    if (i %% 100 == 0) {
      print(i)
      plot(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,3],type="l")
    }
    ARvec = numeric(d)

    #  print(c(MCMC_obj$coalLog,MCMC_obj$logMultiNorm))
    #MCMC_obj = step1$MCMC_obj
    if(updateVec[1] == 1){
      step1 = updateAlphas_SIRS(MCMC_obj,MCMC_setting,i)
      ARvec[1] = step1$AR
      MCMC_obj = step1$MCMC_obj
    }
    if(updateVec[2] == 1){
      step2 = updateS1_SIRS(MCMC_obj,MCMC_setting,i)
      ARvec[4] = step2$AR
      MCMC_obj = step2$MCMC_obj
    }
    if(updateVec[3] == 1){
       stepRecover = updateS2_SIRS(MCMC_obj,MCMC_setting,i)
      ARvec[5] = stepRecover$AR
       MCMC_obj = stepRecover$MCMC_obj
    }
    if(updateVec[4] == 1){
      stepReS = updateReSus_SIRS(MCMC_obj,MCMC_setting,i)
      ARvec[6] = stepReS$AR
      MCMC_obj = stepReS$MCMC_obj
    }
    if(updateVec[5] == 1){
      stepA = update_Scale_SIRS(MCMC_obj, MCMC_setting,i)
      ARvec[7] = stepA$AR
      MCMC_obj = stepA$MCMC_obj
    }
    if(updateVec[6] == 1){
    MCMC_obj = updateTraj(MCMC_obj,MCMC_setting,i)$MCMC_obj
    }
    if(updateVec[7] == 1){
      MCMC_obj = updateLambda_SIRS(MCMC_obj,MCMC_setting,i)$MCMC_obj
    }
    # step4 = updateLambda_SIRS(MCMC_obj,MCMC_setting,i)
  # ARvec[8] = step4$AR
   #MCMC_obj = step4$MCMC_obj
   ARMS[i,] = ARvec
    tjs = abind(tjs,MCMC_obj$LatentTraj,along = 3)
    params[i,] = MCMC_obj$par
    l[i] =  MCMC_obj$logMultiNorm #+ MCMC_obj$LogAlpha1 + MCMC_obj$LogAlpha2
    l1[i] = MCMC_obj$coalLog
    l2[i] = MCMC_obj$logMultiNorm + MCMC_obj$coalLog
    l3[i] = MCMC_obj$LogAlpha1 + MCMC_obj$logMultiNorm + MCMC_obj$LogS2 + MCMC_obj$coalLog + MCMC_obj$LogLambda
  }
  return(list(par = params,Trajectory = tjs,l=l,l1=l1,l2 = l2, l3 =l3,AR = ARMS))
}


