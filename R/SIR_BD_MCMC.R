# SIRLNA_Period
updateAlphas_BD = function(MCMC_obj,MCMC_setting,i){
  #alpha1 = MCMC_obj$par[1] / (MCMC_setting$N -  MCMC_obj$par[1] -  MCMC_obj$par[2])
  #alpha2 = MCMC_obj$par[2] / (MCMC_setting$N -  MCMC_obj$par[1] -  MCMC_obj$par[2])
  alpha1 = MCMC_obj$par[1] / MCMC_obj$par[2]
  alpha1_new = alpha1 * exp(runif(1,- MCMC_setting$pa, MCMC_setting$pa))

  #alpha2_new = alpha2 * exp(runif(1,-0.2,0.2))
  #state_new = c(X = MCMC_setting$N * alpha1_new / (alpha1_new + alpha2_new + 1),
  #             Y = MCMC_setting$N * alpha2_new / (alpha1_new + alpha2_new + 1))
  state_new = c(X = MCMC_setting$N * alpha1_new/(alpha1_new + 1),
                Y = MCMC_setting$N/(alpha1_new + 1))

  Ode_Traj_thin_new <- SIR_BD_ODE(state_new, MCMC_setting$times,MCMC_obj$par[3:6],MCMC_setting$N,MCMC_setting$period)


  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = SIR_BD_KOM_Filter(Ode_Traj_thin_new,MCMC_obj$par[3:6],MCMC_setting$gridsize,MCMC_setting$N,period = MCMC_setting$period)


    LatentTraj_new =cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                            Ode_Traj_coarse_new[,2:3])
    logMultiNorm_new = log_like_trajSIR_BD(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

    if(MCMC_setting$likelihood == "volz"){
      coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),betaDyn(MCMC_obj$par[3],MCMC_obj$par[6],Ode_Traj_coarse_new[,1],period = MCMC_setting$period) * MCMC_setting$N,
                                   MCMC_setting$t_correct,
                                   MCMC_setting$gridsize)
    }else{
      coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                                MCMC_obj$par[7],MCMC_setting$gridsize)
    }
    if (is.nan(logMultiNorm_new)) {
    logMultiNorm_new = -Inf
    # countInf = countInf + 1
  }
  a = dnorm(log(alpha1_new),MCMC_setting$b1,MCMC_setting$a1,log = T) + coalLog_new + #dgamma(log(alpha2_new),MCMC_setting$b2,MCMC_setting$a2,log = T) +
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
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$FT = FT_new
    MCMC_obj$LogAlpha1 = dnorm(log(alpha1_new),MCMC_setting$b1,MCMC_setting$a1,log = T)
    #    MCMC_obj$LogAlpha2 = dgamma(log(alpha2_new),MCMC_setting$b2,MCMC_setting$a2,log = T)
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}



updateS1_SIR_BD = function(MCMC_obj, MCMC_setting, i){
  s1 = MCMC_obj$par[3] / MCMC_obj$par[4] * MCMC_setting$N
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  s1_new = s1 + runif(1,-MCMC_setting$ps1, MCMC_setting$ps1)
  if(s1_new <1 || s1_new > 100){
    # theta1_new = s1_new * MCMC_obj$par[4] / MCMC_setting$N
    # MCMC_obj$par[3] = theta1_new
    return(list(MCMC_obj = MCMC_obj, AR = 0))
  }
  theta1_new = s1_new * MCMC_obj$par[4] / MCMC_setting$N
  param_new = c(theta1_new, MCMC_obj$par[4:6])
  #    Ode_Traj_thin_new <- ode(y = log(state), times = times,
  #                        func = SIR.log.SIR_BD_ODE, parms = param_new)
  Ode_Traj_thin_new <- SIR_BD_ODE(MCMC_obj$par[1:2], MCMC_setting$times,param_new,MCMC_setting$N,MCMC_setting$period)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = SIR_BD_KOM_Filter(Ode_Traj_thin_new,param_new,MCMC_setting$gridsize,MCMC_setting$N,period = MCMC_setting$period)

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                           Ode_Traj_coarse_new[,2:3])
  logMultiNorm_new = log_like_trajSIR_BD(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),betaDyn(theta1_new,MCMC_obj$par[6],Ode_Traj_coarse_new[,1],MCMC_setting$period) * MCMC_setting$N,
                                 MCMC_setting$t_correct,
                              MCMC_setting$gridsize)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[7],MCMC_setting$gridsize)
  }

  if(is.nan(logMultiNorm_new)){
    a = -Inf
    #print("NA")
  }else{
    a = logMultiNorm_new - MCMC_obj$logMultiNorm +
                     coalLog_new - MCMC_obj$coalLog
  }
  AR = 0
  if(is.na(a)){
    AR = 0
    print("NA appears when update s1")
    print(s1_new)
    print(MCMC_obj$par[3] / MCMC_obj$par[4] * MCMC_setting$N)
  }else if (log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[3] = theta1_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    # s1=s1_new
    MCMC_obj$FT = FT_new

    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }

  return(list(MCMC_obj = MCMC_obj, AR = AR))

}
###########################

updateS2_SIR_BD = function(MCMC_obj, MCMC_setting, i){

  s2_new = MCMC_obj$par[4] * exp(runif(1,-MCMC_setting$ps2,MCMC_setting$ps2))
  theta1_new = MCMC_obj$par[3] / MCMC_obj$par[4] * s2_new
  theta2_new = s2_new

  param_new = c(theta1 = theta1_new, theta2 = theta2_new, MCMC_obj$par[5:6])


  Ode_Traj_thin_new <- SIR_BD_ODE(MCMC_obj$par[1:2], MCMC_setting$times,
                            param_new, MCMC_setting$N,MCMC_setting$period)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = SIR_BD_KOM_Filter(Ode_Traj_thin_new,param_new,MCMC_setting$gridsize, MCMC_setting$N,period = MCMC_setting$period)

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                           Ode_Traj_coarse_new[,2:3])
  logMultiNorm_new = log_like_trajSIR_BD(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),
                                 betaDyn(theta1_new,MCMC_obj$par[6],Ode_Traj_coarse_new[,1],MCMC_setting$period) * MCMC_setting$N,
                                 MCMC_setting$t_correct,
                                 MCMC_setting$gridsize)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[7],MCMC_setting$gridsize)

  }

  if (is.nan(logMultiNorm_new)) {
    a = -Inf
  }else{
    a = dnorm(log(s2_new),MCMC_setting$c1,MCMC_setting$c2,log = T) + log(s2_new) + logMultiNorm_new + coalLog_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog -
                    MCMC_obj$LogS2 - log(MCMC_obj$par[4])
  }
  # print(theta2_new)
  AR = 0
  if(is.na(a)){
    AR = 0
    print("NA appears")
  }else if(log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[3] = theta1_new
    MCMC_obj$par[4] = theta2_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$LogS2 = dnorm(log(s2_new),MCMC_setting$c1,MCMC_setting$c2,log = T)
    MCMC_obj$FT = FT_new
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new

  }

  return(list(MCMC_obj = MCMC_obj, AR = AR))
}

##############
updateS1S2_SIR_BD = function(MCMC_obj, MCMC_setting, i){
  s1 = MCMC_obj$par[3] / MCMC_obj$par[4] * MCMC_setting$N
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  s1_new = s1 + runif(1,-MCMC_setting$ps1, MCMC_setting$ps1)
  if(s1_new <1 || s1_new > 100){
    # theta1_new = s1_new * MCMC_obj$par[4] / MCMC_setting$N
    # MCMC_obj$par[3] = theta1_new
    return(list(MCMC_obj = MCMC_obj, AR = 0))
  }
  s2_new = MCMC_obj$par[4] * exp(runif(1,-MCMC_setting$ps2,MCMC_setting$ps2))


  theta1_new =  s1_new * s2_new / MCMC_setting$N
  theta2_new = s2_new

  param_new = c(theta1 = theta1_new, theta2 = theta2_new, MCMC_obj$par[5:6])


  Ode_Traj_thin_new <- SIR_BD_ODE(MCMC_obj$par[1:2], MCMC_setting$times,
                                  param_new, MCMC_setting$N,MCMC_setting$period)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = SIR_BD_KOM_Filter(Ode_Traj_thin_new,param_new,MCMC_setting$gridsize, MCMC_setting$N,period = MCMC_setting$period)

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                           Ode_Traj_coarse_new[,2:3])
  logMultiNorm_new = log_like_trajSIR_BD(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),
                                 betaDyn(MCMC_obj$par[3],MCMC_obj$par[6],Ode_Traj_coarse_new[,1],MCMC_setting$period) * MCMC_setting$N,
                                 MCMC_setting$t_correct,
                                 MCMC_setting$gridsize)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[7],MCMC_setting$gridsize)

  }

  if (is.nan(logMultiNorm_new)) {
    a = -Inf
  }else{
    a = dlnorm(s2_new,MCMC_setting$c1,MCMC_setting$c2,log = T) + log(s2_new) + logMultiNorm_new + coalLog_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog -
      MCMC_obj$LogS2 - log(MCMC_obj$par[4])
  }
  # print(theta2_new)
  AR = 0
  if(is.na(a)){
    AR = 0
    if(is.na(coalLog_new)){
      print("NA appears in likelihood")
      print(LatentTraj_new)
    }else{
    print("NA appears in S1 S2")
    }
  }else if(log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[3] = theta1_new
    MCMC_obj$par[4] = theta2_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$LogS2 = dnorm(log(s2_new),MCMC_setting$c1,MCMC_setting$c2,log = T)
    MCMC_obj$FT = FT_new
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new

  }

  return(list(MCMC_obj = MCMC_obj, AR = AR))
}

#####################



updateReSus_SIR_BD = function(MCMC_obj, MCMC_setting, i){
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  # prior for theta3 log(theta3) ~ N(-29.8,0.5)
  theta3_new = MCMC_obj$par[5] * exp(runif(1,-MCMC_setting$pga, MCMC_setting$pga))

  param_new = c(MCMC_obj$par[3:4], theta3_new, MCMC_obj$par[7])

  #    Ode_Traj_thin_new <- ode(y = log(state), times = times,
  #                        func = SIR.log.SIR_BD_ODE, parms = param_new)
  Ode_Traj_thin_new <- SIR_BD_ODE(MCMC_obj$par[1:2], MCMC_setting$times,param_new, MCMC_setting$N,
                                  MCMC_setting$period)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = SIR_BD_KOM_Filter(Ode_Traj_thin_new,param_new,MCMC_setting$gridsize, MCMC_setting$N,period = MCMC_setting$period)

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                           Ode_Traj_coarse_new[,2:3])

  logMultiNorm_new = log_like_trajSIR_BD(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,
                                         MCMC_setting$t_correct)

  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),
                                 betaDyn(MCMC_obj$par[3],MCMC_obj$par[6],Ode_Traj_coarse_new[,1],MCMC_setting$period) * MCMC_setting$N,
                                 MCMC_setting$t_correct,
                                 MCMC_setting$gridsize)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[7],MCMC_setting$gridsize)
  }
 # coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),
  #                          MCMC_setting$t_correct,MCMC_obj$par[7],MCMC_setting$gridsize)

  if(is.nan(logMultiNorm_new)){
    a = - Inf
    #print("NA")
  }else{
    a = logMultiNorm_new - MCMC_obj$logMultiNorm + coalLog_new - MCMC_obj$coalLog + dnorm(theta3_new,MCMC_setting$e1,MCMC_setting$e2,T) -
                     MCMC_obj$LogMu
                     #dnorm(MCMC_obj$par[6],MCMC_setting$e1,MCMC_setting$e2,T))
  }
  AR = 0
  if(is.na(a)){
    AR = 0
    print("NA appears when update rsus")
    print(theta3_new)
  }else if (log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[5] = theta3_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    # s1=s1_new
    MCMC_obj$FT = FT_new
    MCMC_obj$LogMu = dnorm(theta3_new,MCMC_setting$e1,MCMC_setting$e2,T)
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  #  MCMC_obj$LatentTraj = ESlice(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,log(MCMC_obj$par[1:2]),
  #                           MCMC_setting$Init,MCMC_setting$t_correct,MCMC_obj$par[5],reps = MCMC_setting$reps,MCMC_setting$gridsize)
  #q =  MCMC_obj$logMultiNorm
  # MCMC_obj$logMultiNorm = log_like_trajSIR_BD(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
  #print(MCMC_obj$logMultiNorm - q)
  # MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,MCMC_obj$LatentTraj,MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)

  return(list(MCMC_obj = MCMC_obj, AR = AR))
}



update_Scale_SIR_BD = function(MCMC_obj, MCMC_setting, i){
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  # prior for theta4 theta ~ beta(2,2)
  A_new = MCMC_obj$par[6] + (runif(1,-MCMC_setting$pA,MCMC_setting$pA))
  if(A_new <=0 || A_new >= 1){
    # theta1_new = s1_new * MCMC_obj$par[4] / MCMC_setting$N
    # MCMC_obj$par[3] = theta1_new
    return(list(MCMC_obj = MCMC_obj, AR = 0))
  }
  param_new = c(MCMC_obj$par[3:5],A_new)
  #    Ode_Traj_thin_new <- ode(y = log(state), times = times,
  #                        func = SIR.log.SIR_BD_ODE, parms = param_new)
  Ode_Traj_thin_new <- SIR_BD_ODE(MCMC_obj$par[1:2], MCMC_setting$times,param_new, MCMC_setting$N,MCMC_setting$period)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = SIR_BD_KOM_Filter(Ode_Traj_thin_new,param_new,MCMC_setting$gridsize, MCMC_setting$N,period = MCMC_setting$period)

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                           Ode_Traj_coarse_new[,2:3])
  logMultiNorm_new = log_like_trajSIR_BD(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),
                                 betaDyn(MCMC_obj$par[3],A_new,Ode_Traj_coarse_new[,1],MCMC_setting$period) * MCMC_setting$N,
                                 MCMC_setting$t_correct,
                                 MCMC_setting$gridsize)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[7],MCMC_setting$gridsize)
  }


  if(is.nan(logMultiNorm_new)){
    a = -1
    #print("NA")
  }else{
    a = logMultiNorm_new - MCMC_obj$logMultiNorm + coalLog_new - MCMC_obj$coalLog +  dbeta(A_new,2,2,log=T) - dbeta(MCMC_obj$par[6],2,2,log = T)
  }
  AR = 0
  if(is.na(a)){
    AR = 0
    print("NA appears when update A")
    print(A_new)
  }else if (log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[6] = A_new
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




##################

updateLambda_SIR_BD = function(MCMC_obj,MCMC_setting, i){
  lambda_new = MCMC_obj$par[7] * exp(runif(1,-0.3,0.3))
  coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(MCMC_obj$LatentTraj),MCMC_setting$t_correct,lambda_new,MCMC_setting$gridsize)
  a = coalLog_new - MCMC_obj$coalLog + dgamma(log(lambda_new),MCMC_setting$d1,MCMC_setting$d2,log = T) -
                  MCMC_obj$LogLambda
  AR = 0

  #print(coalLog_new - MCMC_obj$coalLog + dgamma(log(lambda_new),MCMC_setting$d1,MCMC_setting$d2,log = T) -
  #       MCMC_obj$LogLambda)
  if(is.nan(a)){
    print(coalLog_new)
  }else if(log(runif(1,0,1)) < a){
    # rec[i,5] = 1
    AR = 1
    MCMC_obj$LogLambda = dgamma(log(lambda_new),MCMC_setting$d1,MCMC_setting$d2,log = T)
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$par[7] = lambda_new
  }

  return(list(MCMC_obj = MCMC_obj, AR = AR))
}




updateTraj_BD = function(MCMC_obj,MCMC_setting,i){
  #print(c(MCMC_obj$par,MCMC_obj$coalLog + MCMC_obj$logMultiNorm))
  #print(MCMC_obj$par[3])
  #print(betaDyn(MCMC_obj$par[3],MCMC_obj$par[6],MCMC_obj$LatentTraj[,1]))
  MCMC_obj$LatentTraj = ESlice_SIR_BD(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,
                                      MCMC_obj$FT,MCMC_obj$par[1:2], MCMC_setting$Init,betaN =
                                    betaDyn(MCMC_obj$par[3],MCMC_obj$par[6],MCMC_obj$LatentTraj[,1],MCMC_setting$period) * MCMC_setting$N,
                                    MCMC_setting$t_correct,lambda = MCMC_obj$par[7],
                                    reps = MCMC_setting$reps,MCMC_setting$gridsize,
                                    volz = (MCMC_setting$likelihood == "volz"))
  #q =  MCMC_obj$logMultiNorm
  MCMC_obj$logMultiNorm = log_like_trajSIR_BD(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,
                                              MCMC_obj$FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
  #print(MCMC_obj$logMultiNorm - q)

  if(MCMC_setting$likelihood == "volz"){
    MCMC_obj$coalLog = volz_loglik_nh(MCMC_setting$Init, LogTraj(MCMC_obj$LatentTraj),
                                 betaDyn(MCMC_obj$par[3],MCMC_obj$par[6],MCMC_obj$LatentTraj[,1],MCMC_setting$period) * MCMC_setting$N,
                                 MCMC_setting$t_correct,
                                 MCMC_setting$gridsize)
  }else{
    MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,LogTraj(MCMC_obj$LatentTraj ),MCMC_setting$t_correct,
                              MCMC_obj$par[7],MCMC_setting$gridsize)
  }

  #MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,LogTraj(MCMC_obj$LatentTraj),MCMC_setting$t_correct,MCMC_obj$par[8],MCMC_setting$gridsize)
  # print(MCMC_obj$coalLog)
  return(list(MCMC_obj=MCMC_obj))
}




MCMC_setup = function(coal_obs,times,t_correct,N,gridsize=50,niter = 1000,burn = 500,thin = 5, period = 40,
                      a1 = 10, a2 = 20,b1 = 60, b2= 60, c1=-2.3,c2 = 0.4,d1 = 200, d2 =40, e1 = -2.8, e2 = 0.5,
                      pa = 0.1, ps1 = 0.25, ps2 = 0.5, pga = 0, pA = 0.18, control = list(), likelihood = "volz"){
  gridset = seq(1,length(times),by=gridsize)
  grid = times[gridset]
  Init = coal_lik_init(coal_obs$samp_times, coal_obs$n_sampled, coal_obs$coal_times, grid, t_correct)
  MCMC_setting = list(Init = Init,times = times,t_correct = t_correct,N = N,
                      gridsize=gridsize,gridset = gridset, niter = niter,burn = burn,thin = thin,period = period,
                      a1 = a1, a2 = a2,b1 =b1, b2 = b2, c1= c1, c2 = c2,d1 = d1, d2 = d2,e1 = e1, e2 = e2,
                      pa = pa, ps1 = ps1, ps2 = ps2, pga = pga, pA = pA,control = control,
                      reps=1, likelihood = likelihood)
  cat("MCMC set up ready \n")

  return(MCMC_setting)
}


MCMC_initialize2 = function(MCMC_setting){ #, prior_par = c(10,20,-2.3,200,40)){

  logMultiNorm = NaN
  coalLog = NaN
  ########
  while(is.nan(logMultiNorm)||is.nan(coalLog)){
   # print(MCMC_setting$control)
    if(is.null(MCMC_setting$control$alpha)){
      alpha1 = exp(rnorm(1,MCMC_setting$b1,MCMC_setting$a1))
    }else{
      alpha1 = MCMC_setting$control$alpha
    }

    S = MCMC_setting$N * alpha1 / (alpha1  + 1)
    I = MCMC_setting$N / (alpha1 + 1)
    state = c(X = S, Y = I)

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
    Ode_Traj_thin = SIR_BD_ODE(state, MCMC_setting$times,
                               param,MCMC_setting$N, period = MCMC_setting$period)

    Ode_Traj_coarse = Ode_Traj_thin[MCMC_setting$gridset,]


    FT = SIR_BD_KOM_Filter(Ode_Traj_thin,param,MCMC_setting$gridsize, MCMC_setting$N,period = MCMC_setting$period)



    if(is.null(MCMC_setting$control$traj)){
      Latent = Traj_sim_SIR_BD(state,Ode_Traj_coarse,FT,MCMC_setting$t_correct)
      LatentTraj = Latent$SimuTraj
      logMultiNorm = Latent$loglike
    }else{
      LatentTraj = MCMC_setting$control$traj
      if( sum(abs(LatentTraj[1,c(2,3)]) - c(S,I)) > 1){
        print("not consistent")
      }
      logMultiNorm = log_like_trajSIR_BD(LatentTraj,Ode_Traj_coarse,
                                   FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
    }

    if(is.null(MCMC_setting$control$lambda)){
      lambda = exp(rgamma(1,MCMC_setting$d1,MCMC_setting$d2))
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
  LogAlpha1 = dnorm(log(alpha1),MCMC_setting$b1,MCMC_setting$a1,log = T)
  LogMu = dnorm(log(theta3), MCMC_setting$e1, MCMC_setting$e2, log = T)
  #LogAlpha2 = dgamma(log(alpha2),MCMC_setting$b1,MCMC_setting$a2,log = T)
  LogS2 = dnorm(log(s2),MCMC_setting$c1,MCMC_setting$c2,log = T)
  LogLambda = dgamma(log(lambda),MCMC_setting$d1,MCMC_setting$d2,log = T)

  #print(log_like_trajSIR_BD(LatentTraj,Ode_Traj_coarse,FT,MCMC_setting$gridsize,90))
  #print()
  #plot(Ode_Traj_coarse[,3])
   plot(LatentTraj[,1],LatentTraj[,3],type="l")
  if(MCMC_setting$likelihood == "volz"){
    paras =  c(S,I,theta1,theta2,theta3,theta4)
  }else{
    paras =  c(S,I,theta1,theta2,theta3,theta4,lambda)
  }
  MCMC_obj = list(par = paras,LatentTraj = LatentTraj, logMultiNorm = logMultiNorm,
                  Ode_Traj_coarse = Ode_Traj_coarse, FT = FT, coalLog = coalLog,
                  LogAlpha1 = LogAlpha1, LogS2 = LogS2, LogMu,LogLambda = LogLambda)
  ##########
  # MCMC_para = matrix(nrow = niter,ncol = 2)
  cat("Initialize MCMC \n")
  print(paste("size = ", MCMC_setting$N))
  print(paste("S0 = ",S," I0 = ", I))
  print(paste("R0 = ",s1," gamma = ", s2, " beta = ", theta1))
  print(paste("mu = ", theta3, " A = ", theta4))
  return(MCMC_obj)
}



SIR_BD_LNA_MCMC = function(coal_obs,times,t_correct,N,gridsize=1000, niter = 1000, burn = 0, thin = 5,period = 40,
                                  a1 = 10, a2 = 20, b1 = 60 , b2 = 60, c1=-2.3,c2 = 0.4,d1 = 250, d2 =40, e1 = -2.8, e2 = 0.5,
                           pa = 0.1, ps1 = 0.25, ps2 = 0.5, pga = 0, pA = 0.18, control = list(), updateVec = c(1,1,1,1,1), likelihood = "volz"){
  MCMC_setting = MCMC_setup(coal_obs,times,t_correct,N,gridsize,niter,burn,thin,period = period,
                            a1, a2,b1,b2,c1,c2,d1, d2,e1,e2,
                            pa,ps1,ps2,pga,pA,control = control,likelihood = likelihood)
  MCMC_obj = MCMC_initialize2(MCMC_setting)
  if(MCMC_setting$likelihood == "volz"){
    params = matrix(nrow = niter, ncol = 6)
    ARMS = matrix(nrow = niter, ncol = 6)
  }else{
    params = matrix(nrow = niter, ncol = 7)
    ARMS = matrix(nrow = niter, ncol = 7)
  }
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
    ARvec = numeric(dim(ARMS)[2])
    if(updateVec[1] == 1){
      step1 = updateAlphas_BD(MCMC_obj,MCMC_setting,i)
     # print(c(MCMC_obj$coalLog,MCMC_obj$logMultiNorm))
      MCMC_obj = step1$MCMC_obj
      ARvec[1] = step1$AR
    }
    if(updateVec[2] == 1 && updateVec[3] == 0){
      step2 = updateS1_SIR_BD(MCMC_obj,MCMC_setting,i)
      ARvec[3] = step2$AR
      MCMC_obj = step2$MCMC_obj
    }else if(updateVec[3] == 1 && updateVec[2] == 0){
      stepRecover = updateS2_SIR_BD(MCMC_obj,MCMC_setting,i)
      ARvec[4] = stepRecover$AR
      MCMC_obj = stepRecover$MCMC_obj
    }else if(updateVec[3] == 1 && updateVec[2] == 1){
      stepJoint = updateS1S2_SIR_BD(MCMC_obj,MCMC_setting,i)
      ARvec[4] = stepJoint$AR
      MCMC_obj = stepJoint$MCMC_obj
    }

    if(updateVec[4] == 1){
      stepReS = updateReSus_SIR(MCMC_obj,MCMC_setting,i)
      ARvec[5] = stepReS$AR
      MCMC_obj = stepReS$MCMC_obj
    }
    if(updateVec[5] == 1){
      stepA = update_Scale_SIR_BD(MCMC_obj, MCMC_setting,i)
      ARvec[6] = stepA$AR
      MCMC_obj = stepA$MCMC_obj
    }
    if(updateVec[6] == 1){
      MCMC_obj = updateTraj_BD(MCMC_obj,MCMC_setting,i)$MCMC_obj
    }
    if(length(updateVec)>6 && updateVec[7] == 1){
      steplambda = updateLambda_SIR_BD(MCMC_obj,MCMC_setting,i)
      ARvec[7] = steplambda$AR
      MCMC_obj = steplambda$MCMC_obj
    }
    #step4 = updateLambda_SIRS(MCMC_obj,MCMC_setting,i)
    #ARvec[8] = step4$AR
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


effpopfun = function(Traj,beta=0,lambda=1, volz = FALSE){
  if(volz){
    return(1 /(2 * Traj[,3] * beta / Traj[,2]))
  }else{
    return(Traj[,3] / lambda)
  }
}


