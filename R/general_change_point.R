  # SIRLNA_Period

updateAlphas_General = function(MCMC_obj,MCMC_setting,i){
  #alpha1 = MCMC_obj$par[1] / (MCMC_setting$N -  MCMC_obj$par[1] -  MCMC_obj$par[2])
  #alpha2 = MCMC_obj$par[2] / (MCMC_setting$N -  MCMC_obj$par[1] -  MCMC_obj$par[2])
  alpha1 = MCMC_obj$par[1] / MCMC_obj$par[2]
  alpha1_new = alpha1 * exp(runif(1,- MCMC_setting$pa, MCMC_setting$pa))

  #alpha2_new = alpha2 * exp(runif(x1,-0.2,0.2))
  #state_new = c(X = MCMC_setting$N * alpha1_new / (alpha1_new + alpha2_new + 1),
  #             Y = MCMC_setting$N * alpha2_new / (alpha1_new + alpha2_new + 1))
  state_new = c(X = MCMC_setting$x_r[1] * alpha1_new/(alpha1_new + 1),
                Y = MCMC_setting$x_r[1] / (alpha1_new + 1))

  Ode_Traj_thin_new <- General_ODE_rk45(state_new,MCMC_setting$times,MCMC_obj$par[3:(sum(MCMC_setting$x_i)+2)],
                                        MCMC_setting$x_r,MCMC_setting$x_i)


  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = General_KOM_Filter(Ode_Traj_thin_new,MCMC_obj$par[3:(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+2)],
                              MCMC_setting$gridsize,MCMC_setting$x_r,MCMC_setting$x_i)


  LatentTraj_new =cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                          Ode_Traj_coarse_new[,2:3])
  logMultiNorm_new = log_like_traj_general(LatentTraj_new, Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)


  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),
                                 betafs(Ode_Traj_coarse_new[,1],MCMC_obj$par[3:(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+2)], MCMC_setting$x_r,MCMC_setting$x_i),
                                 MCMC_setting$t_correct,
                                 MCMC_setting$gridsize)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[5],MCMC_setting$gridsize)
  }
  if (is.nan(logMultiNorm_new)) {
    logMultiNorm_new = -Inf

  }
  a = dlnorm(alpha1_new,MCMC_setting$b1,MCMC_setting$a1,log = T) - 2*log(1+alpha1_new) + coalLog_new + #dgamma(log(alpha2_new),MCMC_setting$b2,MCMC_setting$a2,log = T) +
    logMultiNorm_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog -
    ( MCMC_obj$LogAlpha1) +  2*log(1+alpha1)

  #print(logMultiNorm_new-log_like_traj2(MCMC_obj$LatentTraj,MCMC_setting$times,log(state_new),MCMC_obj$par[3],MCMC_obj$par[4],MCMC_setting$gridsize,MCMC_setting$t_correct ))

  #print(logMultiNorm_new - MCMC_obj$logMultiNorm)
  if(is.na(a)){a = - Inf}
  # print(c(logMultiNorm_new,MCMC_obj$logMultiNorm,dgamma(log(alpha1_new),60,MCMC_setting$a1,log = T), dgamma(log(alpha2_new),60,MCMC_setting$a2,log = T)))
  AR = 0
  if (log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[1] = state_new[1]
    MCMC_obj$par[2] = state_new[2]
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$FT = FT_new
    MCMC_obj$LogAlpha1 = dlnorm(alpha1_new,MCMC_setting$b1,MCMC_setting$a1,log = T)
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}



##############
updateR0gamma_general = function(MCMC_obj, MCMC_setting, i){
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  R0_new = MCMC_obj$par[3] + runif(1,-MCMC_setting$ps1, MCMC_setting$ps1)
  if(R0_new <1 || R0_new > 10){
    return(list(MCMC_obj = MCMC_obj, AR = 0))
  }
  gamma_new = MCMC_obj$par[4] * exp(runif(1,-MCMC_setting$ps2,MCMC_setting$ps2))


  param_new = c(R0_new, gamma_new, MCMC_obj$par[5:(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+2)])


  Ode_Traj_thin_new <- General_ODE_rk45(MCMC_obj$par[1:2],MCMC_setting$times,param_new,MCMC_setting$x_r,MCMC_setting$x_i)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = General_KOM_Filter(Ode_Traj_thin_new,param_new,
                              MCMC_setting$gridsize,MCMC_setting$x_r,MCMC_setting$x_i)

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                           Ode_Traj_coarse_new[,2:3])
  logMultiNorm_new = log_like_traj_general(LatentTraj_new, Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)
  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),
                                 betafs(Ode_Traj_coarse_new[,1],param_new, MCMC_setting$x_r,MCMC_setting$x_i),
                                 MCMC_setting$t_correct,
                                 MCMC_setting$gridsize)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[5],MCMC_setting$gridsize)

  }

  if (is.nan(logMultiNorm_new)) {
    a = -Inf
  }else{
    a = dlnorm(gamma_new,MCMC_setting$c1,MCMC_setting$c2,log = T) + log(gamma_new) + logMultiNorm_new + coalLog_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog -
      MCMC_obj$LogS2 - log(MCMC_obj$par[4])
  }
  AR = 0
  if(is.na(a)){
    AR = 0
    if(is.na(coalLog_new)){
      print("NA appears in likelihood R0 gamma")
    }else{
      print("NA appears in S1 S2")
    }
  }else if(log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[3] = R0_new
    MCMC_obj$par[4] = gamma_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$LogS2 = dlnorm(gamma_new,MCMC_setting$c1,MCMC_setting$c2,log = T)
    MCMC_obj$FT = FT_new
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}
############

updateR0_general = function(MCMC_obj, MCMC_setting, i){
  R0_new = MCMC_obj$par[3] + runif(1,-MCMC_setting$ps1, MCMC_setting$ps1)
  if(R0_new <1 || R0_new > 10){
    return(list(MCMC_obj = MCMC_obj, AR = 0))
  }


  param_new = c(R0_new, MCMC_obj$par[4], MCMC_obj$par[5:(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+2)])


  Ode_Traj_thin_new <- General_ODE_rk45(MCMC_obj$par[1:2],MCMC_setting$times,param_new,MCMC_setting$x_r,MCMC_setting$x_i)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = General_KOM_Filter(Ode_Traj_thin_new,param_new,
                              MCMC_setting$gridsize,MCMC_setting$x_r,MCMC_setting$x_i)

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                           Ode_Traj_coarse_new[,2:3])

  logMultiNorm_new = log_like_traj_general(LatentTraj_new, Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)
  betaN_new = betafs(Ode_Traj_coarse_new[,1],param_new, MCMC_setting$x_r,MCMC_setting$x_i)
  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),
                                 betaN_new,
                                 MCMC_setting$t_correct,
                                 MCMC_setting$gridsize)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[5],MCMC_setting$gridsize)

  }

  if (is.nan(logMultiNorm_new)) {
    a = -Inf
  }else{
    #a = logMultiNorm_new + coalLog_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog
    a = coalLog_new - MCMC_obj$coalLog
  }
  # print(theta2_new)
  AR = 0
  if(is.na(a)){
    AR = 0
    if(is.na(coalLog_new)){
      print("NA appears in likelihood R0")
      # print(LatentTraj_new)
    }else{
      print("NA appears in S1 S2")
    }
  }else if(log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[3] = R0_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$FT = FT_new
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}





##################
updategamma_general = function(MCMC_obj, MCMC_setting, i){
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  gamma_new = MCMC_obj$par[4] * exp(runif(1,-MCMC_setting$ps2,MCMC_setting$ps2))


  param_new = c(MCMC_obj$par[3], gamma_new, MCMC_obj$par[5:(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+2)])


  Ode_Traj_thin_new <- General_ODE_rk45(MCMC_obj$par[1:2],MCMC_setting$times,param_new,MCMC_setting$x_r,MCMC_setting$x_i)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = General_KOM_Filter(Ode_Traj_thin_new,param_new,
                              MCMC_setting$gridsize,MCMC_setting$x_r,MCMC_setting$x_i)

  LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                           Ode_Traj_coarse_new[,2:3])
  logMultiNorm_new = log_like_traj_general(LatentTraj_new, Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)
  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),
                                 betafs(Ode_Traj_coarse_new[,1],param_new, MCMC_setting$x_r,MCMC_setting$x_i),
                                 MCMC_setting$t_correct,
                                 MCMC_setting$gridsize)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[5],MCMC_setting$gridsize)

  }

  if (is.nan(logMultiNorm_new)) {
    a = -Inf
  }else{
    a = dlnorm(gamma_new,MCMC_setting$c1,MCMC_setting$c2,log = T) + log(gamma_new) + coalLog_new - MCMC_obj$coalLog -
       MCMC_obj$LogS2 - log(MCMC_obj$par[4])
  }
  # print(theta2_new)
  AR = 0
  if(is.na(a)){
    AR = 0
    if(is.na(coalLog_new)){
      print("NA appears in  gamma")
      # print(LatentTraj_new)
    }else{
      print("NA appears in S1 S2")
    }
  }else if(log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[4] = gamma_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$LogS2 = dlnorm(gamma_new,MCMC_setting$c1,MCMC_setting$c2,log = T)
    MCMC_obj$FT = FT_new
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}



###############
updateLambdaGeneral = function(MCMC_obj,MCMC_setting, i){
  lambda_new = MCMC_obj$par[5] * exp(runif(1,-0.3,0.3))
  coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(MCMC_obj$LatentTraj),MCMC_setting$t_correct,lambda_new,MCMC_setting$gridsize)
  a = coalLog_new - MCMC_obj$coalLog + dnorm(lambda_new,MCMC_setting$d1,MCMC_setting$d2,log = T) -
    MCMC_obj$LogLambda
  AR = 0

  #print(coalLog_new - MCMC_obj$coalLog + dgamma(log(lambda_new),MCMC_setting$d1,MCMC_setting$d2,log = T) -
  #       MCMC_obj$LogLambda)
  if(is.nan(a)){
    print(coalLog_new)
  }else if(log(runif(1,0,1)) < a){
    # rec[i,5] = 1
    AR = 1
    MCMC_obj$LogLambda = dnorm(lambda_new,MCMC_setting$d1,MCMC_setting$d2,log = T)
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$par[5] = lambda_new
  }

  return(list(MCMC_obj = MCMC_obj, AR = AR))
}

update_hyper = function(MCMC_obj, MCMC_setting, i){
  p = MCMC_setting$p
  x_i = MCMC_setting$x_i
  param_id = (p+1):(p + x_i[1] + x_i[2] + 1)
  param = MCMC_obj$par[param_id]
  hyper_scale = param[x_i[1] + x_i[2] + 1]
  u = runif(1, - MCMC_setting$proposal$hyper_prop, MCMC_setting$proposal$hyper_prop)
  hyper_scale_new = hyper_scale * exp(u)

  param_new = param
  param_new[x_i[1] + x_i[2] + 1] = hyper_scale_new
  param_new[(x_i[2] + 1):(x_i[1] + x_i[2])] = exp(log(param_new[(x_i[2] + 1):(x_i[1] + x_i[2])]) *
                                                      hyper_scale / hyper_scale_new)

  prior_proposal_offset = dgamma(hyper_scale_new,MCMC_setting$prior$hyper_pr[1], MCMC_setting$prior$hyper_pr[2], log = T) -
      dgamma(hyper_scale, MCMC_setting$prior$hyper_pr[1], MCMC_setting$prior$hyper_pr[2], log = T) + u

  update_res = Update_Param(param_new, MCMC_obj$par[1:p],MCMC_setting$times, MCMC_obj$OriginTraj,
                 MCMC_setting$x_r, MCMC_setting$x_i, MCMC_setting$Init, MCMC_setting$gridsize, MCMC_obj$coalLog, prior_proposal_offset,
                 MCMC_setting$t_correct, model = MCMC_setting$model,
                 volz = MCMC_setting$likelihood == "volz")

  if(update_res$accept){
      MCMC_obj$par[param_id] = param_new
      MCMC_obj$FT = update_res$FT_new
      MCMC_obj$Ode_Traj_coarse = update_res$Ode
      MCMC_obj$betaN = update_res$betaN
      MCMC_obj$coalLog = update_res$coalLog
      MCMC_obj$LatentTraj = update_res$LatentTraj
      MCMC_obj$LogHyper =  dgamma(hyper_scale_new,MCMC_setting$prior$hyper_pr[1], MCMC_setting$prior$hyper_pr[2], log = T)
  }
  return(list(MCMC_obj = MCMC_obj))
}


update_ChangePoint_ESlice = function(MCMC_obj, MCMC_setting,i){
  p = MCMC_setting$p
  param_id = (p+1):(p + MCMC_setting$x_i[1] + MCMC_setting$x_i[2]+1)
  param = MCMC_obj$par[param_id]
  ESlice_Result = ESlice_change_points(param, MCMC_obj$par[1:p],MCMC_setting$times, MCMC_obj$OriginTraj,
                                       MCMC_setting$x_r, MCMC_setting$x_i, MCMC_setting$Init, MCMC_setting$gridsize, MCMC_obj$coalLog,
                                       MCMC_setting$t_correct, model = MCMC_setting$model,
                                       volz = MCMC_setting$likelihood == "volz")

  MCMC_obj$par[param_id] = ESlice_Result$param[1:length(param_id)]
  MCMC_obj$LatentTraj = ESlice_Result$LatentTraj
  MCMC_obj$betaN = ESlice_Result$betaN
  MCMC_obj$FT = ESlice_Result$FT
  MCMC_obj$coalLog = ESlice_Result$CoalLog
  MCMC_obj$Ode_Traj_coarse = ESlice_Result$OdeTraj
  return(MCMC_obj)
}

update_ChangePoint_general = function(MCMC_obj, MCMC_setting, i){
  chpxid =  (MCMC_setting$x_i[2]+3):(sum(MCMC_setting$x_i)+2)
  for( i in chpxid){
    newpara = MCMC_obj$par[c(3,4,5,chpxid)]
    newpara[i-2] = MCMC_obj$par[i] * exp(runif(1,-0.1,0.1))

    Ode_Traj_thin_new <- General_ODE_rk45(MCMC_obj$par[1:2],MCMC_setting$times,newpara,MCMC_setting$x_r,MCMC_setting$x_i)
    Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

    FT_new = General_KOM_Filter(Ode_Traj_thin_new,newpara,
                                MCMC_setting$gridsize,MCMC_setting$x_r,MCMC_setting$x_i)

    LatentTraj_new = cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                             Ode_Traj_coarse_new[,2:3])
    logMultiNorm_new = log_like_traj_general(LatentTraj_new, Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

    if(MCMC_setting$likelihood == "volz"){
      coalLog_new = volz_loglik_nh(MCMC_setting$Init, LogTraj(LatentTraj_new),
                                   betafs(Ode_Traj_coarse_new[,1],newpara, MCMC_setting$x_r,MCMC_setting$x_i),
                                   MCMC_setting$t_correct,
                                   MCMC_setting$gridsize)
    }else{
      coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                                MCMC_obj$par[5],MCMC_setting$gridsize)
    }

    if (is.nan(logMultiNorm_new)) {
      a = -Inf
    }else{
      #a = dlnorm(newpara[i-2],0,MCMC_setting$chpr) + logMultiNorm_new + coalLog_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog - dlnorm(MCMC_obj$par[i],0,MCMC_setting$chpr)
      a = dlnorm(newpara[i-2],0,MCMC_setting$chpr) + coalLog_new - MCMC_obj$coalLog - dlnorm(MCMC_obj$par[i],0,MCMC_setting$chpr)
      }
    # print(theta2_new)
    AR = 0
    if(is.na(a)){
      AR = 0
      if(is.na(coalLog_new)){
        print("NA appears in likelihood changepoint")
        # print(LatentTraj_new)
      }else{
        print("NA appears in changpoint update")
      }
    }else if(log(runif(1,0,1)) < a) {
      AR = 1
      MCMC_obj$par[i] = newpara[i-2]
      MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
      MCMC_obj$logMultiNorm = logMultiNorm_new
      MCMC_obj$FT = FT_new
      MCMC_obj$coalLog = coalLog_new
      MCMC_obj$LatentTraj = LatentTraj_new
    }
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}

updateTraj_general = function(MCMC_obj,MCMC_setting,i){
  #print(c(MCMC_obj$par,MCMC_obj$coalLog + MCMC_obj$logMultiNorm))
  #print(MCMC_obj$par[3])
  #print(betaDyn(MCMC_obj$par[3],MCMC_obj$par[6],MCMC_obj$LatentTraj[,1]))
  MCMC_obj$LatentTraj = ESlice_general(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,
                                      MCMC_obj$FT,MCMC_obj$par[1:2], MCMC_setting$Init,betaN = betafs(MCMC_obj$LatentTraj[,1],MCMC_obj$par[3:(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+2)], MCMC_setting$x_r,MCMC_setting$x_i),
                                      MCMC_setting$t_correct,lambda = MCMC_obj$par[5],
                                      reps = MCMC_setting$reps,MCMC_setting$gridsize,
                                      volz = (MCMC_setting$likelihood == "volz"))
  #q =  MCMC_obj$logMultiNorm
  MCMC_obj$logMultiNorm = log_like_traj_general(MCMC_obj$LatentTraj, MCMC_obj$Ode_Traj_coarse,
                                                MCMC_obj$FT,MCMC_setting$gridsize,MCMC_setting$t_correct)

  if(MCMC_setting$likelihood == "volz"){
    MCMC_obj$coalLog = volz_loglik_nh(MCMC_setting$Init, LogTraj(MCMC_obj$LatentTraj),
                                      betafs(MCMC_obj$LatentTraj[,1],MCMC_obj$par[3:(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+2)], MCMC_setting$x_r,MCMC_setting$x_i),
                                      MCMC_setting$t_correct,
                                      MCMC_setting$gridsize)
  }else{
    MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,LogTraj(MCMC_obj$LatentTraj ),MCMC_setting$t_correct,
                                   MCMC_obj$par[5],MCMC_setting$gridsize)
  }

  return(list(MCMC_obj=MCMC_obj))
}


MCMC_setup_general = function(coal_obs,times,t_correct,N,gridsize=50,niter = 1000,burn = 500,thin = 5, changetime,DEMS = c("E", "I"),
                              prior=list(pop_pr=c(1,10), R0_pr=c(1,7), mu_pr = c(3,0.2), gamma_pr = c(3,0.2), ch_pr = 1,hyper_pr=c(0.01,0.01)),
                              proposal = list(pop_prop = 1, R0_prop = c(0.01), mu_prop=0.1, gamma_prop = 0.2, ch_prop=0.05),
                              control = list(), likelihood = "volz", model = "SIR", Index = c(0,1), nparam = 2,PCOV = NULL){

  gridset = seq(1,length(times),by=gridsize)
  grid = times[gridset]
  Init = coal_lik_init(coal_obs$samp_times, coal_obs$n_sampled, coal_obs$coal_times, grid)
  if(is.null(t_correct)){
    t_correct = max(coal_obs$coal)
  }
  x_i = c(length(changetime),nparam,Index[1],Index[2])
  MCMC_setting = list(Init = Init, times = times,t_correct = t_correct,x_r = c(N,changetime),
                      gridsize = gridsize,gridset = gridset, niter = niter,burn = burn,thin = thin,x_i = x_i,
                      prior = prior, proposal = proposal, control = control, p = length(DEMS),
                      reps = 1, likelihood = likelihood, model = model,PCOV = PCOV)
  cat("MCMC set up ready \n")

  return(MCMC_setting)
}



MCMC_initialize_general = function(MCMC_setting){

  logMultiNorm = NaN
  coalLog = NaN

  ########
  priorLog = numeric(MCMC_setting$p+MCMC_setting$x_i[2])
  while(is.nan(logMultiNorm)||is.nan(coalLog)){
    state = numeric(MCMC_setting$p)
    state[1] = MCMC_setting$x_r[1]
    for(i in 2:MCMC_setting$p){

      if(is.null(MCMC_setting$control$alpha)){
        state[i] = rlnorm(1, MCMC_setting$prior$pop_pr[2*i-3],MCMC_setting$prior$pop_pr[2 * i-2])
      }else{
        state[i] = MCMC_setting$control$alpha[i-1] * MCMC_setting$x_r[1]
      }
      priorLog[i] = dnorm(log(state[i]), MCMC_setting$prior$pop_pr[2 * i - 3],MCMC_setting$prior$pop_pr[2 * i - 2],log = T)
    }

    if(is.null(MCMC_setting$control$R0)){
      R0 = runif(1,MCMC_setting$prior$R0_pr[1],MCMC_setting$prior$R0_pr[2])
    }else{
      R0 = MCMC_setting$control$R0
    }

    if(MCMC_setting$model == "SEIR" || MCMC_setting$model == "SEIR2"){
      if(is.null(MCMC_setting$control$mu)){
        mu = exp(rnorm(1,MCMC_setting$prior$mu_pr[1],MCMC_setting$prior$mu_pr[2]))
      }else{
        mu = MCMC_setting$control$mu
      }
      priorLog[MCMC_setting$p+2] = dnorm(log(mu),MCMC_setting$prior$mu_pr[1],MCMC_setting$prior$mu_pr[2],log = T)
    }

    if(is.null(MCMC_setting$control$gamma)){
      gamma = exp(rnorm(1,MCMC_setting$prior$gamma_pr[1], MCMC_setting$prior$gamma_pr[2]))
    }else{
      gamma = MCMC_setting$control$gamma
    }
    priorLog[MCMC_setting$p+MCMC_setting$x_i[4]+1] = dnorm(log(gamma),MCMC_setting$prior$gamma_pr[1],MCMC_setting$prior$gamma_pr[2],log = T)

    ch=c()
    if(is.null(MCMC_setting$control$hyper)){
      hyper = rgamma(1, MCMC_setting$prior$hyper_pr[1],MCMC_setting$prior$hyper_pr[1])
    }else{
      hyper = MCMC_setting$control$hyper
    }
    if(length(MCMC_setting$x_r) > 1){
      if(is.null(MCMC_setting$control$ch)){
        ch = rlnorm(MCMC_setting$x_i[1],0,1/hyper)
      }else{
        ch = MCMC_setting$control$ch
      }
    }
    if(MCMC_setting$model == "SEIR" || MCMC_setting$model == "SEIR2"){
      param = c(R0,mu,gamma,ch,hyper)
    }else if(MCMC_setting$model == "SIR"){
      param = c(R0, gamma,ch,hyper)
    }

    paramlist = New_Param_List(param, state, MCMC_setting$gridsize, MCMC_setting$times, MCMC_setting$x_r,
                               MCMC_setting$x_i, model = MCMC_setting$model)

   # Ode_Traj_thin = ODE_rk45(state, MCMC_setting$times, param, MCMC_setting$x_r, MCMC_setting$x_i,model = MCMC_setting$model)

    Ode_Traj_coarse = paramlist$Ode
    FT = paramlist$FT
    betaN = paramlist$betaN
    # FT = KF_param_chol(Ode_Traj_thin,param,MCMC_setting$gridsize, MCMC_setting$x_r,MCMC_setting$x_i,model = MCMC_setting$model)

    #betaN = beta
    if(is.null(MCMC_setting$control$traj)){
      Latent = Traj_sim_general_noncentral(Ode_Traj_coarse,FT,MCMC_setting$t_correct)
      #LatentTraj = Latent$SimuTraj
      LatentTraj = Latent$SimuTraj
      OriginTraj = Latent$OriginTraj
      logMultiNorm = Latent$logMultiNorm
      logOrigin = Latent$logOrigin
    }else{
      OriginTraj = MCMC_setting$control$traj
      LatentTraj = TransformTraj(Ode_Traj_coarse, OriginTraj,FT)
      logMultiNorm = log_like_traj_general_adjust(OriginTraj,Ode_Traj_coarse,FT,MCMC_setting$gridsize, MCMC_setting$t_correct)
      logOrigin = - sum(OriginTraj * OriginTraj) / 2
      if( sum(abs(LatentTraj[1,2:(MCMC_setting$p+1)] - state)) > 1){
        print("not consistent")
      }
      logMultiNorm = log_like_traj_general_adjust(LatentTraj,Ode_Traj_coarse,
                                                  FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
    }

    #coalLog = Structural_Coal_lik(MCMC_setting$Init_Detail,LatentTraj,param, MCMC_setting$x_r, MCMC_setting$x_i,
     #                             model = "SEIR2")
    coalLog = volz_loglik_nh2(MCMC_setting$Init,LatentTraj,
                              betaN, MCMC_setting$t_correct, MCMC_setting$x_i[3:4])

    print(paste("coalescent likelihood after initialization ", coalLog))

    if(!is.nan((coalLog))){
      coalLog = ifelse(coalLog <= - 10000000,NaN, coalLog)
    }
    plot(LatentTraj[,1],LatentTraj[,3],type="l",xlab = "time",col="blue", lwd = 2)
    #lines(LatentTraj[,1],LatentTraj[,3],col="red", lwd = 2)
  }


  chprob = 0

  for(i in 1:MCMC_setting$x_i[1]){
    chprob = chprob + dnorm(log(param[MCMC_setting$x_i[2] + i]), 0, 1/hyper,log = T)
  }
  priorLog[MCMC_setting$p + MCMC_setting$x_i[2] +1] = chprob
  #LogGamma = dlnorm(gamma,MCMC_setting$prior$gamma_pr[1],MCMC_setting$prior$gamma_pr[2],log = T)
  #LogE = dlnorm(mu, MCMC_setting$prior$mu_pr[1], MCMC_setting$prior$mu_pr[2],log = T)
  plot(LatentTraj[,1],LatentTraj[,3],type="l",xlab = "time",col="blue", lwd = 2)

  paras =  c(state,param)
  MCMC_obj = list(par = paras,LatentTraj = LatentTraj, logMultiNorm = logMultiNorm,p = MCMC_setting$p,
                  Ode_Traj_coarse = Ode_Traj_coarse, FT = FT, coalLog = coalLog, OriginTraj = OriginTraj,logOrigin = logOrigin,
                  priorLog = priorLog, betaN = betaN)
  ##########
  cat("Initialize MCMC \n")
  print(paste("population size = ", MCMC_setting$N))
  print(paste(paras))
  return(MCMC_obj)
}


updateR0gamma_general_NC = function(MCMC_obj, MCMC_setting, i){
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  p = MCMC_obj$p
  idx1 = MCMC_setting$x_i[3] + 1
  idx2 = MCMC_setting$x_i[4] + 1
  #gamma_new = MCMC_obj$par[p+idx] * exp(runif(1,-MCMC_setting$proposal$gamma_prop, MCMC_setting$proposal$gamma_prop))
  R0logG = mvrnorm(1,c(MCMC_obj$par[p + idx1],log(MCMC_obj$par[p + idx2])),MCMC_setting$PCOV)
  param_new = MCMC_obj$par[(p+1):(MCMC_setting$x_i[1] + MCMC_setting$x_i[2]+p)]
  param_new[idx1] = R0logG[1]
  param_new[idx2] = exp(R0logG[2])

  Ode_Traj_thin_new <- ODE_rk45(MCMC_obj$par[1:p],MCMC_setting$times,param_new,MCMC_setting$x_r,MCMC_setting$x_i, model = MCMC_setting$model)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = KF_param_chol(Ode_Traj_thin_new,param_new,
                         MCMC_setting$gridsize,MCMC_setting$x_r,MCMC_setting$x_i,model = MCMC_setting$model)

  LatentTraj_new = TransformTraj(Ode_Traj_coarse_new, MCMC_obj$OriginTraj, FT_new)

  logMultiNorm_new = log_like_traj_general_adjust(LatentTraj_new, Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)
  betaN_new = betaTs(param_new, Ode_Traj_coarse_new[,1], MCMC_setting$x_r,MCMC_setting$x_i)

  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh2(MCMC_setting$Init, LatentTraj_new,
                                  betaN_new,
                                  MCMC_setting$t_correct,MCMC_setting$x_i[3:4])
  }else if(MCMC_setting$likelihood == "structural"){
    coalLog_new = Structural_Coal_lik(MCMC_setting$Init_Detail, LatentTraj = LatentTraj_new, param = param_new,
                                      MCMC_setting$x_r, MCMC_setting$x_i,model = MCMC_setting$model)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[5],MCMC_setting$gridsize)

  }

  if (is.nan(logMultiNorm_new)) {
    a = -Inf
  }else{
    #a = logMultiNorm_new - MCMC_obj$logMultiNorm + dlnorm(gamma_new,MCMC_setting$b1,MCMC_setting$b2,log = T) -
    #  dlnorm(MCMC_obj$par[p+idx],MCMC_setting$b1,MCMC_setting$b2,log = T) + coalLog_new - MCMC_obj$coalLog
    a = dnorm(R0logG[2],MCMC_setting$prior$gamma_pr[1],MCMC_setting$prior$gamma_pr[2],log = T) -
      MCMC_obj$priorLog[MCMC_setting$p + 1 + MCMC_setting$x_i[4]] + coalLog_new - MCMC_obj$coalLog
  }

  AR = 0

  if(is.na(a)){
    AR = 0
    if(is.na(coalLog_new)){
      print("NA appears in  gamma")
    }else{
      print("NA appears in S1 S2")
    }
  }else if(log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[p+idx1] = R0logG[1]
    MCMC_obj$par[p+idx2] = exp(R0logG[2])
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$priorLog[MCMC_setting$p + 1 + MCMC_setting$x_i[4]] = dnorm(R0logG[2], MCMC_setting$prior$gamma_pr[1],MCMC_setting$prior$gamma_pr[2],log = T)
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$FT = FT_new
    MCMC_obj$LatentTraj = LatentTraj_new
    MCMC_obj$betaN = betaN_new
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}




update_All_general_NC = function(MCMC_obj, MCMC_setting, i){
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  if(i %% 1000 == 0){
    print(MCMC_setting$PCOV)
  }
  p = MCMC_obj$p
  idx1 = MCMC_setting$x_i[3] + 1
  idx2 = MCMC_setting$x_i[4] + 1
  otherId = (p + idx1 + 1):(p + MCMC_setting$x_i[1] + MCMC_setting$x_i[2])
  changePid = (p + MCMC_setting$x_i[2] + 1):(p + MCMC_setting$x_i[1] + MCMC_setting$x_i[2])
  #idx2 = MCMC_setting$x_i[4] + 1
  #gamma_new = MCMC_obj$par[p+idx] * exp(runif(1,-MCMC_setting$proposal$gamma_prop, MCMC_setting$proposal$gamma_prop))
  R0logGlogch = mvrnorm(1,c(MCMC_obj$par[p + idx1],log(MCMC_obj$par[otherId])),
                            MCMC_setting$PCOV)

  param_new = MCMC_obj$par[(p+1):(MCMC_setting$x_i[1] + MCMC_setting$x_i[2]+p)]
  param_new[idx1] = R0logGlogch[1]
  param_new[-idx1] = exp(R0logGlogch[-1])
  if(param_new[2] > 100){
    print(MCMC_obj$par)
    print(MCMC_setting$PCOV)
  }
  Ode_Traj_thin_new <- ODE_rk45(MCMC_obj$par[1:p],MCMC_setting$times,param_new,MCMC_setting$x_r,MCMC_setting$x_i, model = MCMC_setting$model)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = KF_param_chol(Ode_Traj_thin_new,param_new,
                         MCMC_setting$gridsize,MCMC_setting$x_r,MCMC_setting$x_i,model = MCMC_setting$model)

  LatentTraj_new = TransformTraj(Ode_Traj_coarse_new, MCMC_obj$OriginTraj, FT_new)

  logMultiNorm_new = log_like_traj_general_adjust(LatentTraj_new, Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh2(MCMC_setting$Init, LatentTraj_new,
                                  betaTs(param_new, Ode_Traj_coarse_new[,1], MCMC_setting$x_r,MCMC_setting$x_i),
                                  MCMC_setting$t_correct,MCMC_setting$x_i[3:4])
  }else if(MCMC_setting$likelihood == "structural"){
    coalLog_new = Structural_Coal_lik(MCMC_setting$Init_Detail, LatentTraj = LatentTraj_new, param = param_new,
                                      MCMC_setting$x_r, MCMC_setting$x_i,model = MCMC_setting$model)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[5],MCMC_setting$gridsize)

  }
  ChProb = 0
  a = 0
  if (is.nan(logMultiNorm_new)) {
    a = -Inf
  }else{
    #a = logMultiNorm_new - MCMC_obj$logMultiNorm + dlnorm(gamma_new,MCMC_setting$b1,MCMC_setting$b2,log = T) -
    #  dlnorm(MCMC_obj$par[p+idx],MCMC_setting$b1,MCMC_setting$b2,log = T) + coalLog_new - MCMC_obj$coalLog
    # for SIR model only
    ChProb = sum(dnorm(R0logGlogch[changePid-p],0, MCMC_setting$prior$ch_pr,log = T))
    a = a + ChProb - MCMC_obj$priorLog[MCMC_setting$p + MCMC_setting$x_i[2] + 1]
    # SIR model only
    a = a + dnorm(R0logGlogch[2],MCMC_setting$prior$gamma_pr[1],MCMC_setting$prior$gamma_pr[2],log = T) -
      MCMC_obj$priorLog[MCMC_setting$p + 1 + MCMC_setting$x_i[4]] + coalLog_new - MCMC_obj$coalLog
  }

  AR = 0

  if(is.na(a)){
    AR = 0
    if(is.na(coalLog_new)){
      print("NA appears in  gamma")
    }else{
      print("NA appears in S1 S2")
    }
  }else if(log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[p+idx1] = R0logGlogch[1]
    MCMC_obj$par[otherId] = exp(R0logGlogch[-1])
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    if(MCMC_setting$model == "SEIR" || MCMC_setting$model == "SEIR2"){
      MCMC_obj$priorLog[MCMC_setting$p + 2] = dnorm(R0logGlogch[2], MCMC_setting$prior$mu_pr[1],MCMC_setting$prior$mu_pr[2],log = T)
    }
    MCMC_obj$priorLog[p + idx2] = dnorm(R0logGlogch[idx2], MCMC_setting$prior$gamma_pr[1],MCMC_setting$prior$gamma_pr[2],log = T)
    MCMC_obj$priorLog[MCMC_setting$p + 1 + MCMC_setting$x_i[2]] = ChProb
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$FT = FT_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}

update_All_general_NC_Adaptive = function(MCMC_obj, MCMC_setting, i, beta = 0.95, M1, M2){
  p = MCMC_obj$p
  nparam = MCMC_setting$x_i[1] + MCMC_setting$x_i[2]
  #param_new = MCMC_obj$par[(p+1):(MCMC_setting$x_i[1] + MCMC_setting$x_i[2]+p)]
 # if(i == 200){
  #  SigmaN = (2.38^2) * cov(cbind(params[100:199,p+1], log(params[100:199,(p+2):(p+nparam)]))) / nparam
 #   MCMC_setting$PCOV = beta * SigmaN + (1 - beta) * diag(rep(1,nparam)) * 0.01 / nparam
    #}else if(i > 2000 && i %% 1000 == 0){
    #  SigmaN = (2.38^2) * params[1000:i, (p+1):(p+nparam)] / nparam
    # MCMC_setting$PCOV = beta * SigmaN + (1 - beta) * diag(rep(1,nparam)) * 0.01 / nparam
    #}#
#  }else{
    # MCMC_obj$M2 = MCMC_obj$M2 * ((i-2) / (i-1)) + param_new %*% t(param_new) / (i-1)
    # MCMC_obj$M1 = MCMC_obj$M1 * ((i - 1) / i) + param_new / i
    SigmaN = M2 - (M1) %*% t(M1)
    #SigmaN = (2.38^2) * cov(cbind(params[100:199,p+1], log(params[100:199,(p+2):(p+nparam)]))) / nparam
    SigmaN = (2.38^2) * SigmaN / nparam
    MCMC_setting$PCOV = beta * SigmaN + (1 - beta) * diag(rep(1,nparam)) * 0.01 / nparam
#  }

  step = update_All_general_NC(MCMC_obj, MCMC_setting, i)
  return(list(MCMC_obj = step$MCMC_obj, AR = step$AR))
}




General_MCMC = function(coal_obs,times,t_correct,N,gridsize=1000, niter = 1000, burn = 0, thin = 5,changetime,DEMS=c("S","I"),
                        prior=list(pop_pr=c(1,10,1,10), R0_pr=c(1,7), mu_pr = c(3,0.2), gamma_pr = c(3,0.2), ch_pr = 1),
                        proposal = list(pop_prop = 1, R0_prop = c(0.01), mu_prop=0.1, gamma_prop = 0.2, ch_prop=0.05),
                        control = list(), updateVec = c(1,1,1,1,1,1), likelihood = "volz",model = "SIR",
                        Index = c(0,2), nparam=3, joint = F, PCOV = NULL,beta = 0.05, burn1 = 5000){

  MCMC_setting = MCMC_setup_general(coal_obs, times,t_correct,N,gridsize,niter,burn,
                                    thin,changetime, DEMS,prior,proposal,
                                    control,likelihood,model,Index,nparam, PCOV)

  nparam = sum(MCMC_setting$x_i[1:2])
  MCMC_obj = MCMC_initialize_general(MCMC_setting)

  if(MCMC_setting$likelihood == "volz"){
    params = matrix(nrow = niter, ncol = nparam + MCMC_obj$p)
    ARMS = matrix(nrow = niter, ncol =  nparam + MCMC_obj$p)
  }else if(MCMC_setting$likelihood == "structural"){
    params = matrix(nrow = niter, ncol =  nparam + MCMC_obj$p)
    ARMS = matrix(nrow = niter, ncol =  nparam + MCMC_obj$p)
  }else{
    params = matrix(nrow = niter, ncol = sum(MCMC_setting$x_i) + MCMC_obj$p)
    ARMS = matrix(nrow = niter, ncol = sum(MCMC_setting$x_i) + MCMC_obj$p)
  }
  M1 = numeric(sum(MCMC_setting$x_i[1:2]))
  M2 = matrix(rep(0, sum(MCMC_setting$x_i[1:2])^2), ncol = sum(MCMC_setting$x_i[1:2]))
  l = numeric(niter)
  l1 = l
  l2 = l
  l3 = l
  l0 = l
  tjs = array(dim = c(dim(MCMC_obj$LatentTraj),niter))

  for (i in 1:MCMC_setting$niter) {
    if (i %% 100 == 0) {
      print(i)
      print(MCMC_obj$par)
      print(l1[i-1])
      print(l2[i-1])
      print(l3[i-1])
      plot(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,3],type="l")
      lines(MCMC_obj$Ode_Traj_coarse[,1],MCMC_obj$Ode_Traj_coarse[,3],col="red",lty=2)
    }
    ARvec = numeric(dim(ARMS)[2])
      if(updateVec[1] == 1){
        step1 = updateAlphas_General(MCMC_obj,MCMC_setting,i)
        MCMC_obj = step1$MCMC_obj
        ARvec[1] = step1$AR
      }
   # if(i == burn1){
    #  idx = floor(burn1/2) : (burn1-1)
    #  SigmaN = (2.38^2) * cov(cbind(params[idx,MCMC_obj$p + MCMC_setting$x_i[3] + 1],
    #                 log(params[idx,-(1:(MCMC_obj$p+1))]))) / nparam
    #  MCMC_setting$PCOV = beta * SigmaN + (1 - beta) * diag(rep(1,nparam)) * 0.01 / nparam
    #  print(MCMC_setting$PCOV)
  #  }

    if(i == burn1){
      idx = floor(burn1/2) : (burn1-1)
      SigmaN = (2.38^2) * cov(cbind(params[idx,MCMC_obj$p + MCMC_setting$x_i[3] + 1],
                                log(params[idx,MCMC_obj$p + MCMC_setting$x_i[4] + 1])))/2

      MCMC_setting$PCOV = beta * SigmaN + (1 - beta) * diag(rep(1,2)) * 0.01 / 2

    }
    if(i < burn1){
      if(updateVec[2] == 1){
        step2 = updateR0_general_NC(MCMC_obj,MCMC_setting,i)
        # print(c(MCMC_obj$coalLog,MCMC_obj$logMultiNorm))
        MCMC_obj = step2$MCMC_obj
        ARvec[2] = step2$AR
      }
      if(updateVec[3] == 1){
        step3 = updategamma_general_NC(MCMC_obj,MCMC_setting,i)
        # print(c(MCMC_obj$coalLog,MCMC_obj$logMultiNorm))
        MCMC_obj = step3$MCMC_obj
        ARvec[3] = step3$AR
      }

      if(updateVec[5] == 1){
        #MCMC_obj = update_ChangePoint_general_NC(MCMC_obj,MCMC_setting,i)$MCMC_obj
        MCMC_obj = update_ChangePoint_ESlice(MCMC_obj,MCMC_setting,i)
      }
      if(updateVec[7] == 1){
        MCMC_obj = update_hyper(MCMC_obj, MCMC_setting, i)$MCMC_obj
      }
    }else{

      step2 = updateR0gamma_general_NC(MCMC_obj, MCMC_setting, i)
      MCMC_obj = step2$MCMC_obj
      ARvec[2] = step2$AR
      if(updateVec[5] == 1){
        #MCMC_obj = update_ChangePoint_general_NC(MCMC_obj,MCMC_setting,i)$MCMC_obj
        MCMC_obj = update_ChangePoint_ESlice(MCMC_obj,MCMC_setting,i)
      }

      if(updateVec[6] == 1){
        MCMC_obj = updateTraj_general_NC(MCMC_obj,MCMC_setting,i)$MCMC_obj
      }
    }
      # if(updateVec[3] == 1 && updateVec[2] == 1){
      #  stepJoint = updateR0gamma_general(MCMC_obj,MCMC_setting,i)
      #  ARvec[4] = stepJoint$AR
      #  MCMC_obj = stepJoint$MCMC_obj
      #}

     # if(updateVec[4] == 1){
     #   steplambda = updateLambdaGeneral(MCMC_obj,MCMC_setting,i)
     #   ARvec[5] = steplambda$AR
     #   MCMC_obj = steplambda$MCMC_obj
     # }

    #step4 = updateLambda_SIRS(MCMC_obj,MCMC_setting,i)
    #ARvec[8] = step4$AR
    #MCMC_obj = step4$MCMC_obj
    ARMS[i,] = ARvec
    #tjs = abind(tjs,MCMC_obj$LatentTraj,along = 3)
    tjs[,,i] = MCMC_obj$LetentTra
    params[i,] = MCMC_obj$par
    l[i] =  MCMC_obj$priorLog[MCMC_obj$p + MCMC_setting$x_i[4] + 1] #+ MCMC_obj$LogAlpha1 + MCMC_obj$LogAlpha2
    l1[i] = MCMC_obj$coalLog
    l2[i] = sum(MCMC_obj$priorLog) + MCMC_obj$coalLog + MCMC_obj$logOrigin
    l3[i] =  MCMC_obj$logOrigin
    newSample = c(params[i,MCMC_obj$p + MCMC_setting$x_i[3] + 1],
           log(params[i,-(1:(MCMC_obj$p+1))]))
    M1 = M1 * (i - 1) / i + newSample / i
    M2 = M2 * (i - 1) / i + newSample %*% t(newSample) / i
  }
  return(list(par = params,Trajectory = tjs,l=l,l1=l1,l2 = l2,l3=l3,AR = ARMS))
}


#######
####################
updateR0_general_NC = function(MCMC_obj, MCMC_setting, i){
  p = MCMC_obj$p
  R0_new = MCMC_obj$par[p+1] + runif(1,-MCMC_setting$proposal$R0_prop, MCMC_setting$proposal$R0_prop)

  if(R0_new <1 || R0_new > 10){
    return(list(MCMC_obj = MCMC_obj, AR = 0))
  } # reject if out of bound

  param_new = MCMC_obj$par[(p+1):(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+p)]
  param_new[1] = R0_new

  Ode_Traj_thin_new <- ODE_rk45(MCMC_obj$par[1:p],MCMC_setting$times,param_new,MCMC_setting$x_r,MCMC_setting$x_i, model = MCMC_setting$model)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = KF_param_chol(Ode_Traj_thin_new,param_new,
                         MCMC_setting$gridsize,MCMC_setting$x_r,MCMC_setting$x_i,model = MCMC_setting$model)

  LatentTraj_new = TransformTraj(Ode_Traj_coarse_new, MCMC_obj$OriginTraj, FT_new)

  logMultiNorm_new = log_like_traj_general_adjust(LatentTraj_new, Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)
  betaN_new = betaTs(param_new, Ode_Traj_coarse_new[,1],MCMC_setting$x_r,MCMC_setting$x_i)

  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh2(MCMC_setting$Init, LatentTraj_new,
                                  betaN_new,
                                  MCMC_setting$t_correct,MCMC_setting$x_i[3:4])
  }else if(MCMC_setting$likelihood == "structural"){
    coalLog_new = Structural_Coal_lik(MCMC_setting$Init_Detail, LatentTraj = LatentTraj_new, param = param_new,
                                      MCMC_setting$x_r, MCMC_setting$x_i,model = MCMC_setting$model)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[5],MCMC_setting$gridsize)
  }

  if (is.nan(logMultiNorm_new)) {
    a = -Inf
  }else{
    a =  coalLog_new -  MCMC_obj$coalLog
  }
  AR = 0
  #a = logMultiNorm_new - MCMC_obj$logMultiNorm + coalLog_new - MCMC_obj$coalLog
  if(is.na(a)){

    if(is.na(coalLog_new)){
      print("NA appears in likelihood R0")
      # print(LatentTraj_new)
    }else{
      print("NA Trajectory")
    }
  }else if(log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[p+1] = R0_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$FT = FT_new
    MCMC_obj$LatentTraj = LatentTraj_new
    MCMC_obj$betaN = betaN_new
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}



##################
updategamma_general_NC = function(MCMC_obj, MCMC_setting, i){
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  p = MCMC_obj$p
  idx = MCMC_setting$x_i[4]+1
  gamma_new = MCMC_obj$par[p+idx] * exp(runif(1,-MCMC_setting$proposal$gamma_prop, MCMC_setting$proposal$gamma_prop))

  param_new = MCMC_obj$par[(p+1):(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+p)]
  param_new[idx] = gamma_new

  Ode_Traj_thin_new <- ODE_rk45(MCMC_obj$par[1:p],MCMC_setting$times,param_new,MCMC_setting$x_r,MCMC_setting$x_i, model = MCMC_setting$model)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = KF_param_chol(Ode_Traj_thin_new,param_new,
                         MCMC_setting$gridsize,MCMC_setting$x_r,MCMC_setting$x_i,model = MCMC_setting$model)


  LatentTraj_new = TransformTraj(Ode_Traj_coarse_new, MCMC_obj$OriginTraj, FT_new)

  logMultiNorm_new = log_like_traj_general_adjust(LatentTraj_new, Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

  betaN_new = betaTs(param_new, Ode_Traj_coarse_new[,1], MCMC_setting$x_r,MCMC_setting$x_i)
  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh2(MCMC_setting$Init, LatentTraj_new,
                                  betaN_new,
                                  MCMC_setting$t_correct,MCMC_setting$x_i[3:4])
  }else if(MCMC_setting$likelihood == "structural"){
    coalLog_new = Structural_Coal_lik(MCMC_setting$Init_Detail, LatentTraj = LatentTraj_new, param = param_new,
                                      MCMC_setting$x_r, MCMC_setting$x_i,model = MCMC_setting$model)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[5],MCMC_setting$gridsize)

  }

  if (is.nan(logMultiNorm_new)) {
    a = -Inf
  }else{
    #a = logMultiNorm_new - MCMC_obj$logMultiNorm + dlnorm(gamma_new,MCMC_setting$b1,MCMC_setting$b2,log = T) -
    #  dlnorm(MCMC_obj$par[p+idx],MCMC_setting$b1,MCMC_setting$b2,log = T) + coalLog_new - MCMC_obj$coalLog
    a = dnorm(log(gamma_new),MCMC_setting$prior$gamma_pr[1],MCMC_setting$prior$gamma_pr[2],log = T) -
      MCMC_obj$priorLog[MCMC_setting$p + 1 + MCMC_setting$x_i[4]] + coalLog_new - MCMC_obj$coalLog
  }

  AR = 0

  if(is.na(a)){
    AR = 0
    if(is.na(coalLog_new)){
      print("NA appears in  gamma")
    }else{
      print("NA appears in S1 S2")
    }
  }else if(log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[p+idx] = gamma_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$priorLog[MCMC_setting$p + 1 + MCMC_setting$x_i[4]] = dnorm(log(gamma_new),MCMC_setting$prior$gamma_pr[1],MCMC_setting$prior$gamma_pr[2],log = T)
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$FT = FT_new
    MCMC_obj$LatentTraj = LatentTraj_new
    MCMC_obj$betaN = betaN_new
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}


updateE_general_NC = function(MCMC_obj, MCMC_setting, i,idx = 2){
  #s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 3),6)
  p = MCMC_obj$p
  idx = 2
  gamma_new = MCMC_obj$par[p+idx] * exp(runif(1, - MCMC_setting$proposal$mu_prop, MCMC_setting$proposal$mu_prop))

  param_new = MCMC_obj$par[(p+1):(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+p)]
  param_new[idx] = gamma_new

  Ode_Traj_thin_new <- ODE_rk45(MCMC_obj$par[1:p],MCMC_setting$times,param_new,MCMC_setting$x_r,MCMC_setting$x_i, model = MCMC_setting$model)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = KF_param_chol(Ode_Traj_thin_new,param_new,
                         MCMC_setting$gridsize,MCMC_setting$x_r,MCMC_setting$x_i,model = MCMC_setting$model)


  LatentTraj_new = TransformTraj(Ode_Traj_coarse_new, MCMC_obj$OriginTraj, FT_new)

  logMultiNorm_new = log_like_traj_general_adjust(LatentTraj_new, Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

  if(MCMC_setting$likelihood == "volz"){
    coalLog_new = volz_loglik_nh2(MCMC_setting$Init, LatentTraj_new,
                                  betaTs(param_new, Ode_Traj_coarse_new[,1], MCMC_setting$x_r,MCMC_setting$x_i),
                                  MCMC_setting$t_correct,MCMC_setting$x_i[3:4])
  }else if(MCMC_setting$likelihood == "structural"){
    coalLog_new = Structural_Coal_lik(MCMC_setting$Init_Detail, LatentTraj = LatentTraj_new, param = param_new,
                                      MCMC_setting$x_r, MCMC_setting$x_i,model = MCMC_setting$model)
  }else{
    coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                              MCMC_obj$par[5],MCMC_setting$gridsize)
  }

  if (is.nan(logMultiNorm_new)) {
    a = -Inf
  }else{
    #a = logMultiNorm_new - MCMC_obj$logMultiNorm + dlnorm(gamma_new,MCMC_setting$c1,MCMC_setting$c2,log = T) -
    # dlnorm(MCMC_obj$par[p+idx],MCMC_setting$c1,MCMC_setting$c2,log = T) + coalLog_new - MCMC_obj$coalLog
    a = dnorm(log(gamma_new),MCMC_setting$prior$mu_pr[1],MCMC_setting$prior$mu_pr[2],log = T) -
      MCMC_obj$priorLog[MCMC_setting$p + 2] + coalLog_new - MCMC_obj$coalLog
  }

  AR = 0

  if(is.na(a)){
    AR = 0
    if(is.na(coalLog_new)){
      print("NA appears in  gamma")
    }else{
      print("NA appears in S1 S2")
    }
  }else if(log(runif(1,0,1)) < a) {
    AR = 1
    MCMC_obj$par[p+idx] = gamma_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$priorLog[MCMC_setting$p + 2] = dnorm(log(gamma_new),MCMC_setting$prior$mu_pr[1],MCMC_setting$prior$mu_pr[2],log = T)
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$FT = FT_new
    MCMC_obj$LatentTraj = LatentTraj_new
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}





update_ChangePoint_general_NC = function(MCMC_obj, MCMC_setting, i){
  p = MCMC_obj$p
  chpxid =  (MCMC_setting$x_i[2]+p+1):(sum(MCMC_setting$x_i[1:2])+p)
  prob = 0
  for(i in chpxid){
    newpara = MCMC_obj$par[c((p + 1):(MCMC_setting$x_i[2] + p),chpxid)]
    newpara[i-p] = MCMC_obj$par[i] * exp(runif(1,-MCMC_setting$proposal$ch_prop, MCMC_setting$proposal$ch_prop))

    Ode_Traj_thin_new <- ODE_rk45(MCMC_obj$par[1:p],MCMC_setting$times,newpara,MCMC_setting$x_r,MCMC_setting$x_i, model = MCMC_setting$model)

    Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

    FT_new = KF_param_chol(Ode_Traj_thin_new,newpara,
                           MCMC_setting$gridsize,MCMC_setting$x_r,MCMC_setting$x_i,model = MCMC_setting$model)

    LatentTraj_new = TransformTraj(Ode_Traj_coarse_new, MCMC_obj$OriginTraj, FT_new)
    logMultiNorm_new = log_like_traj_general_adjust(LatentTraj_new, Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)


    if(MCMC_setting$likelihood == "volz"){
      coalLog_new = volz_loglik_nh2(MCMC_setting$Init, LatentTraj_new,
                                    betaTs(newpara,Ode_Traj_coarse_new[,1],MCMC_setting$x_r,MCMC_setting$x_i),
                                    MCMC_setting$t_correct,MCMC_setting$x_i[3:4])
    }else if(MCMC_setting$likelihood == "structural"){
      coalLog_new = Structural_Coal_lik(MCMC_setting$Init_Detail, LatentTraj = LatentTraj_new, param = newpara,
                                        MCMC_setting$x_r, MCMC_setting$x_i,model = MCMC_setting$model)
    }else{
      coalLog_new = coal_loglik(MCMC_setting$Init,LogTraj(LatentTraj_new),MCMC_setting$t_correct,
                                MCMC_obj$par[5],MCMC_setting$gridsize)
    }

    if (is.nan(logMultiNorm_new)){
      a = -Inf
    }else{
      #a = dlnorm(newpara[i-p],0,MCMC_setting$chpr,log = T) + logMultiNorm_new + coalLog_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog - dlnorm(MCMC_obj$par[i],0,MCMC_setting$chpr,log = T)
      a = dnorm(log(newpara[i-2]),0,MCMC_setting$prior$ch_pr,log = T) + coalLog_new -
        MCMC_obj$coalLog - dnorm(log(MCMC_obj$par[i]),0,MCMC_setting$prior$ch_pr, log = T)
    }

    AR = 0
    if(is.na(a)){
      AR = 0
      if(is.na(coalLog_new)){
        print("NA appears in likelihood changepoint")
        # print(LatentTraj_new)
      }else{
        print("NA appears in changpoint update trajectory")
      }
    }else if(log(runif(1,0,1)) < a) {
      AR = 1
      MCMC_obj$FT = FT_new
      MCMC_obj$par[i] = newpara[i-p]
      MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
      MCMC_obj$logMultiNorm = logMultiNorm_new
      MCMC_obj$coalLog = coalLog_new
      MCMC_obj$LatentTraj = LatentTraj_new
    }
    prob = prob + dnorm(log(MCMC_obj$par[i]),0,MCMC_setting$prior$ch_pr,log = T)
  }
  MCMC_obj$priorLog[MCMC_setting$p + MCMC_setting$x_i[2] + 1] = prob
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}
