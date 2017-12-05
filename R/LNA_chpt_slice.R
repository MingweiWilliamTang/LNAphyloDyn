
update_Param_joint = function(MCMC_obj, MCMC_setting, method = "admcmc", parId, isLog, priorId, proposeId = NULL){

  if(length(parId) != length(isLog)){
    stop("parameters do not match")
  }
  if(length(parId) != length(priorId)){
    stop("priors do not match")
  }
  d = length(parId)
  p = MCMC_obj$p
  x_i = MCMC_setting$x_i
  par_new = MCMC_obj$par
  initialId = parId[parId <= p]
  paramId = parId[parId > p & parId <= (p + x_i[1] + x_i[2] + 1)]

  initial_new = MCMC_obj$par[1:p]
  param_new = MCMC_obj$par[-(1:p)]
  hyperId = (p + MCMC_setting$x_i[1] + MCMC_setting$x_i[2] + 1)
  if(method == "admcmc"){

    RawTransParam = MCMC_obj$par[parId]
    RawTransParam[isLog == 1] = log(RawTransParam[isLog == 1])
    RawTransParam_new = mvrnorm(1,RawTransParam, MCMC_setting$PCOV)


    prior_proposal_offset = 0
    for(i in 1:length(parId)){
      newdiff = 0

      if(priorId[i] %in% c(1,3,4)){
        pr = MCMC_setting$prior[[ priorId[i] ]]
        newdiff = dnorm(RawTransParam_new[i], pr[1], pr[2], log = T) - dnorm(RawTransParam[i], pr[1], pr[2], log = T)
      }
      if(priorId[i] == 2){
        pr = MCMC_setting$prior[[priorId[i]]]
        newdiff = dunif(RawTransParam_new[i], pr[1], pr[2], log = T) - dunif(RawTransParam[i], pr[1], pr[2], log = T)
      }
      if(priorId[i] == 5){
        pr = MCMC_setting$prior[[priorId[i]]]
        newdiff = dgamma(RawTransParam_new[i], pr[1], pr[2], log = T) - dgamma(RawTransParam[i], pr[1], pr[2], log = T)
      }

      prior_proposal_offset = prior_proposal_offset + newdiff
    }

    RawTransParam_new[isLog == 1] = exp(RawTransParam_new[isLog == 1])

    if(hyperId %in% parId){
      par_new[(p + x_i[2] + 1): (hyperId - 1)] = exp(log(par_new[(p + x_i[2] + 1): (hyperId - 1)]) *
                                                       RawTransParam[d] / RawTransParam_new[d])
    }

    par_new[parId] = RawTransParam_new
    initial_new[initialId] = par_new[initialId]
    param_new = par_new[-(1:p)]

    update_res = Update_Param(param_new, initial_new, MCMC_setting$times, MCMC_obj$OriginTraj,
                              MCMC_setting$x_r, MCMC_setting$x_i, MCMC_setting$Init, MCMC_setting$gridsize, MCMC_obj$coalLog,prior_proposal_offset,
                              MCMC_setting$t_correct, model = MCMC_setting$model,
                              volz = MCMC_setting$likelihood == "volz")

    if(update_res$accept){
      MCMC_obj$par = par_new
      MCMC_obj$FT = update_res$FT_new
      MCMC_obj$Ode_Traj_coarse = update_res$Ode
      MCMC_obj$betaN = update_res$betaN
      MCMC_obj$coalLog = update_res$coalLog
      MCMC_obj$LatentTraj = update_res$LatentTraj
    }

    return(list(MCMC_obj = MCMC_obj))
  }else if(method == "jointProp"){

    if(length(parId) != length(proposeId)){
      stop("propsals do not match")
    }

    RawTransParam = MCMC_obj$par[parId]
    RawTransParam[isLog == 1] = log(RawTransParam[isLog == 1])
    RawTransParam_new = RawTransParam
    prior_proposal_offset = 0
    for(i in 1:length(parId)){
      newdiff = 0

      if(priorId[i] %in% c(1,3,4)){
        pr = MCMC_setting$prior[[ priorId[i] ]]
        po = MCMC_setting$proposal[[proposeId[i]]]
        RawTransParam_new[i] = RawTransParam[i] + runif(1,-po, po)
        newdiff = dnorm(RawTransParam_new[i], pr[1], pr[2], log = T) - dnorm(RawTransParam[i], pr[1], pr[2], log = T)
      }
      if(priorId[i] == 2){
        pr = MCMC_setting$prior[[priorId[i]]]
        po = MCMC_setting$proposal[[proposeId[i]]]
        RawTransParam_new[i] = RawTransParam[i] + runif(1,-po, po)
        newdiff = dunif(RawTransParam_new[i], pr[1], pr[2], log = T) - dunif(RawTransParam[i], pr[1], pr[2], log = T)
      }
      if(priorId[i] == 5){
        pr = MCMC_setting$prior[[priorId[i]]]
        po = MCMC_setting$proposal[[proposeId[i]]]
        u = runif(1,-po, po)
        RawTransParam_new[i] = RawTransParam[i] * exp(u)
        newdiff = dgamma(RawTransParam_new[i], pr[1], pr[2], log = T) - u - dgamma(RawTransParam[i], pr[1], pr[2], log = T)
      }

      prior_proposal_offset = prior_proposal_offset + newdiff
    }

    RawTransParam_new[isLog == 1] = exp(RawTransParam_new[isLog == 1])

    if(hyperId %in% parId){
      par_new[(p + x_i[2] + 1): (hyperId - 1)] = exp(log(par_new[(p + x_i[2] + 1): (hyperId - 1)]) *
                                                       RawTransParam[d] / RawTransParam_new[d])
    }

    par_new[parId] = RawTransParam_new
    initial_new[initialId] = par_new[initialId]
    param_new = par_new[-(1:p)]

    update_res = Update_Param(param_new, initial_new, MCMC_setting$times, MCMC_obj$OriginTraj,
                              MCMC_setting$x_r, MCMC_setting$x_i, MCMC_setting$Init, MCMC_setting$gridsize, MCMC_obj$coalLog,prior_proposal_offset,
                              MCMC_setting$t_correct, model = MCMC_setting$model,
                              volz = MCMC_setting$likelihood == "volz")

    if(update_res$accept){
      MCMC_obj$par = par_new
      MCMC_obj$FT = update_res$FT_new
      MCMC_obj$Ode_Traj_coarse = update_res$Ode
      MCMC_obj$betaN = update_res$betaN
      MCMC_obj$coalLog = update_res$coalLog
      MCMC_obj$LatentTraj = update_res$LatentTraj
    }

    return(list(MCMC_obj = MCMC_obj))
  }

}

check_list_eq = function(list1, list2){
  q1 = (length(list1) == length(list2))
  q2 = (all.equal(unlist(lapply(list1, length)), unlist(lapply(list2,length)),
            check.attributes = FALSE) == TRUE)
  return(q1 && q2)
}


Update_Param_each = function(MCMC_obj, MCMC_setting, parIdlist, isLoglist, priorIdlist, proposeIdlist){

  if(!check_list_eq(parIdlist, isLoglist)){
    stop("parameters do not match")
  }
  if(!check_list_eq(parIdlist, priorIdlist)){
    stop("priors do not match")
  }
  if(!check_list_eq(parIdlist, proposeIdlist)){
    stop("priors do not match")
  }

  d = length(parIdlist)
  p = MCMC_obj$p
  x_i = MCMC_setting$x_i
  hyperId = (p + MCMC_setting$x_i[1] + MCMC_setting$x_i[2] + 1)

  for(i in 1:d){
    par_new = MCMC_obj$par
    initial_new = MCMC_obj$par[1:p]
    param_new = MCMC_obj$par[-(1:p)]
    subd = length(parIdlist[[i]])
    parId = parIdlist[[i]]
    isLog = isLoglist[[i]]
    priorId = priorIdlist[[i]]
    proposeId = proposeIdlist[[i]]

    initialId = parId[parId <= p]
    paramId = parId[parId > p & parId <= (p + x_i[1] + x_i[2] + 1)]

    RawTransParam = MCMC_obj$par[parId]
    RawTransParam[isLog == 1] = log(RawTransParam[isLog == 1])
    RawTransParam_new = RawTransParam
    prior_proposal_offset = 0

    for(i in 1:length(parId)){

      newdiff = 0

      if(priorId[i] %in% c(1,3,4)){
        pr = MCMC_setting$prior[[ priorId[i] ]]
        po = MCMC_setting$proposal[[proposeId[i]]]
        RawTransParam_new[i] = RawTransParam[i] + runif(1,-po, po)
        newdiff = dnorm(RawTransParam_new[i], pr[1], pr[2], log = T) - dnorm(RawTransParam[i], pr[1], pr[2], log = T)
      }
      if(priorId[i] == 2){
        pr = MCMC_setting$prior[[priorId[i]]]
        po = MCMC_setting$proposal[[proposeId[i]]]
        RawTransParam_new[i] = RawTransParam[i] + runif(1,-po, po)
        newdiff = dunif(RawTransParam_new[i], pr[1], pr[2], log = T) - dunif(RawTransParam[i], pr[1], pr[2], log = T)
      }
      if(priorId[i] == 5){
        pr = MCMC_setting$prior[[priorId[i]]]
        po = MCMC_setting$proposal[[proposeId[i]]]
        u = runif(1,-po, po)
        RawTransParam_new[i] = RawTransParam[i] * exp(u)
        newdiff = dgamma(RawTransParam_new[i], pr[1], pr[2], log = T) - u - dgamma(RawTransParam[i], pr[1], pr[2], log = T)
      }

      prior_proposal_offset = prior_proposal_offset + newdiff
    }

    RawTransParam_new[isLog == 1] = exp(RawTransParam_new[isLog == 1])

    if(hyperId %in% parId){
      par_new[(p + x_i[2] + 1): (hyperId - 1)] = exp(log(par_new[(p + x_i[2] + 1): (hyperId - 1)]) *
                                                       RawTransParam[subd] / RawTransParam_new[subd])
    }

    par_new[parId] = RawTransParam_new
    initial_new[initialId] = par_new[initialId]
    param_new = par_new[-(1:p)]

    update_res = Update_Param(param_new, initial_new, MCMC_setting$times, MCMC_obj$OriginTraj,
                              MCMC_setting$x_r, MCMC_setting$x_i, MCMC_setting$Init, MCMC_setting$gridsize, MCMC_obj$coalLog,prior_proposal_offset,
                              MCMC_setting$t_correct, model = MCMC_setting$model,
                              volz = MCMC_setting$likelihood == "volz")


    if(update_res$accept){
      MCMC_obj$par = par_new
      MCMC_obj$FT = update_res$FT_new
      MCMC_obj$Ode_Traj_coarse = update_res$Ode
      MCMC_obj$betaN = update_res$betaN
      MCMC_obj$coalLog = update_res$coalLog
      MCMC_obj$LatentTraj = update_res$LatentTraj
    }

  }
  return(list(MCMC_obj = MCMC_obj))
}



updateTraj_general_NC = function(MCMC_obj,MCMC_setting,i){

  if(MCMC_setting$likelihood == "structural"){

    Res = ESlice_general_NC_Structural(MCMC_obj$OriginTraj, MCMC_obj$Ode_Traj_coarse,
                                       MCMC_obj$FT, MCMC_setting$Init_Detail,
                                       MCMC_obj$par[(MCMC_obj$p+1):(MCMC_obj$p + MCMC_setting$x_i[1] + MCMC_setting$x_i[2])],
                                       MCMC_setting$x_r, MCMC_setting$x_i,
                                       coal_log = MCMC_obj$coalLog,
                                       model = MCMC_setting$model)
    MCMC_obj$coalLog = Res$CoalLog
  }else{

    Res = ESlice_general_NC(MCMC_obj$OriginTraj,MCMC_obj$Ode_Traj_coarse,
                            MCMC_obj$FT, MCMC_obj$par[1:MCMC_obj$p], MCMC_setting$Init,
                            betaN = betaTs(MCMC_obj$par[(MCMC_obj$p+1):(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+MCMC_obj$p)],MCMC_obj$LatentTraj[,1], MCMC_setting$x_r,MCMC_setting$x_i),
                            MCMC_setting$t_correct,lambda = 1,
                            coal_log = MCMC_obj$coalLog, MCMC_setting$gridsize,
                            volz = (MCMC_setting$likelihood == "volz"), model = MCMC_setting$model)


    if(MCMC_setting$likelihood == "volz"){
      MCMC_obj$coalLog = volz_loglik_nh2(MCMC_setting$Init, Res$LatentTraj,
                                         betaTs(MCMC_obj$par[(MCMC_obj$p+1):(MCMC_setting$x_i[1]+MCMC_setting$x_i[2]+MCMC_obj$p)],MCMC_obj$LatentTraj[,1], MCMC_setting$x_r,MCMC_setting$x_i),
                                         MCMC_setting$t_correct,
                                         index = MCMC_setting$x_i[3:4])
    }else{
      MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,LogTraj(Res$LatentTraj ),MCMC_setting$t_correct,
                                     MCMC_obj$par[5],MCMC_setting$gridsize)
    }
  }

  MCMC_obj$LatentTraj = Res$LatentTraj
  MCMC_obj$OriginTraj = Res$OriginTraj
  MCMC_obj$logOrigin = Res$logOrigin


  return(list(MCMC_obj = MCMC_obj))
}


#'
#'
#' nparam number of parameters in infectious disease model

General_MCMC2 = function(coal_obs,times,t_correct,N,gridsize=1000, niter = 1000, burn = 0, thin = 5,changetime, DEMS=c("S","I"),
                         prior=list(pop_pr=c(1,1,1,10), R0_pr=c(1,7), mu_pr = c(3,0.2), gamma_pr = c(3,0.2), hyper_pr = c(0.001,0.001)),
                         proposal = list(pop_prop = 0.5, R0_prop = c(0.01), mu_prop=0.1, gamma_prop = 0.2, hyper_prop=0.05),
                         control = list(), updateVec = c(1,1,1,1,1,1,1), likelihood = "volz",model = "SIR",
                         Index = c(0,2), nparam=2, method = "seq",options = list(joint = F, PCOV = NULL,beta = 0.05, burn1 = 5000, parIdlist = NULL, priorIdlist = NULL), verbose = T){

  MCMC_setting = MCMC_setup_general(coal_obs, times,t_correct,N,gridsize,niter,burn,
                                    thin,changetime, DEMS,prior,proposal,
                                    control,likelihood,model,Index,nparam, options$PCOV)
  p = MCMC_setting$p
  nparam = sum(MCMC_setting$x_i[1:2]) + 1

  MCMC_obj = MCMC_initialize_general(MCMC_setting)

  if (MCMC_setting$likelihood == "volz") { # volz likelihood model

    params = matrix(nrow = niter, ncol = nparam + MCMC_obj$p)

  } else if (MCMC_setting$likelihood == "structural") { # structured coalescent model

    params = matrix(nrow = niter, ncol =  nparam + MCMC_obj$p)

  }else{ # other models
    params = matrix(nrow = niter, ncol = sum(MCMC_setting$x_i) + MCMC_obj$p)
  }

  l = numeric(niter)
  l1 = l
  l2 = l
  tjs = array(dim = c(dim(MCMC_obj$LatentTraj),niter))



  #' updateVec
  #' 1 parameter for intial state
  #' 2 R0
  #' 3 gamma
  #' 4 mu
  #' 5 hyper-parameter controls the smoothness of changepoints
  #' 6 changepoints
  #' 7 LNA noise in trajectory

  logIndexAll = c(1,0,1,0)
  parIndexALL = c(p:(p + MCMC_setting$x_i[2]), nparam+p)

  parId = parIndexALL[updateVec[c(1:3,5)] == 1]
  logId = logIndexAll[updateVec[c(1:3,5)] == 1]
  priorId = which(updateVec[c(1:4,5)] == 1)
  proposeId = which(updateVec[c(1:4,5)] == 1)

  print(paste("parameter ID", parId))
  print(paste("logId", logId))
  print(paste("priorId", priorId))
  print(paste("proposeId", proposeId))
  print("Warming up")
  for (i in 1:MCMC_setting$niter) {
    if(i %% 10000 == 0){
      print(i)
    }
    if (i %% 100 == 0 && verbose == T) {
      print(i)
      print(MCMC_obj$par)
      print(l1[i-1])
      print(l2[i-1])

      plot(MCMC_obj$LatentTraj[,1], MCMC_obj$LatentTraj[, MCMC_setting$x_i[4] + 2],type="l",xlab = "time", ylab = "Infected")
      lines(MCMC_obj$Ode_Traj_coarse[,1],MCMC_obj$Ode_Traj_coarse[, MCMC_setting$x_i[4] + 2],col="red",lty=2)
    }

   if(i == options$burn1){
      print("end warm up phase")
      idx = floor(options$burn1/2) : (options$burn1-1)
      if(!is.null(MCMC_setting$PCOV)){
        SigmaN = NULL
        for(i in 1:length(parId)){
          if(logId[i] == 1){
          SigmaN = cbind(SigmaN, log(params[idx,parId[i]]))
          }else{
            SigmaN = cbind(SigmaN, params[idx, parId[i]])
          }
        }
        SigmaN = (2.38^2) * cov(SigmaN) / length(parId)
        MCMC_setting$PCOV = options$beta * SigmaN + (1 - options$beta) * diag(rep(1,length(parId))) * 0.01 / length(parId)
      }
      print(MCMC_setting$PCOV)
    }

    if(i < options$burn1){

     # if(method == "seq"){
        MCMC_obj = Update_Param_each(MCMC_obj, MCMC_setting, options$parIdlist[-4],
                                     options$isLoglist[-4], options$priorIdlist[-4],options$priorIdlist[-4])$MCMC_obj
      #}else{
      #  MCMC_obj = update_Param_joint(MCMC_obj, MCMC_setting, method = "jointProp", parId, logId, priorId, proposeId)$MCMC_obj
      #}
      if(updateVec[6] == 1){
        #MCMC_obj = update_ChangePoint_general_NC(MCMC_obj,MCMC_setting,i)$MCMC_obj
        MCMC_obj = tryCatch({update_ChangePoint_ESlice(MCMC_obj,MCMC_setting,i)},
                            error = function(cond){
                              message(cond)
                              # Choose a return value in case of error
                              return(MCMC_obj)
                            })
      }
    }else{

      #MCMC_obj = update_Param_joint(MCMC_obj, MCMC_setting, method = "jointProp", parId = unlist(options$parIdlist),
       #                             unlist(options$isLoglist), unlist(options$priorIdlist), unlist(options$priorIdlist))$MCMC_obj
      if(method == "seq"){
        MCMC_obj = Update_Param_each(MCMC_obj, MCMC_setting, options$parIdlist,
                                     options$isLoglist, options$priorIdlist,options$priorIdlist)$MCMC_obj
      }else{
        MCMC_obj = update_Param_joint(MCMC_obj, MCMC_setting, method, parId, logId, priorId, proposeId)$MCMC_obj
      }



      if(updateVec[6] == 1){
        #MCMC_obj = update_ChangePoint_general_NC(MCMC_obj,MCMC_setting,i)$MCMC_obj
        MCMC_obj = tryCatch({update_ChangePoint_ESlice(MCMC_obj,MCMC_setting,i)},
                            error = function(cond){
                              message(cond)
                              # Choose a return value in case of error
                              return(MCMC_obj)
                            })
      }
      if(updateVec[7] == 1){
        MCMC_obj = tryCatch({updateTraj_general_NC(MCMC_obj,MCMC_setting,i)$MCMC_obj},
                            error = function(cond){
                              message(cond)
                              # Choose a return value in case of error
                              return(MCMC_obj)
                            })
      }
    }

    tjs[,,i] = MCMC_obj$LatentTraj
    params[i,] = MCMC_obj$par
    l[i] =  MCMC_obj$priorLog[MCMC_obj$p + MCMC_setting$x_i[4] + 1] #+ MCMC_obj$LogAlpha1 + MCMC_obj$LogAlpha2
    l1[i] = MCMC_obj$coalLog
    l2[i] = sum(MCMC_obj$priorLog) + MCMC_obj$coalLog + MCMC_obj$logOrigin
  }
  return(list(par = params,Trajectory = tjs,l=l,l1=l1,l2 = l2))
}

