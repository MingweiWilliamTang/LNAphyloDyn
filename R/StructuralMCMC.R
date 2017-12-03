# General structural coalescent model

#' changetime: the time point for changing infection rate
#'
#'
#'
#'
#'


MCMC_setup_structural = function(bdt,times,t_correct,N,gridsize=50,niter = 1000,burn = 500,thin = 5, changetime,DEMS = c("E", "I"),
                                 prior=list(pop_pr=c(1,10), R0_pr=c(1,7), mu_pr = c(3,0.2), gamma_pr = c(3,0.2), ch_pr = 1),
                                 proposal = list(pop_prop = 1, R0_prop = c(0.01), mu_prop=0.1, gamma_prop = 0.2, ch_prop=0.05),
                                 control = list(), likelihood = "structural",model = "SEIR2",Index = c(0,2),nparam=3, PCOV = NULL){


  gridset = seq(1,length(times),by=gridsize)
  grid = times[gridset]
  datedphylo = as.DatePhylo(bdt,endTime = max(times), matrix(rep(c(0,1),each=bdt$Nnode + 1),ncol=2))
  #if(max(time) < bdt$maxHeight){
   # stop("Tree height is larger than the interval length")
  #}

  Init_Detail = colikcpp_prepare(datedphylo, grid, DEMS)
  if(is.null(t_correct)){
    t_correct = bdt$maxHeight
  }

  x_i = c(length(changetime),nparam,Index[1],Index[2])
  MCMC_setting = list(Init_Detail = Init_Detail, times = times,t_correct = t_correct,x_r = c(N,changetime),
                      gridsize = gridsize,gridset = gridset, niter = niter,burn = burn,thin = thin,x_i = x_i,
                      prior = prior, proposal = proposal, control = control, p = length(DEMS),
                      reps = 1, likelihood = likelihood, model = model, PCOV = PCOV)
  cat("MCMC set up ready \n")

  return(MCMC_setting)
}


MCMC_initialize_stuctural = function(MCMC_setting){

  logMultiNorm = NaN
  coalLog = NaN

  ########
  priorLog = numeric(MCMC_setting$p+4)
  while(is.nan(logMultiNorm)||is.nan(coalLog)){

    # initialize the parameters
    alpha = numeric(MCMC_setting$p)
    if(is.null(MCMC_setting$control$alpha)){
      alpha = numeric(MCMC_setting$p)
      for(i in 1:MCMC_setting$p){
        alpha[i] = exp(rnorm(1,MCMC_setting$prior$pop_pr[2 * i - 1],MCMC_setting$prior$pop_pr[2 * i]))
      }
    }else{
      for(i in 1:MCMC_setting$p){
        alpha[i] = MCMC_setting$control$alpha[i]
      }
    }
    state = numeric(MCMC_setting$p)
    for(i in 1:MCMC_setting$p){
      state[i] = MCMC_setting$x_r[1] * alpha[i]
      priorLog[i] = dnorm(log(alpha[i]), MCMC_setting$prior$pop_pr[2 * i - 1],MCMC_setting$prior$pop_pr[2 * i],log = T)
    }

    if(is.null(MCMC_setting$control$R0)){
      R0 = runif(1,MCMC_setting$prior$R0_pr[1],MCMC_setting$prior$R0_pr[2])
    }else{
      R0 = MCMC_setting$control$R0
    }


    if(is.null(MCMC_setting$control$mu)){
      mu = exp(rnorm(1,MCMC_setting$prior$mu_pr[1],MCMC_setting$prior$mu_pr[2]))
    }else{
      mu = MCMC_setting$control$mu
    }
    priorLog[MCMC_setting$p+2] = dnorm(log(mu),MCMC_setting$prior$mu_pr[1],MCMC_setting$prior$mu_pr[2],log = T)

    if(is.null(MCMC_setting$control$gamma)){
      gamma = exp(rnorm(1,MCMC_setting$prior$gamma_pr[1], MCMC_setting$prior$gamma_pr[2]))
    }else{
      gamma = MCMC_setting$control$gamma
    }
    priorLog[MCMC_setting$p+3] = dnorm(log(gamma),MCMC_setting$prior$gamma_pr[1],MCMC_setting$prior$gamma_pr[2],log = T)
    ch=c()
    if(length(MCMC_setting$x_r) > 1){
      if(is.null(MCMC_setting$control$ch)){
        ch = rlnorm(MCMC_setting$x_i[1],0,MCMC_setting$chpr)
      }else{
        ch = MCMC_setting$control$ch
      }
    }
    if(MCMC_setting$model == "SEIR" || MCMC_setting$model == "SEIR2"){
      param = c(R0,mu,gamma,ch)
    }else{
      param = c(R0, gamma, ch)
    }


    Ode_Traj_thin = ODE_rk45(state, MCMC_setting$times, param, MCMC_setting$x_r, MCMC_setting$x_i,model = MCMC_setting$model)

    Ode_Traj_coarse = Ode_Traj_thin[MCMC_setting$gridset,]

    FT = KF_param_chol(Ode_Traj_thin,param,MCMC_setting$gridsize, MCMC_setting$x_r,MCMC_setting$x_i,model = MCMC_setting$model)


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

    coalLog = Structural_Coal_lik(MCMC_setting$Init_Detail,LatentTraj,param, MCMC_setting$x_r, MCMC_setting$x_i,
                                  model = "SEIR2")
    print(paste("coalescent likelihood after initialization ", coalLog))

    if(!is.nan((coalLog))){
      coalLog = ifelse(coalLog <= - 10000000,NaN, coalLog)
    }
    plot(LatentTraj[,1],LatentTraj[,2],type="l",xlab = "time",col="blue", lwd = 2)
    lines(LatentTraj[,1],LatentTraj[,3],col="red", lwd = 2)
  }


  chprob = 0

  for(i in 1:MCMC_setting$x_i[1]){
    chprob = chprob + dnorm(log(param[MCMC_setting$x_i[2] + i]), 0, MCMC_setting$prior$ch_pr,log = T)
  }
  priorLog[MCMC_setting$p + MCMC_setting$x_i[2] +1] = chprob
  #LogGamma = dlnorm(gamma,MCMC_setting$prior$gamma_pr[1],MCMC_setting$prior$gamma_pr[2],log = T)
  #LogE = dlnorm(mu, MCMC_setting$prior$mu_pr[1], MCMC_setting$prior$mu_pr[2],log = T)
  plot(LatentTraj[,1],LatentTraj[,2],type="l",xlab = "time",col="blue", lwd = 2)
  lines(LatentTraj[,1],LatentTraj[,3],col="red", lwd = 2)

  paras =  c(state,param)
  MCMC_obj = list(par = paras,LatentTraj = LatentTraj, logMultiNorm = logMultiNorm,p = MCMC_setting$p,
                  Ode_Traj_coarse = Ode_Traj_coarse, FT = FT, coalLog = coalLog, OriginTraj = OriginTraj,logOrigin = logOrigin,
                  priorLog = priorLog)
  ##########
  cat("Initialize MCMC \n")
  print(paste("population size = ", MCMC_setting$N))
  print(paste(paras))
  return(MCMC_obj)
}


####################


################







SEIR_general_MCMC_Structural = function(bdt,times,t_correct,N,gridsize=1000, niter = 1000, burn = 0, thin = 5,
                                        changetime,DEMS = c("E", "I"),
                                        prior=list(pop_pr=c(1,10,1,10), R0_pr=c(1,7), mu_pr = c(3,0.2), gamma_pr = c(3,0.2), ch_pr = 1),
                                        proposal = list(pop_prop = 1, R0_prop = c(0.01), mu_prop=0.1, gamma_prop = 0.2, ch_prop=0.05),
                                        control = list(), updateVec = c(1,1,1,1,1,1), likelihood = "structural", model = "SEIR2",
                                        Index = c(0,2), nparam=3, joint = F, PCOV = NULL,beta = 0.05, burn1 = 5000){

  MCMC_setting = MCMC_setup_structural(bdt, times,t_correct,N,gridsize,niter,burn,
                                       thin,changetime, DEMS,prior,proposal,
                                       control,likelihood,model,Index,nparam,PCOV)
  nparams = length(changetime) + nparam
  MCMC_obj = MCMC_initialize_stuctural(MCMC_setting)
  if(MCMC_setting$likelihood == "volz"){
    params = matrix(nrow = niter, ncol =  sum(MCMC_setting$x_i[1:2]) + MCMC_obj$p)
    ARMS = matrix(nrow = niter, ncol =  sum(MCMC_setting$x_i[1:2]) + MCMC_obj$p)
  }else if(MCMC_setting$likelihood == "structural"){
    params = matrix(nrow = niter, ncol =  sum(MCMC_setting$x_i[1:2]) + MCMC_obj$p)
    ARMS = matrix(nrow = niter, ncol =  sum(MCMC_setting$x_i[1:2]) + MCMC_obj$p)
  }else{
    params = matrix(nrow = niter, ncol = sum(MCMC_setting$x_i) + MCMC_obj$p)
    ARMS = matrix(nrow = niter, ncol = sum(MCMC_setting$x_i) + MCMC_obj$p)
  }

  l = numeric(niter)
  l0 = l
  l1 = l
  l2 = l
  tjs = NULL

  for (i in 1:MCMC_setting$niter) {

    if (i %% 100 == 0) {
      print(i)
      print(MCMC_obj$par)
      print(l2[i-1])
      print(l0[i-1])
      plot(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2],type="l",col = "blue",lwd=2)
      lines(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,3],col = "red",lwd=2)
      lines(MCMC_obj$Ode_Traj_coarse[,1], MCMC_obj$Ode_Traj_coarse[,2], col="blue",lty=2)
      lines(MCMC_obj$Ode_Traj_coarse[,1], MCMC_obj$Ode_Traj_coarse[,3], col="red",lty=2)
    }
    ARvec = numeric(dim(ARMS)[2])
    if(i == burn1){
      idx = floor(burn1/2) : (burn1-1)
      SigmaN = (2.38^2) * cov(cbind(params[idx,MCMC_obj$p + MCMC_setting$x_i[3] + 1],
                                    log(params[idx,-(1:(MCMC_obj$p+1))]))) / nparams
      MCMC_setting$PCOV = beta * SigmaN + (1 - beta) * diag(rep(1,nparams)) * 0.01 / nparams
      print(MCMC_setting$PCOV)
    }

    if(i < burn1 && joint == F){
      if(updateVec[1] == 1){
        step1 = updateAlphas_General(MCMC_obj,MCMC_setting,i)
        # print(c(MCMC_obj$coalLog,MCMC_obj$logMultiNorm))
        MCMC_obj = step1$MCMC_obj
        ARvec[1] = step1$AR
      }
      if(updateVec[2] == 1){
        #step2 = updateR0_general2(MCMC_obj,MCMC_setting,i)
        step2 = updateR0_general_NC(MCMC_obj, MCMC_setting, i)
        MCMC_obj = step2$MCMC_obj
        ARvec[2] = step2$AR
      }
      if(updateVec[3] == 1){
        step3 = updategamma_general_NC(MCMC_obj,MCMC_setting,i)
        MCMC_obj = step3$MCMC_obj
        ARvec[3] = step3$AR
      }

      if(updateVec[4] == 1){
        step3 = updateE_general_NC(MCMC_obj,MCMC_setting,i,idx = 2)
        MCMC_obj = step3$MCMC_obj
        ARvec[4] = step3$AR
      }
      if(updateVec[5] == 1){
        step4 = update_ChangePoint_general_NC(MCMC_obj,MCMC_setting,i)
        MCMC_obj = step4$MCMC_obj
      }
    }else if(i < burn1 && joint == T){
        stepN = update_All_general_NC(MCMC_obj, MCMC_setting, i)
        MCMC_obj = stepN$MCMC_obj
      }else{
      stepN = update_All_general_NC(MCMC_obj, MCMC_setting, i)
      MCMC_obj = stepN$MCMC_obj
      if(updateVec[6] == 1){
        MCMC_obj = updateTraj_general_NC(MCMC_obj, MCMC_setting, i)$MCMC_obj
      }
    }
    #step4 = updateLambda_SIRS(MCMC_obj,MCMC_setting,i)
    #ARvec[8] = step4$AR
    #MCMC_obj = step4$MCMC_obj
    ARMS[i,] = ARvec
    tjs = abind(tjs,MCMC_obj$LatentTraj,along = 3)
    params[i,] = MCMC_obj$par
    l[i] =  MCMC_obj$priorLog[3]
    l0[i] = MCMC_obj$logOrigin#+ MCMC_obj$LogAlpha1 + MCMC_obj$LogAlpha2
    l1[i] = MCMC_obj$coalLog
    l2[i] = MCMC_obj$coalLog + sum(MCMC_obj$priorLog) + l0[i]
  }
  return(list(par = params,Trajectory = tjs,l = l,l1 = l1,l2 = l2,l0 =l0,AR = ARMS))
}


