
SIR_LNA_ESlice_MCMC4 =function(Init,times,t_correct,N,gridsize=50,niter = 1000,burn = 500,thin = 5,
                               a1 = 10, a2 = 20, c1=-2.3,d1 = 200, d2 =40){


  # build a List to Store MCMC results
  #' theta1
  #' theta2
  #' Trajectory matrix: in log-scale

  #' state Initial state in 2D (X,Y)
  #' OdeTraj_thin: ode trajectory with high resoluton, used to generate the filter
  #' OdeTraj_coarse: ode trajectory with low resultion, used to generate the path
  #' LatentTraj: the latent SIR 3d trajectory, has the same resolution as OdeTraj_coarse
  #' logMultNorm: loglikelihood of Ode trajectory
  #'
  #' Infectious rate and recovery rate: theta1, theta2
  #' Scale parameter: lambda
  #'
  #' alpha1: log(alpha1) ~ gamma(60,a1)
  #' alpha2: log(alpha2) ~ gamma(60,a2)
  #' prior
  #' s1:  s1 = theta1*N/theta2  U[1,10]
  #' propose function random walk Uniform[-0.5,0.5] on log space
  #'
  #' s2: theta2 log(s2) ~ norm(c1,0.4)
  #'
  #'
  #' lambda log(lambda) ~ gamma(d1,d2) propose random walk uniform [-0.5,0.5] on logspace
  #'
  #'
  #' log(p1/p3) ~ gamma(60,10)
  #' log(p2/p3) ~ gamma(60,20)


  # fix the inital state
  # initialize mcmc
  q = 1

  alpha1 = exp(rgamma(1,60,a1))
  alpha2 = exp(rgamma(1,60,a2))
  S = N * alpha1 / (alpha1 + alpha2 + 1)
  I = N * alpha2 / (alpha1 + alpha2 + 1)
  state = c(X = S, Y = I)

  countInf = 0
  difs = NULL
  rec = matrix(rep(0,5*niter),ncol=5,nrow=niter)
  logMultiNorm = NaN
  ########
  s1 = runif(1,1,10)
  s2 = exp(rnorm(1,c1,0.4))
  theta1 = s1 * s2 / N
  theta2 = s2
  param = c(theta1 = theta1, theta2 = theta2)
  gridset=seq(1,length(times),by=gridsize)
  #Ode_Traj_thin <- ODE(log(state),times,param)

  Ode_Traj_thin = ODE(log(state), times,
                      param)

  Ode_Traj_coarse = Ode_Traj_thin[gridset,]


  FT = SIR_log_KOM_Filter2(Ode_Traj_thin,theta1,theta2,gridsize)


  Latent = Traj_sim(log(state),Ode_Traj_coarse,FT)
  LatentTraj = Latent$SimuTraj

  # lambda = exp(rgamma(1,16,4))

  logMultiNorm = Latent$loglike

  ##########
  while(is.nan(logMultiNorm)){
    s1 = runif(1,1,10)
    s2 = exp(rnorm(1,c1,0.4))
    theta1 = s1 * s2 / N
    theta2 = s2
    param = c(theta1 = theta1, theta2 = theta2)
    gridset=seq(1,length(times),by=gridsize)
    #Ode_Traj_thin <- ODE(log(state),times,param)

    Ode_Traj_thin = ODE(log(state), times,
                        param)

    Ode_Traj_coarse = Ode_Traj_thin[gridset,]


    FT = SIR_log_KOM_Filter2(Ode_Traj_thin,theta1,theta2,gridsize)


    Latent = Traj_sim(log(state),Ode_Traj_coarse,FT)
    LatentTraj = Latent$SimuTraj

    # lambda = exp(rgamma(1,16,4))

    logMultiNorm = Latent$loglike
  }
  lambda = 100
  coalLog = coal_loglik(Init,LatentTraj,t_correct,lambda,gridsize)
  MCMC_para = NULL
  MCMC_Traj = NULL
  # MCMC_para = matrix(nrow = niter,ncol = 2)
  cat("begin MCMC")
  L_joint = numeric(niter + 1)
  L_joint[1] = coalLog + logMultiNorm
  #+ dgamma(-log(s1),50,5,log = T) +
  #  dnorm(s2,2500,500,log = T)
  #dgamma(log(alpha1),60,10,log = T) +
  #dgamma(log(alpha2),60,20,log = T)

  for(i in 1:niter){

    ############
    # update the state
    alpha1_new = alpha1 * exp(runif(1,-0.2,0.2))
    alpha2_new = alpha2 * exp(runif(1,-0.2,0.2))
    state_new = c(X = N * alpha1_new / (alpha1_new + alpha2_new + 1),
                  Y = N * alpha2_new / (alpha1_new + alpha2_new + 1))
    Ode_Traj_thin_new <- ODE(log(state_new),times,
                             param)

    Ode_Traj_coarse_new = Ode_Traj_thin_new[gridset,]

    FT_new = SIR_log_KOM_Filter2(Ode_Traj_thin_new,theta1,theta2,gridsize)

    if(i<burn){
      print(1)
      LatentTraj_new = ESlice(LatentTraj,Ode_Traj_coarse_new,FT_new,log(state_new),
                              Init,t_correct,lambda,reps = 1,gridsize)

      logMultiNorm_new = log_like_traj(LatentTraj_new,Ode_Traj_coarse_new,FT_new,gridsize)

      coalLog_new = coal_loglik(Init,LatentTraj_new,t_correct,lambda,gridsize)
    }else{
      logMultiNorm_new = log_like_traj2(LatentTraj,times,log(state_new),theta1,theta2,gridsize)
      coalLog_new = coalLog
    }
    if(is.nan(logMultiNorm_new)){
      logMultiNorm_new = -Inf
      countInf = countInf + 1
    }

    a = min(exp(dgamma(log(alpha1_new),60,a1,log = T) + dgamma(log(alpha2_new),60,a2,log = T) +
                  coalLog_new + (logMultiNorm_new - logMultiNorm)/q -
                  ( coalLog+ dgamma(log(alpha1),60,a1,log = T) +
                      dgamma(log(alpha2),60,a2,log = T))
    ),1)
    if(is.na(a)){a = -1}
    if(runif(1,0,1) < a){
      rec[i,1:2] = c(1,1)
      state=state_new
      alpha1 = alpha1_new
      alpha2 = alpha2_new
      Ode_Traj_coarse = Ode_Traj_coarse_new
      logMultiNorm = logMultiNorm_new
      FT = FT_new
      if(i<burn){
        coalLog = coalLog_new
        LatentTraj = LatentTraj_new
      }
    }

    LatentTraj= ESlice(LatentTraj,Ode_Traj_coarse,FT,log(state),
                       Init,t_correct,lambda,reps = 1,gridsize)

    logMultiNorm = log_like_traj2(LatentTraj,times,log(state),theta1,theta2,gridsize)


    ######################




    # Sample s1
    # propose theta1_new condition on theta1

    s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 1),10)
    theta1_new = s1_new * s2 / N
    param_new = c(theta1 = theta1_new, theta2 = theta2)
    #   print(param_new)
    #    Ode_Traj_thin_new <- ode(y = log(state), times = times,
    #                        func = SIR.log.ode2, parms = param_new)
    Ode_Traj_thin_new <- ODE(log(state),times,
                             param_new)

    Ode_Traj_coarse_new = Ode_Traj_thin_new[gridset,]

    FT_new = SIR_log_KOM_Filter2(Ode_Traj_thin_new,theta1_new,theta2,gridsize)
    if(i < burn){

      print(2)
      LatentTraj_new = ESlice(LatentTraj,Ode_Traj_coarse_new,FT_new,log(state),
                              Init,t_correct,lambda,reps = 1,gridsize)

      logMultiNorm_new = log_like_traj(LatentTraj_new,Ode_Traj_coarse_new,FT_new,gridsize)

      coalLog_new = coal_loglik(Init,LatentTraj_new,t_correct,lambda,gridsize)
    }else{
      #ttt = Traj_sim(log(state),Ode_Traj_coarse_new,FT_new)
      #LatentTraj_new = ttt$SimuTraj
      #logMultiNorm_new = ttt$loglike
      #coalLog_new = coal_loglik(Init,LatentTraj_new,t_correct,lambda)
      logMultiNorm_new = log_like_traj2(LatentTraj,times,log(state),theta1_new,theta2,gridsize)
      coalLog_new = coalLog
    }
    if(is.nan(logMultiNorm_new)){
      countInf = countInf + 1
      a = -1
    }else{
      a = min(exp(coalLog_new + (logMultiNorm_new - logMultiNorm)/q -  coalLog
      ),1)
    }


    if(runif(1,0,1) < a){
      rec[i,3] = 1
      theta1 = theta1_new
      Ode_Traj_coarse = Ode_Traj_coarse_new
      logMultiNorm = logMultiNorm_new
      param = param_new
      s1=s1_new
      FT = FT_new

      if(i<burn){
        coalLog = coalLog_new
        LatentTraj = LatentTraj_new
      }
    }

    LatentTraj= ESlice(LatentTraj,Ode_Traj_coarse,FT,log(state),
                       Init,t_correct,lambda,reps = 1,gridsize)
    logMultiNorm = log_like_traj2(LatentTraj,times,log(state),theta1,theta2,gridsize)


    # Sample theta2
    # propose theta_2_new condition on theta2

    s2_new = s2 * exp(runif(1,-0.3,0.3))
    theta1_new = s1 * s2_new / N
    theta2_new = s2_new

    param_new = c(theta1 = theta1_new, theta2 = theta2_new)


    Ode_Traj_thin_new <- ODE(log(state), times,
                             param_new)

    Ode_Traj_coarse_new = Ode_Traj_thin_new[gridset,]
    FT_new = SIR_log_KOM_Filter2(Ode_Traj_thin_new,theta1_new,theta2_new,gridsize)
    if(i < burn){

      LatentTraj_new = ESlice(LatentTraj,Ode_Traj_coarse_new,FT_new,log(state),
                              Init,t_correct,lambda,reps = 1,gridsize)

      logMultiNorm_new = log_like_traj(LatentTraj_new,Ode_Traj_coarse_new,FT_new,gridsize)

      coalLog_new = coal_loglik(Init,LatentTraj_new,t_correct,lambda,gridsize)
    }else{

      logMultiNorm_new = log_like_traj2(LatentTraj,times,log(state),theta1_new,theta2_new,gridsize)
      coalLog_new = coalLog

    }
    if(is.nan(logMultiNorm_new)){
      countInf = countInf + 1
      a = -1
    }else{
      a = min(c(exp(dnorm(log(s2_new),c1,0.4,log = T) +coalLog_new +(logMultiNorm_new - logMultiNorm)/q -
                      ( coalLog + dnorm(log(s2),c1,0.4,log = T))),1))
    }


    if(runif(1,0,1) < a){
      rec[i,4] = 1
      theta1 = theta1_new
      theta2 = theta2_new
      Ode_Traj_coarse = Ode_Traj_coarse_new
      logMultiNorm = logMultiNorm_new
      param = param_new
      s2 = s2_new
      FT = FT_new
      if(i<burn){
        coalLog = coalLog_new
        LatentTraj = LatentTraj_new
      }
    }
    LatentTraj= ESlice(LatentTraj,Ode_Traj_coarse,FT,log(state),
                       Init,t_correct,lambda,reps = 1,gridsize)
    logMultiNorm = log_like_traj2(LatentTraj,times,log(state),theta1,theta2,gridsize)



    # Sample the trajectory using Elipitical Slice Sampling
    # LatentTraj = ESlice(LatentTraj,Ode_Traj_coarse,FT,log(state),Init,t_correct,lambda,reps = 10)

    #logMultiNorm = log_like_traj2(LatentTraj,times,log(state),theta1,theta2,gridsize)
    # Sample the scale parameter lambda

    lambda_new = lambda * exp(runif(1,-0.3,0.3))
    coalLog_new = coal_loglik(Init,LatentTraj,t_correct,lambda_new,gridsize)
    a = min(c(exp(coalLog_new - coalLog + dgamma(log(lambda_new),d1,d2,log = T) -
                    dgamma(log(lambda),d1,d2,log=T)), 1))
    if(runif(1,0,1) < a){
      rec[i,5] = 1
      coalLog = coalLog_new
      lambda = lambda_new
    }
    # print(lambda)

    # store the MCMC samples

    L_joint[i + 1] = coalLog + dnorm(log(s2),c1,0.4,log = T) +  logMultiNorm + dgamma(log(alpha1),60,a1,log = T) +
      dgamma(log(alpha2),60,a2,log = T) + dgamma(log(lambda),d1,d2,log=T)
    #+ dgamma(-log(s1),60,5,log = T) +
    #  dnorm(s2,2500,500,log = T)
    #    dgamma(log(alpha1),60,10,log = T) +
    #   dgamma(log(alpha2),60,20,log = T)
    if(i %% 100 == 0) {
      cat(floor(i)," iterations complete","\n")
    }
    if( i %% thin ==1 && i >=burn){
      MCMC_para = cbind(MCMC_para,c(state,theta1,theta2,lambda))
      MCMC_Traj = abind(MCMC_Traj,LatentTraj,along = 3)
    }
  }
  print(paste("num of Inf",countInf))
  return(list(par=MCMC_para,Traj=MCMC_Traj,OdeTraj = Ode_Traj_coarse,
              FT=FT,L=L_joint,
              records = rec,dif=difs,priorinfo = c(a1,a2,c1,d1,d2)))
}


posterior_hist = function(MCMC_para,t_para = NULL){
  dev.off()
  par(mfrow = c(3,2))
  hist(MCMC_para[1,], main = "S")
  abline(v = t_para[1],col = "red", lwd=2)
  hist(MCMC_para[2,], main = "I")
  abline(v = t_para[2],col = "red", lwd=2)
  hist(MCMC_para[3,], main = "infection rate beta")
  abline(v = t_para[3],col = "red", lwd=2)
  hist(MCMC_para[4,], main = "recover rate")
  abline(v = t_para[4],col = "red", lwd=2)
  hist(MCMC_para[5,], main = "lambda")
  abline(v = t_para[5],col = "red", lwd=2)
  plot(0,type='n',axes=FALSE,ann=FALSE)
}

prior_dist = function(N,pr_par){
  alpha1 = rgamma(1000,60,pr_par[1])
  alpha2 = rgamma(1000,60,pr_par[2])
  s1 = runif(1000,1,10)
  s2 = exp(rnorm(1000,pr_par[3],0.4))
  lambda = exp(rgamma(1000,pr_par[1],pr_par[2]))
  S = N * alpha1/(alpha1 + alpha2 + 1)
  I = N * alpha2/(alpha1 + alpha2 + 1)
  theta1 =  s1 * s2 / N
  theta2 = s2
  par(mfrow=c(3,2))
  hist(S,main = "S")
  hist(I,main = "I")
  hist(theta1,main = "infection rate beta")
  hist(theta2,main = "recover rate gamma")
  hist(lambda,main = "scale parameter")
  plot(0,type='n',axes=FALSE,ann=FALSE)
}


posterior_par_line = function(MCMC_para,index){
  dev.off()
  par(mfrow = c(3,2))
  plot(MCMC_para[1,index], type = "l", main = "S")
  abline(v = t_para[1],col = "red", lwd=2)
  plot(MCMC_para[2,], type = "l", main = "I")
  abline(v = t_para[2],col = "red", lwd=2)
  plot(MCMC_para[3,], type= "l", main = "infection rate beta")
  abline(v = t_para[3],col = "red", lwd=2)
  plot(MCMC_para[4,], type= l,main = "recover rate")
  abline(v = t_para[4],col = "red", lwd=2)
  plot(MCMC_para[5,], type = "l", main = "lambda")
  abline(v = t_para[5],col = "red", lwd=2)
  plot(0,type='n',axes=FALSE,ann=FALSE)
}




updateAlphas = function(MCMC_obj,MCMC_setting,i){
  #alpha1 = MCMC_obj$par[1] / (MCMC_setting$N -  MCMC_obj$par[1] -  MCMC_obj$par[2])
  #alpha2 = MCMC_obj$par[2] / (MCMC_setting$N -  MCMC_obj$par[1] -  MCMC_obj$par[2])
  alpha1 = MCMC_obj$par[1] / MCMC_obj$par[2]
  alpha1_new = alpha1 * exp(runif(1,-0.2,0.2))

  #alpha2_new = alpha2 * exp(runif(1,-0.2,0.2))
  #state_new = c(X = MCMC_setting$N * alpha1_new / (alpha1_new + alpha2_new + 1),
  #             Y = MCMC_setting$N * alpha2_new / (alpha1_new + alpha2_new + 1))
  state_new = c(X = MCMC_setting$N * alpha1_new / (alpha1_new + 1),
                Y = MCMC_setting$N  / (alpha1_new + 1))

  Ode_Traj_thin_new <- ODE(log(state_new),MCMC_setting$times,
                           MCMC_obj$par[3:4])

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = SIR_log_KOM_Filter2(Ode_Traj_thin_new,MCMC_obj$par[3],MCMC_obj$par[4],MCMC_setting$gridsize)

  if(i < MCMC_setting$burn){

    LatentTraj_new = ESlice(MCMC_obj$LatentTraj,Ode_Traj_coarse_new,FT_new,log(state_new),
                            MCMC_setting$Init,MCMC_setting$t_correct,MCMC_obj$par[5],reps = MCMC_setting$reps,MCMC_setting$gridsize)

    logMultiNorm_new = log_like_traj(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

    coalLog_new = coal_loglik(MCMC_setting$Init,LatentTraj_new,MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)
  }else{

    LatentTraj_new =cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                            Ode_Traj_coarse_new[,2:3])
    logMultiNorm_new = log_like_traj(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)
    coalLog_new = coal_loglik(MCMC_setting$Init,LatentTraj_new,MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)

  }
  if(is.nan(logMultiNorm_new)){
    logMultiNorm_new = -Inf
    # countInf = countInf + 1
  }
  a = min(c(exp(dnorm(log(alpha1_new),MCMC_setting$b1,MCMC_setting$a1,log = T) + coalLog_new + #dgamma(log(alpha2_new),MCMC_setting$b2,MCMC_setting$a2,log = T) +
                  logMultiNorm_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog -
                  ( MCMC_obj$LogAlpha1)
  ),1))

  #print(logMultiNorm_new-log_like_traj2(MCMC_obj$LatentTraj,MCMC_setting$times,log(state_new),MCMC_obj$par[3],MCMC_obj$par[4],MCMC_setting$gridsize,MCMC_setting$t_correct ))

  #print(logMultiNorm_new - MCMC_obj$logMultiNorm)
  if(is.na(a)){a = -1}
  # print(c(logMultiNorm_new,MCMC_obj$logMultiNorm,dgamma(log(alpha1_new),60,MCMC_setting$a1,log = T), dgamma(log(alpha2_new),60,MCMC_setting$a2,log = T)))
  AR=0
  if(runif(1,0,1) < a){
    AR = 1
    state=state_new
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
  # MCMC_obj$LatentTraj = ESlice(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,log(MCMC_obj$par[1:2]),
  #                              MCMC_setting$Init,MCMC_setting$t_correct,MCMC_obj$par[5],reps = MCMC_setting$reps,MCMC_setting$gridsize)
  # q =  MCMC_obj$logMultiNorm
  #  MCMC_obj$logMultiNorm = log_like_traj(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
  #print(MCMC_obj$logMultiNorm - q)
  #  MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,MCMC_obj$LatentTraj,MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)

  return(list(MCMC_obj = MCMC_obj, AR = AR))
}



updateS1 = function(MCMC_obj, MCMC_setting, i){
  s1 = MCMC_obj$par[3] / MCMC_obj$par[4] * MCMC_setting$N
  s1_new = pmin(pmax(s1 + runif(1,-0.5,0.5), 1),10)
  theta1_new = s1_new * MCMC_obj$par[4] / MCMC_setting$N
  param_new = c(theta1 = theta1_new, theta2 = MCMC_obj$par[4])
  #   print(param_new)
  #    Ode_Traj_thin_new <- ode(y = log(state), times = times,
  #                        func = SIR.log.ode2, parms = param_new)
  Ode_Traj_thin_new <- ODE(log(MCMC_obj$par[1:2]),MCMC_setting$times,
                           param_new)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]

  FT_new = SIR_log_KOM_Filter2(Ode_Traj_thin_new,theta1_new,MCMC_obj$par[4],MCMC_setting$gridsize)

  if(i < MCMC_setting$burn){


    LatentTraj_new = ESlice(MCMC_obj$LatentTraj,Ode_Traj_coarse_new,FT_new,log(MCMC_obj$par[1:2]),
                            MCMC_setting$Init,MCMC_setting$t_correct,MCMC_obj$par[5],reps = MCMC_setting$reps,MCMC_setting$gridsize)

    logMultiNorm_new = log_like_traj(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

    coalLog_new = coal_loglik(MCMC_setting$Init,LatentTraj_new,MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)
  }else{
    #ttt = Traj_sim(log(state),Ode_Traj_coarse_new,FT_new)
    #LatentTraj_new = ttt$SimuTraj
    #logMultiNorm_new = ttt$loglike
    #coalLog_new = coal_loglik(Init,LatentTraj_new,t_correct,lambda)
    LatentTraj_new =cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                            Ode_Traj_coarse_new[,2:3])
    logMultiNorm_new = log_like_traj(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)
    coalLog_new = coal_loglik(MCMC_setting$Init,LatentTraj_new,MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)
  }

  if(is.nan(logMultiNorm_new)){
    a = -1
  }else{
    a = min(c(exp((logMultiNorm_new - MCMC_obj$logMultiNorm + coalLog_new - MCMC_obj$coalLog)
    ),1))
  }
  AR = 0
  if(runif(1,0,1) < a){
    AR = 1
    MCMC_obj$par[3] = theta1_new
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
  # MCMC_obj$logMultiNorm = log_like_traj(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
  #print(MCMC_obj$logMultiNorm - q)
  # MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,MCMC_obj$LatentTraj,MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)

  return(list(MCMC_obj = MCMC_obj, AR = AR))

}



updateS2 = function(MCMC_obj, MCMC_setting, i){
  s2_new = MCMC_obj$par[4] * exp(rnorm(1,0,0.1))
  theta1_new = MCMC_obj$par[3] / MCMC_obj$par[4] * s2_new
  theta2_new = s2_new

  param_new = c(theta1 = theta1_new, theta2 = theta2_new)


  Ode_Traj_thin_new <- ODE(log(MCMC_obj$par[1:2]), MCMC_setting$times,
                           param_new)

  Ode_Traj_coarse_new = Ode_Traj_thin_new[MCMC_setting$gridset,]
  FT_new = SIR_log_KOM_Filter2(Ode_Traj_thin_new,theta1_new,theta2_new,MCMC_setting$gridsize)

  if(i < MCMC_setting$burn){

    LatentTraj_new = ESlice(MCMC_obj$LatentTraj,Ode_Traj_coarse_new,FT_new,log(MCMC_obj$par[1:2]),
                            MCMC_setting$Init,MCMC_setting$t_correct,MCMC_obj$par[5],reps = MCMC_setting$reps,MCMC_setting$gridsize)

    logMultiNorm_new = log_like_traj(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)

    coalLog_new = coal_loglik(MCMC_setting$Init,LatentTraj_new,MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)
  }else{

    LatentTraj_new =cbind(MCMC_obj$LatentTraj[,1],MCMC_obj$LatentTraj[,2:3] - MCMC_obj$Ode_Traj_coarse[,2:3] +
                            Ode_Traj_coarse_new[,2:3])
    logMultiNorm_new = log_like_traj(LatentTraj_new,Ode_Traj_coarse_new,FT_new,MCMC_setting$gridsize,MCMC_setting$t_correct)
    coalLog_new = coal_loglik(MCMC_setting$Init,LatentTraj_new,MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)

  }
  if(is.nan(logMultiNorm_new)){
    a = -1
  }else{
    a = min(c(exp(dnorm(log(s2_new),MCMC_setting$c1,0.4,log = T) + logMultiNorm_new + coalLog_new - MCMC_obj$logMultiNorm - MCMC_obj$coalLog -
                    MCMC_obj$LogS2 ),1))
  }
  # print(theta2_new)
  AR = 0
  if(runif(1,0,1) < a){
    AR = 1
    MCMC_obj$par[3] = theta1_new
    MCMC_obj$par[4] = theta2_new
    MCMC_obj$Ode_Traj_coarse = Ode_Traj_coarse_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$LogS2 = dnorm(log(s2_new),MCMC_setting$c1,0.4,log = T)
    MCMC_obj$FT = FT_new
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$LatentTraj = LatentTraj_new

  }
  #MCMC_obj$LatentTraj = ESlice(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,log(MCMC_obj$par[1:2]),
  #                  MCMC_setting$Init,MCMC_setting$t_correct,MCMC_obj$par[5],reps = MCMC_setting$reps,MCMC_setting$gridsize)
  #q =  MCMC_obj$logMultiNorm
  #MCMC_obj$logMultiNorm = log_like_traj(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
  #print(MCMC_obj$logMultiNorm - q)

  #MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,MCMC_obj$LatentTraj,MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)

  return(list(MCMC_obj = MCMC_obj, AR = AR))
}


updateLambda = function(MCMC_obj,MCMC_setting, i){
  lambda_new = MCMC_obj$par[5] * exp(runif(1,-0.3,0.3))
  coalLog_new = coal_loglik(MCMC_setting$Init,MCMC_obj$LatentTraj,MCMC_setting$t_correct,lambda_new,MCMC_setting$gridsize)
  a = min(c(exp(coalLog_new - MCMC_obj$coalLog + dgamma(log(lambda_new),MCMC_setting$d1,MCMC_setting$d2,log = T) -
                  MCMC_obj$LogLambda), 1))
  AR = 0
  #print(coalLog_new - MCMC_obj$coalLog + dgamma(log(lambda_new),MCMC_setting$d1,MCMC_setting$d2,log = T) -
  #       MCMC_obj$LogLambda)
  if(runif(1,0,1) < a){
    # rec[i,5] = 1
    AR = 1
    MCMC_obj$LogLambda = dgamma(log(lambda_new),MCMC_setting$d1,MCMC_setting$d2,log = T)
    MCMC_obj$coalLog = coalLog_new
    MCMC_obj$par[5] = lambda_new
  }
  return(list(MCMC_obj = MCMC_obj, AR = AR))
}


#' likelihood free MCMC
#'
#'
#'

Like_Free = function(MCMC_obj,MCMC_setting,i){
  Sim_res = Traj_sim(log(MCMC_obj$par[1:2]),MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT, MCMC_setting$t_correct)
  LatentTraj_new = Sim_res$SimuTraj
  logMultiNorm_new = Sim_res$loglike
  coalLog_new = coal_loglik(MCMC_setting$Init, LatentTraj_new, MCMC_setting$t_correct, MCMC_obj$par[5], MCMC_setting$gridsize)
  a = min(c(exp(coalLog_new - MCMC_obj$coalLog),1))
  if(runif(1,0,1) < a){
    MCMC_obj$LatentTraj = LatentTraj_new
    MCMC_obj$logMultiNorm = logMultiNorm_new
    MCMC_obj$coalLog = coalLog_new
  }
  return(list(MCMC_obj = MCMC_obj))
}


post_like_eval = function(param,LatentTraj,MCMC_setting){
  alpha1 = param[1] / (MCMC_setting$N - param[1] -param[2])
  alpha2 = param[2] / (MCMC_setting$N - param[1] -param[2])
  s1 = param[3] * MCMC_setting$N / param[4]
  s2 = param[4]
  LogAlpha1 = dnorm(log(alpha1),MCMC_setting$b1,MCMC_setting$a1,log = T)
  #LogAlpha2 = dgamma(log(alpha2),MCMC_setting$b2,MCMC_setting$a2,log = T)
  LogS2 = dnorm(log(s2),MCMC_setting$c1,0.4)
  LogLambda = dnorm(log(param[5]),MCMC_setting$d1,MCMC_setting$d2,log = T)
  LogMultiNorm = log_like_traj2(LatentTraj,MCMC_setting$times,log(param[1:2]),param[1],param[2],MCMC_setting$gridsize,MCMC_setting$t_correct)
  coalLog = coal_loglik(MCMC_setting$Init, LatentTraj, MCMC_setting$t_correct, param[5],MCMC_setting$gridsize)
  l = LogAlpha1 + LogAlpha2 + LogS2 + LogLambda + coalLog + LogMultiNorm
  return(l)
}

updateTraj = function(MCMC_obj,MCMC_setting,i){
  #print(c(MCMC_obj$par,MCMC_obj$coalLog + MCMC_obj$logMultiNorm))
  MCMC_obj$LatentTraj = ESlice(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,log(MCMC_obj$par[1:2]),
                               MCMC_setting$Init,MCMC_setting$t_correct,MCMC_obj$par[5],reps = MCMC_setting$reps,MCMC_setting$gridsize)
  #q =  MCMC_obj$logMultiNorm
  MCMC_obj$logMultiNorm = log_like_traj(MCMC_obj$LatentTraj,MCMC_obj$Ode_Traj_coarse,MCMC_obj$FT,MCMC_setting$gridsize,MCMC_setting$t_correct)
  #print(MCMC_obj$logMultiNorm - q)
  MCMC_obj$coalLog = coal_loglik(MCMC_setting$Init,MCMC_obj$LatentTraj,MCMC_setting$t_correct,MCMC_obj$par[5],MCMC_setting$gridsize)
  return(list(MCMC_obj=MCMC_obj))
}



#' @title initialiaze MCMC
#' @return an MCMC_object
#'
#'


MCMC_initialize = function(MCMC_setting){ #, prior_par = c(10,20,-2.3,200,40)){
  alpha1 = exp(rgamma(1,MCMC_setting$b1,MCMC_setting$a1))
  alpha2 = exp(rgamma(1,MCMC_setting$b2,MCMC_setting$a2))
  S = MCMC_setting$N * alpha1 / (alpha1 + alpha2 + 1)
  I = MCMC_setting$N * alpha2 / (alpha1 + alpha2 + 1)
  state = c(X = S, Y = I)

  logMultiNorm = NaN
  ########
  s1 = runif(1,1,10)
  s2 = exp(rnorm(1,MCMC_setting$c1,0.4))
  theta1 = s1 * s2 / N
  theta2 = s2
  param = c(theta1 = theta1, theta2 = theta2)

  #Ode_Traj_thin <- ODE(log(state),times,param)

  Ode_Traj_thin = ODE(log(state), MCMC_setting$times,
                      param)

  Ode_Traj_coarse = Ode_Traj_thin[MCMC_setting$gridset,]


  FT = SIR_log_KOM_Filter2(Ode_Traj_thin,theta1,theta2,MCMC_setting$gridsize)


  Latent = Traj_sim(log(state),Ode_Traj_coarse,FT)
  LatentTraj = Latent$SimuTraj
  lambda = 500
  coalLog = coal_loglik(MCMC_setting$Init,LatentTraj,MCMC_setting$t_correct,lambda,MCMC_setting$gridsize)
  # lambda = exp(rgamma(1,16,4))

  logMultiNorm = Latent$loglike
  LogAlpha1 = dgamma(log(alpha1),MCMC_setting$b1,MCMC_setting$a1,log = T)
  LogAlpha1 = dgamma(log(alpha1),MCMC_setting$b1,MCMC_setting$a2,log = T)
  LogS2 = dnorm(log(s2),MCMC_setting$c1,0.4,log = T)
  LogLambda = dgamma(log(lambda),MCMC_setting$d1,MCMC_setting$d2,log = T)
  MCMC_obj = list(par = c(S,I,theta1,theta2,lambda),LatentTraj = LatentTraj, logMultiNorm = logMultiNorm,
                  Ode_Traj_coarse = Ode_Traj_coarse, FT = FT, coal_loglik = coalLog)
  ##########
  while(is.nan(logMultiNorm)){
    s1 = runif(1,1,10)
    s2 = exp(rnorm(1,MCMC_setting$c1,0.4))
    theta1 = s1 * s2 / N
    theta2 = s2
    state =  state = c(X = S, Y = I)
    param = c(theta1 = theta1, theta2 = theta2)

    Ode_Traj_thin = ODE(log(state), MCMC_setting$times,
                        param)

    Ode_Traj_coarse = Ode_Traj_thin[MCMC_setting$gridset,]


    FT = SIR_log_KOM_Filter2(Ode_Traj_thin,theta1,theta2,MCMC_setting$gridsize)


    Latent = Traj_sim(log(state),Ode_Traj_coarse,FT)
    LatentTraj = Latent$SimuTraj

    # lambda = exp(rgamma(1,16,4))

    logMultiNorm = Latent$loglike
  }

  coalLog = coal_loglik(MCMC_setting$Init,LatentTraj,MCMC_setting$t_correct,lambda,MCMC_setting$gridsize)
  MCMC_obj = list(par = c(S,I,theta1,theta2,lambda),LatentTraj = LatentTraj, logMultiNorm = logMultiNorm,
                  Ode_Traj_coarse = Ode_Traj_coarse, FT = FT, coalLog = coalLog)
  # MCMC_para = matrix(nrow = niter,ncol = 2)
  cat("Initialize MCMC \n")
  return(MCMC_obj)
}


MCMC_setup = function(coal_obs,times,t_correct,N,gridsize=50,niter = 1000,burn = 500,thin = 5,
                      a1 = 10, a2 = 20,b1 = 60, b2= 60, c1=-2.3,d1 = 200, d2 =40){
  gridset = seq(1,length(times),by=gridsize)
  grid = times[gridset]
  Init = coal_lik_init(coal_obs$samp_times, coal_obs$n_sampled, coal_obs$coal_times, times, t_correct)
  MCMC_setting = list(Init = Init,times = times,t_correct = t_correct,N = N,
                      gridsize=gridsize,gridset = gridset, niter = niter,burn = burn,thin = thin,
                      a1 = a1, a2 = a2,b1 =b1, b2 = b2, c1= c1,d1 = d1, d2 = d2,
                      reps=1)
  cat("MCMC set up ready \n")
  return(MCMC_setting)
}


MCMC_initialize2 = function(MCMC_setting){ #, prior_par = c(10,20,-2.3,200,40)){
  #alpha1 = exp(rgamma(1,MCMC_setting$b1,MCMC_setting$a1))
  #alpha2 = exp(rgamma(1,MCMC_setting$b2,MCMC_setting$a2))
  #S = MCMC_setting$N * alpha1 / (alpha1 + alpha2 + 1)
  #I = MCMC_setting$N * alpha2 / (alpha1 + alpha2 + 1)
  alpha1 = exp(rnorm(1,MCMC_setting$b1,MCMC_setting$a1))

  S = MCMC_setting$N * alpha1 / (alpha1  + 1)
  I = MCMC_setting$N / (alpha1 + 1)
  state = c(X = S, Y = I)

  logMultiNorm = NaN
  ########
  while(is.nan(logMultiNorm)){

    s1 = runif(1,1,10)
    s2 = exp(rnorm(1,MCMC_setting$c1,0.4))
    theta1 = s1 * s2 / MCMC_setting$N
    theta2 = s2
    #s2 = theta2
    #s1 = theta1 * MCMC_setting$N / s2
    param = c(theta1 = theta1, theta2 = theta2)

    #Ode_Traj_thin <- ODE(log(state),times,param)

    Ode_Traj_thin = ODE(log(state), MCMC_setting$times,
                        param)

    Ode_Traj_coarse = Ode_Traj_thin[MCMC_setting$gridset,]


    FT = SIR_log_KOM_Filter2(Ode_Traj_thin,theta1,theta2,MCMC_setting$gridsize)


    Latent = Traj_sim(log(state),Ode_Traj_coarse,FT)
    LatentTraj = Latent$SimuTraj
    lambda = 500
    coalLog = coal_loglik(MCMC_setting$Init,LatentTraj,MCMC_setting$t_correct,lambda,MCMC_setting$gridsize)
    # lambda = exp(rgamma(1,16,4))

    logMultiNorm = Latent$loglike
  }
  LogAlpha1 = dnorm(log(alpha1),MCMC_setting$b1,MCMC_setting$a1,log = T)
  #LogAlpha2 = dgamma(log(alpha2),MCMC_setting$b1,MCMC_setting$a2,log = T)
  LogS2 = dnorm(log(s2),MCMC_setting$c1,0.4,log = T)
  LogLambda = dgamma(log(lambda),MCMC_setting$d1,MCMC_setting$d2,log = T)
  MCMC_obj = list(par = c(S,I,theta1,theta2,lambda),LatentTraj = LatentTraj, logMultiNorm = logMultiNorm,
                  Ode_Traj_coarse = Ode_Traj_coarse, FT = FT, coalLog = coalLog,
                  LogAlpha1 = LogAlpha1, LogS2 = LogS2, LogLambda = LogLambda)
  ##########
  # MCMC_para = matrix(nrow = niter,ncol = 2)
  cat("Initialize MCMC \n")
  return(MCMC_obj)
}



SIR_LNA_MCMC = function(coal_obs,times,t_correct,N,gridsize=1000,niter = 1000,burn = 0,thin = 5,
                        a1 = 10, a2 = 20, b1 = 60 , b2 = 60,c1=-2.3,d1 = 200, d2 =40){
  MCMC_setting = MCMC_setup(coal_obs,times,t_correct,N,gridsize,niter,burn,thin,
                            a1, a2,b1,b2,c1,d1, d2)
  MCMC_obj = MCMC_initialize2(MCMC_setting)
  # MCMC_obj$LatentTraj = Traj
  # MCMC_obj$LogS2 =
  # MCMC_obj$LogLambda
  # MCMC_obj$par[1]=9999
  # MCMC_obj$par[2]=1000
  params = matrix(nrow = niter, ncol = 5)
  l = numeric(niter)
  l1 = l
  l2 = l
  l3 = l
  tjs = NULL
  for (i in 1 : MCMC_setting$niter) {
    if (i %% 100 == 0) {
      print(i)
    }
    step1 = updateAlphas(MCMC_obj,MCMC_setting,i)#  MCMC_obj = step1$MCMC_obj
    MCMC_obj = step1$MCMC_obj
    step2 = updateS1(MCMC_obj,MCMC_setting,i)
    # MCMC_obj = Like_Free(step2$MCMC_obj,MCMC_setting,i)$MCMC_obj
    MCMC_obj = step2$MCMC_obj
    step3 = updateS2(MCMC_obj,MCMC_setting,i)
    MCMC_obj = step3$MCMC_obj
    MCMC_obj = updateTraj(MCMC_obj,MCMC_setting,i)$MCMC_obj
    step4 = updateLambda(MCMC_obj,MCMC_setting,i)
    MCMC_obj = step4$MCMC_obj
    tjs = abind(tjs,MCMC_obj$LatentTraj,along = 3)
    params[i,] = MCMC_obj$par
    l[i] =  MCMC_obj$logMultiNorm #+ MCMC_obj$LogAlpha1 + MCMC_obj$LogAlpha2
    l1[i] = MCMC_obj$coalLog
    l2[i] = MCMC_obj$LogAlpha1
    l3[i] = MCMC_obj$LogAlpha1 + MCMC_obj$logMultiNorm + MCMC_obj$LogS2 + MCMC_obj$coalLog + MCMC_obj$LogLambda
  }
  return(list(par = params,Trajectory = tjs,l=l,l1=l1,l2 = l2, l3 =l3))
}
