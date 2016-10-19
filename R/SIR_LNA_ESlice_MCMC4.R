
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
  q = 10
  
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
      
      LatentTraj_new = ESlice(LatentTraj,Ode_Traj_coarse_new,FT_new,log(state_new),
                              Init,t_correct,lambda,reps = 2,gridsize)
      
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
      
      
      LatentTraj_new = ESlice(LatentTraj,Ode_Traj_coarse_new,FT_new,log(state),
                              Init,t_correct,lambda,reps = 5,gridsize)
      
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
                              Init,t_correct,lambda,reps = 5,gridsize)
      
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
      a = min(c(exp(dnorm(s2_new,c1,0.4,log = T) +coalLog_new +(logMultiNorm_new - logMultiNorm)/q -
                      ( coalLog + dnorm(s2,c1,0.4,log = T))),1))
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
    
    L_joint[i + 1] = coalLog 
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