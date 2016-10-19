#' @title Eliptical slice sampling
#' 
#' @param f_cur Current trajectory
#' @param OdeTraj Ode trajectory of the given parameter
#' @param FTs
#' @param state Initial state
#' @param t_correct time correction for coalescent times
#' @param lambda scale parameter for effective population
#' @param reps Repeat times
#' 
ESlice = function(f_cur,OdeTraj,FTs,state,Init,t_correct,lambda=10,reps=1,gridsize = 100){
  for(count in 1:reps){
  # centered old trajectory without time
  f_cur_centered = f_cur[,2:3] - OdeTraj[,2:3]
  
  # simulate a new  trajectory
  v = Traj_sim(state,OdeTraj,FTs) # OdeTraj is not the originall one with high resolution
  # centered new trajectory without time
  v_traj = v$SimuTraj[,2:3] - OdeTraj[,2:3] 

#  v_traj = cbind(OdeTraj[,1],v_traj)

  u = runif(1,0,1)
  logy = coal_loglik(Init,f_cur,t_correct,lambda,gridsize) + log(u)
  
  theta = runif(1,0,2*pi)
  theta_min = theta - 2*pi
  theta_max = theta
  
  f_prime = f_cur_centered * cos(theta) + v_traj * sin(theta)
  newTraj = cbind(OdeTraj[,1],f_prime + OdeTraj[,2:3])
  if(is.na(sum(newTraj))){
    return(f_cur)
  }
  else{
    while(coal_loglik(Init,newTraj,t_correct,lambda,gridsize) <= logy){
   
      # shrink the bracket
     if(theta < 0){
      theta_min = theta
    }else{
      theta_max = theta
    }
    theta = runif(1,theta_min,theta_max)
    
    f_prime = f_cur_centered * cos(theta) + v_traj * sin(theta)
    newTraj = cbind(OdeTraj[,1],f_prime + OdeTraj[,2:3])

      f_cur = newTraj
    }
  }
  return(newTraj)
  }
}
##############


#'
#'
#'
#'
#'
#'
SIR_LNA_ESlice_MCMC =function(state,Init,times,gridsize=100,niter = 1000,burn = 500,thin = 5){
  
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
  #' prior
  #' s1:  s1 = theta1  log(1/s1) ~ gamma(50,5)
  #' propose function random walk Uniform[-0.5,0.5] on log space
  #' 
  #' s2: theta2/theta1   proposal function: random walk Uniform[-0.5,0.5]
  #' Normal(20000,2000)
  #' 
  #' lambda log(lambda) ~ gamma(16,4) propose random walk uniform [-0.5,0.5] on logspace
  
  
  
  # fix the inital state
  # initialize mcmc
  #SIR_prop = rdirichlet(1,c(1,0.1,0.01))
  #state = c(X = SIR_prop[1] * 100050, Y = SIR_prop[2] * 100050)
  #R = SIR_prop[3] * 100050
  countInf = 0
  difs = NULL
  rec = matrix(rep(0,3*niter),ncol=3,nrow=niter)
  s1 = exp(-rgamma(1,60,5))
  s2 = rnorm(1,20000,2000)
   theta1 = s1
  theta2 = s1 * s2
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
  lambda = 100
  logMultiNorm = Latent$loglike
  coalLog = coal_loglik(Init,LatentTraj,t_correct,lambda)
  MCMC_para = NULL
  MCMC_Traj = NULL
 # MCMC_para = matrix(nrow = niter,ncol = 2)
  cat("begin MCMC")
  L_joint = numeric(niter + 1)
  L_joint[1] = coalLog + logMultiNorm + dgamma(-log(s1),60,5,log = T) +
    dnorm(s2,20000,2000,log = T)
  
  for(i in 1:niter){
  
    ############  
    # update the state
#    S_new = state[1] * exp(runif(1,-0.3,0.3))
#    I_new = state[2] * exp(runif(1,-0.3,0.3))
#   
#     SIR_prop_new = matrix(c(S_new,I_new,R)/sum(c(S_new,I_new,R)),nrow=1)
#    state_new = SIR_prop_new[1:2] * 100050
#    R_new = SIR_prop_new[3] * 100050
    
#    Ode_Traj_thin_new <- ODE(log(state_new),times, 
#                             param)
    
#    FT_new = SIR_log_KOM_Filter2(Ode_Traj_thin_new,theta1,theta2,gridsize)
    
#    Ode_Traj_coarse_new = Ode_Traj_thin_new[gridset,]
    
 #   LatentTraj_new = ESlice(LatentTraj,Ode_Traj_coarse_new,FT_new,log(state_new),
 #                           Init,t_correct,lambda,reps = 5)
    
    
#    logMultiNorm_new = log_like_traj(LatentTraj_new,Ode_Traj_coarse_new,FT_new,gridsize)
    
#    coalLog_new = coal_loglik(Init,LatentTraj_new,t_correct,lambda)
    
#    if(is.nan(logMultiNorm_new)){
 #     logMultiNorm_new = -Inf
#      countInf = countInf + 1
#    }
#    a = min(exp(ddirichlet(SIR_prop_new,c(1,0.1,0.01),log = T) +
#                  coalLog_new + logMultiNorm_new -
#                  (logMultiNorm + coalLog+ ddirichlet(SIR_prop,c(1,0.1,0.01),log = T) )
#    ),1)
    
#    if(runif(1,0,1) < a){
 #     state=state_new
  #    LatentTraj = LatentTraj_new
  #    Ode_Traj_coarse = Ode_Traj_coarse_new
  #    logMultiNorm = logMultiNorm_new
  #    FT = FT_new
  #    coalLog = coalLog_new
  #    R = R_new
  #  }
    
    
######################    
    
    
    
    
  # Sample s1
      # propose theta1_new condition on theta1
   
    s1_new = ifelse(i<burn,s1 * exp(runif(1,-1,1)), s1 + runif(1,-1,1) * 0.00000002) 
    theta1_new = s1_new 
    theta2_new = s1_new * s2
    param_new = c(theta1 = theta1_new, theta2 = theta2_new)
 #   print(param_new)
#    Ode_Traj_thin_new <- ode(y = log(state), times = times, 
 #                        func = SIR.log.ode2, parms = param_new)
    Ode_Traj_thin_new <- ODE(log(state),times, 
                       param_new)
      
    FT_new = SIR_log_KOM_Filter2(Ode_Traj_thin_new,theta1_new,theta2_new,gridsize)
    
    Ode_Traj_coarse_new = Ode_Traj_thin_new[gridset,]
    
    if(i < burn){
    LatentTraj_new = ESlice(LatentTraj,Ode_Traj_coarse_new,FT_new,log(state),
                        Init,t_correct,lambda,reps = 1)
    
    
    logMultiNorm_new = log_like_traj(LatentTraj_new,Ode_Traj_coarse_new,FT_new,gridsize)
    
    coalLog_new = coal_loglik(Init,LatentTraj_new,t_correct,lambda)
    }else{
      LatentTraj_new = LatentTraj
      logMultiNorm_new = log_like_traj2(LatentTraj,times,log(state),theta1_new,theta2_new,gridsize)
      coalLog_new = coalLog
    }
    if(is.nan(logMultiNorm_new)){
      countInf = countInf + 1
      a = -1 
    }else{
      if(i<burn){
     a = min(exp(dgamma(-log(s1_new),25,2.5,log = T) +
            coalLog_new + logMultiNorm_new -
    (logMultiNorm + coalLog+ dgamma(-log(s1),25,2.5,log=T) )
    ),1)
      }
      else{
        a = min(exp(dgamma(-log(s1_new),25,2.5,log = T) +  logMultiNorm_new -
                      (logMultiNorm + dgamma(-log(s1),25,2.5,log=T) )
        ),1)
      }
    }
   
    
    if(runif(1,0,1) < a){
      rec[i,1] = 1
      theta1 = theta1_new
      theta2 = theta2_new
      LatentTraj = LatentTraj_new
      Ode_Traj_coarse = Ode_Traj_coarse_new
      logMultiNorm = logMultiNorm_new
      param = param_new
      s1=s1_new
      FT = FT_new
      coalLog = coalLog_new
    }
    
    LatentTraj= ESlice(LatentTraj,Ode_Traj_coarse,FT,log(state),
                       Init,t_correct,lambda,reps = 1)
    logMultiNorm = log_like_traj2(LatentTraj,times,log(state),theta1,theta2,gridsize)
    
    # Sample theta2
      # propose theta_2_new condition on theta2
   
    s2_new = s2 + ifelse(i<burn,rnorm(1,0,2000),rnorm(1,0,200))
    theta2_new = s1 * s2_new

    param_new = c(theta1 = theta1, theta2 = theta2_new)
#    print(param_new)
#    Ode_Traj_thin_new <- ode(y = log(state), times = times, 
#                            func = SIR.log.ode2, parms = param_new)
    Ode_Traj_thin_new <- ODE(log(state), times, 
                      param_new)
        
    FT_new = SIR_log_KOM_Filter2(Ode_Traj_thin_new,theta1,theta2_new,gridsize)
    
    
    Ode_Traj_coarse_new = Ode_Traj_thin_new[gridset,]
    
    if(i < burn){
    LatentTraj_new = ESlice(LatentTraj,Ode_Traj_coarse_new,FT_new,log(state),
                            Init,t_correct,lambda,reps = 1)
    
    logMultiNorm_new = log_like_traj(LatentTraj_new,Ode_Traj_coarse_new,FT_new,gridsize)
    
    coalLog_new = coal_loglik(Init,LatentTraj_new,t_correct,lambda)
    }else{
    LatentTraj_new = LatentTraj    
    logMultiNorm_new = log_like_traj2(LatentTraj,times,log(state),theta1,theta2_new,gridsize)
    coalLog_new = coalLog
    }
     if(is.nan(logMultiNorm_new)){
      countInf = countInf + 1
      a = -1
     }else{
       if(i<burn){
    a = min(c(exp(dnorm(s2_new,20000,2000,log = T) +coalLog_new +logMultiNorm_new - 
                    (logMultiNorm+ coalLog + dnorm(s2,20000,2000,log = T))),1))
       }else{
         a = min(c(exp(dnorm(s2_new,20000,2000,log = T) +logMultiNorm_new - 
                         (logMultiNorm + dnorm(s2,20000,2000,log = T))),1))
         # difs = c(difs,logMultiNorm_new - 
              #       log_like_traj2(LatentTraj,times,log(state),theta1,theta2,gridsize)) 
         }
     }
    
    
    
    if(runif(1,0,1) < a){
      rec[i,2] = 1
      theta2 = theta2_new
      Ode_Traj_coarse = Ode_Traj_coarse_new
      logMultiNorm = logMultiNorm_new
      param = param_new
      FT = FT_new
      s2 = s2_new
     coalLog = coalLog_new
      LatentTraj = LatentTraj_new
    }
    LatentTraj= ESlice(LatentTraj,Ode_Traj_coarse,FT,log(state),
                            Init,t_correct,lambda,reps = 1)
    logMultiNorm = log_like_traj2(LatentTraj,times,log(state),theta1,theta2,gridsize)
  # Sample the trajectory using Elipitical Slice Sampling
#    LatentTraj = ESlice(LatentTraj,Ode_Traj_coarse,FT,log(state),Init,t_correct,lambda,reps = 10)
  
    
  # Sample the scale parameter lambda
    lambda_new = ifelse(i<burn, lambda * exp(runif(1,-0.3,0.3)),lambda + 5 * runif(1,-1,1))
    coalLog_new = coal_loglik(Init,LatentTraj,t_correct,lambda_new)
    a = min(c(exp(coalLog_new - coalLog + dgamma(log(lambda_new),200,40,log = T) - 
                    dgamma(log(lambda),200,40,log=T)), 1))
        if(runif(1,0,1) < a){
          rec[i,3] = 1
      coalLog = coalLog_new
      lambda = lambda_new
    }
   # print(lambda)
  
  # store the MCMC samples 
    
    L_joint[i + 1] = coalLog + logMultiNorm + dgamma(-log(s1),60,5,log = T) +
      dnorm(s2,20000,2000,log = T)
    if(i %% 100 == 0) {
    cat(floor(i)," iterations complete","\n")
    }
    if(i>burn && i %% thin ==1 ){
    MCMC_para = cbind(MCMC_para,c(theta1,theta2,lambda))
    MCMC_Traj = abind(MCMC_Traj,LatentTraj,along = 3)   
    }
  }
  print(paste("num of Inf",countInf))
  return(list(par=MCMC_para,Traj=MCMC_Traj,OdeTraj = Ode_Traj_coarse,
              FT=FT,L=L_joint,
              records = rec,dif=difs))
}
###########################








#################

SIR_LNA_ESlice_MCMC2 =function(state,Init,times,t_correct,gridsize=50,niter = 1000,burn = 500,thin = 5){
  
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
  #' prior
  #' s1:  s1 = theta1  log(1/s1) ~ gamma(50,5)
  #' propose function random walk Uniform[-0.5,0.5] on log space
  #' 
  #' s2: theta2/theta1   proposal function: random walk Uniform[-0.5,0.5]
  #' Normal(20000,2000)
  #' 
  #' lambda log(lambda) ~ gamma(16,4) propose random walk uniform [-0.5,0.5] on logspace
  #'
  #'
  #' log(p1/p3) ~ gamma(60,10)
  #' log(p2/p3) ~ gamma(60,20)
  
  
  # fix the inital state
  # initialize mcmc
  
  alpha1 = exp(rgamma(1,60,10))
  alpha2 = exp(rgamma(1,60,20))
  S = (11000) * alpha1 / (alpha1 + alpha2 + 1)
  I = (11000) * alpha2 / (alpha1 + alpha2 + 1)
  state = c(X = S, Y = I)
  
  countInf = 0
  difs = NULL
  rec = matrix(rep(0,5*niter),ncol=5,nrow=niter)
  logMultiNorm = NaN
 ########
  s1 = exp(-rgamma(1,50,6))
 #s1 = 0.00002
  # s2 = 2500
  s2 = rnorm(1,2500,500)
  theta1 = s1
  theta2 = s1 * s2
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
  s1 = exp(-rgamma(1,50,6))

   s2 = rnorm(1,2500,500)
  theta1 = s1
  theta2 = s1 * s2
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
  coalLog = coal_loglik(Init,LatentTraj,t_correct,lambda)
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
        state_new = c(X = (11000) * alpha1_new / (alpha1_new + alpha2_new + 1),
                  Y = (11000) * alpha2_new / (alpha1_new + alpha2_new + 1))
        Ode_Traj_thin_new <- ODE(log(state_new),times, 
                                param)
    
       Ode_Traj_coarse_new = Ode_Traj_thin_new[gridset,]
      
       FT_new = SIR_log_KOM_Filter2(Ode_Traj_thin_new,theta1,theta2,gridsize)
       if(i<burn){
           
       LatentTraj_new = ESlice(LatentTraj,Ode_Traj_coarse_new,FT_new,log(state_new),
                               Init,t_correct,lambda,reps = 2)
    
    
     logMultiNorm_new = log_like_traj(LatentTraj_new,Ode_Traj_coarse_new,FT_new,gridsize)
    
        coalLog_new = coal_loglik(Init,LatentTraj_new,t_correct,lambda)
       }else{
         logMultiNorm_new = log_like_traj2(LatentTraj,times,log(state_new),theta1,theta2,gridsize)
         coalLog_new = coalLog
         }
        if(is.nan(logMultiNorm_new)){
         logMultiNorm_new = -Inf
          countInf = countInf + 1
        }
        
       a = min(exp(dgamma(log(alpha1_new),60,10,log = T) + dgamma(log(alpha2_new),60,20,log = T) +
                     coalLog_new + (logMultiNorm_new - logMultiNorm)/1 - 
                      ( coalLog+ dgamma(log(alpha1),60,10,log = T) +
                         dgamma(log(alpha2),60,20,log = T))
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
                          Init,t_correct,lambda,reps = 1)
       logMultiNorm = log_like_traj2(LatentTraj,times,log(state),theta1,theta2,gridsize)
       
    
    ######################    
    
    
    
    
    # Sample s1
    # propose theta1_new condition on theta1
    
    s1_new = s1 * exp(runif(1,-0.5,0.5))
    theta1_new = s1_new 
    theta2_new = s1_new * s2
    param_new = c(theta1 = theta1_new, theta2 = theta2_new)
    #   print(param_new)
    #    Ode_Traj_thin_new <- ode(y = log(state), times = times, 
    #                        func = SIR.log.ode2, parms = param_new)
    Ode_Traj_thin_new <- ODE(log(state),times, 
                             param_new)
    
    Ode_Traj_coarse_new = Ode_Traj_thin_new[gridset,]
   
    FT_new = SIR_log_KOM_Filter2(Ode_Traj_thin_new,theta1_new,theta2_new,gridsize)
    if(i < burn){
     
      
      LatentTraj_new = ESlice(LatentTraj,Ode_Traj_coarse_new,FT_new,log(state),
                              Init,t_correct,lambda,reps = 5)
      
      logMultiNorm_new = log_like_traj(LatentTraj_new,Ode_Traj_coarse_new,FT_new,gridsize)
      
      coalLog_new = coal_loglik(Init,LatentTraj_new,t_correct,lambda)
    }else{
      #ttt = Traj_sim(log(state),Ode_Traj_coarse_new,FT_new)
      #LatentTraj_new = ttt$SimuTraj
      #logMultiNorm_new = ttt$loglike
      #coalLog_new = coal_loglik(Init,LatentTraj_new,t_correct,lambda)
      logMultiNorm_new = log_like_traj2(LatentTraj,times,log(state),theta1_new,theta2_new,gridsize)
      coalLog_new = coalLog
      }
    if(is.nan(logMultiNorm_new)){
      countInf = countInf + 1
      a = -1 
    }else{
        a = min(exp(dgamma(-log(s1_new),50,6,log = T) +
                      coalLog_new + (logMultiNorm_new - logMultiNorm)/1 - 
                      ( coalLog+ dgamma(-log(s1),50,6,log=T) )
        ),1)
    }
    
    
    if(runif(1,0,1) < a){
      rec[i,3] = 1
      theta1 = theta1_new
      theta2 = theta2_new
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
                       Init,t_correct,lambda,reps = 1)
    logMultiNorm = log_like_traj2(LatentTraj,times,log(state),theta1,theta2,gridsize)
    
  
    # Sample theta2
    # propose theta_2_new condition on theta2
    
    s2_new = s2 + ifelse(i<burn,rnorm(1,0,500),rnorm(1,0,200))
    theta2_new = s1 * s2_new
    
    param_new = c(theta1 = theta1, theta2 = theta2_new)
  
    
    Ode_Traj_thin_new <- ODE(log(state), times, 
                             param_new)
    
    Ode_Traj_coarse_new = Ode_Traj_thin_new[gridset,]
    FT_new = SIR_log_KOM_Filter2(Ode_Traj_thin_new,theta1,theta2_new,gridsize)
    if(i < burn){
      
      LatentTraj_new = ESlice(LatentTraj,Ode_Traj_coarse_new,FT_new,log(state),
                              Init,t_correct,lambda,reps = 5)
      
     logMultiNorm_new = log_like_traj(LatentTraj_new,Ode_Traj_coarse_new,FT_new,gridsize)
      
      coalLog_new = coal_loglik(Init,LatentTraj_new,t_correct,lambda)
    }else{
      
      logMultiNorm_new = log_like_traj2(LatentTraj,times,log(state),theta1,theta2_new,gridsize)
      coalLog_new = coalLog
    
      }
    if(is.nan(logMultiNorm_new)){
     countInf = countInf + 1
      a = -1
    }else{
        a = min(c(exp(dnorm(s2_new,2500,500,log = T) +coalLog_new +(logMultiNorm_new - logMultiNorm)/1 -
                        ( coalLog + dnorm(s2,2500,500,log = T))),1))
    }
    
    
    if(runif(1,0,1) < a){
      rec[i,4] = 1
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
                       Init,t_correct,lambda,reps = 1)
    logMultiNorm = log_like_traj2(LatentTraj,times,log(state),theta1,theta2,gridsize)
 
    
    
       # Sample the trajectory using Elipitical Slice Sampling
       # LatentTraj = ESlice(LatentTraj,Ode_Traj_coarse,FT,log(state),Init,t_correct,lambda,reps = 10)
    
        #logMultiNorm = log_like_traj2(LatentTraj,times,log(state),theta1,theta2,gridsize)
    # Sample the scale parameter lambda
   
     lambda_new = lambda * exp(runif(1,-0.3,0.3))
    coalLog_new = coal_loglik(Init,LatentTraj,t_correct,lambda_new)
    a = min(c(exp(coalLog_new - coalLog + dgamma(log(lambda_new),200,40,log = T) - 
                    dgamma(log(lambda),200,40,log=T)), 1))
    if(runif(1,0,1) < a){
      rec[i,5] = 1
      coalLog = coalLog_new
      lambda = lambda_new
    }
    # print(lambda)
    
    # store the MCMC samples 
    
    L_joint[i + 1] = coalLog + logMultiNorm 
    #+ dgamma(-log(s1),60,5,log = T) +
    #  dnorm(s2,2500,500,log = T)  
   #    dgamma(log(alpha1),60,10,log = T) + 
   #   dgamma(log(alpha2),60,20,log = T)
    if(i %% 100 == 0) {
      cat(floor(i)," iterations complete","\n")
    }
    if(i>burn && i %% thin ==1 ){
      MCMC_para = cbind(MCMC_para,c(state,theta1,theta2,lambda))
      MCMC_Traj = abind(MCMC_Traj,LatentTraj,along = 3)   
    }
  }
  print(paste("num of Inf",countInf))
  return(list(par=MCMC_para,Traj=MCMC_Traj,OdeTraj = Ode_Traj_coarse,
              FT=FT,L=L_joint,
              records = rec,dif=difs))
}
