coal_init_stan = function(init,Dt,t_correct,gridsize,mu,N,period,prior_para = NULL){

  ng = length(init$grid_idx) - 1
  nc = length(init$D)
  dt = Dt / gridsize
  t0 = t_correct - max(init$t) + Dt
  t_new = seq(t0, t_correct, by = Dt)
  ts = seq(t0, t_correct, by = dt)
  ngrid = length(t_new) - 1
  n = length(ts)
  grid_inv = diff(init$grid_idx)[ng:1]
  tid = NULL
  for(i in 1:ng){
    tid = c(tid, rep(i,grid_inv[i]))
  }
  w = init$D[nc:1]
  y = init$y[nc:1]
  C = init$C[nc:1]
  return(list(n = n, ngrid = ngrid, gridsize = gridsize,period = period,N = N,
              ts = ts, tt = t_new, mu = mu,
              nc = nc, tids = tid, w = w, C = C, y = y))
}



coal_init_stan_changepoint = function(init,Dt,t_correct,gridsize,N,ct,nch,
                                      prior_para = list(r0=0,r1=10, a0 = 8, a1= 1,
                                                        g0 = 15, g1 = 1, l0=10,l1=3, chpr=2)){

  ng = length(init$grid_idx) - 1
  nc = length(init$D)
  dt = Dt / gridsize
  t0 = t_correct - max(init$t) + Dt
  t_new = seq(t0, t_correct, by = Dt)
  ts = seq(t0, t_correct, by = dt)
  ngrid = length(t_new) - 1
  n = length(ts)
  grid_inv = diff(init$grid_idx)[ng:1]
  tid = NULL
  for(i in 1:ng){
    tid = c(tid, rep(i,grid_inv[i]))
  }
  w = init$D[nc:1]
  y = init$y[nc:1]
  C = init$C[nc:1]
  return(list(ngrid = ngrid,tt = t_new,nch = nch,x_r = c(N,ct),
              nc = nc, tids = tid, w = w, C = C, y = y,
              r0=prior_para$r0,r1=prior_para$r1, a0 = prior_para$a0, a1= prior_para$a1,
              g0 = prior_para$g0, g1 = prior_para$g1, l0=prior_para$l0,l1=prior_para$l1,
              chpr=prior_para$chpr)
         )
}



SI_traj_initialize_stan = function(input,R0,gamma,mu,A,Alpha){
  state = numeric(2)
  state[1] = input$N * Alpha / (1 + Alpha)
  state[2] = input$N - state[1]
  beta = R0 * gamma / input$N
  traj = Traj_sim_SIR_BD_ez(state, input$ts, param = c(beta,gamma,0,0), input$gridsize, input$N,
                            max(input$ts), 20)$Simu
  plot(traj[,1],traj[,3],type="l")
  traj_list = list()
  for(i in 1:dim(traj)[1]){
    traj_list[[i]] = traj[i,2:3]
  }
  return(traj_list)
}


SI_traj_initialize_stan2 = function(input,R0,gamma,mu,A,Alpha,N,period = 20){
  state = numeric(2)
  state[1] = N * Alpha / (1 + Alpha)
  state[2] = N - state[1]
  beta = R0 * gamma / N
  traj = SIR_BD_ODE(state,input$tt,param = c(beta,gamma,mu,A),N,period)[2:length(input$tt),]
  plot(traj[,1],traj[,3],type="l")
  return(traj[,2:3])
}


para_init_stan_changepoint = function(input,N,chains=1){
  res = list()
  for(i in 1:chains){
    R0 = runif(1,1,10)
    gamma = exp(rnorm(1,2.5,1))
    Alpha = rnorm(1,8,2)
    ch = 0.2
    trajlist = SI_traj_initialize_stan2(input, R0, gamma, 0, 0, Alpha,N)
    res[[i]] = list(R0=R0,gamma=gamma,Alpha=Alpha,lambda = 5000,ch = ch,SI = trajlist)
  }
  return(res)
}


para_init_stan_changepoint2 = function(input,N,chains=1){
  res = list()
  for(i in 1:chains){
    R0 = runif(1,1,10)
    gamma = exp(rnorm(1,2.5,1))
    Alpha = rnorm(1,8,2)
    ch = rlnorm(input$nch,0,0.5)
    #trajlist = SI_traj_initialize_stan2(input, R0, gamma, 0, 0, Alpha,N)
    res[[i]] = list(R0=R0,gamma=gamma,Alpha=Alpha,lambda = 5000,ch = ch)
  }
  return(res)
}


para_init_stan = function(input,mu = 0.2,chains=1){
  res = list()
  for(i in 1:chains){
  R0 = runif(1,1,7)
  gamma = exp(rnorm(1,-3,0.4))
  Alpha = rnorm(1,4,0.5)
  A = runif(1,0,1)
  trajlist = SI_traj_initialize_stan2(input, R0, gamma, mu, A, Alpha)
  res[[i]] = list(R0=R0,gamma=gamma,Alpha=Alpha,A = A,SI = trajlist)
  }
  return(res)
}

stan_res_convert = function(fit,N){

  R0 = As.mcmc.list(fit, pars = "R0")[[1]]
  nrep = length(R0)
  SI_matrix = matrix(As.mcmc.list(fit,pars = "SI")[[1]], byrow = T, ncol = nrep)
  gamma1 = As.mcmc.list(fit, pars = "gamma")[[1]]
  A = As.mcmc.list(fit, pars = "A")[[1]]
  Alpha = as.numeric(As.mcmc.list(fit, pars = "Alpha")[[1]])
  ngrid = dim(SI_matrix)[1] / 2
  nrep = length(R0)
  LNATraj = array(dim = c(ngrid+1,2,nrep))

  for(i in 1:nrep){

    LNATraj[1,1,i] = N * Alpha[i] /(1 + Alpha[i])
    LNATraj[1,2,i] = N - LNATraj[1,1,i]

    for(j in 1:ngrid){
      LNATraj[j+1,,i] = SI_matrix[(2*j-1):(2*j),i]
    }
  }
  beta = R0 * gamma1 / N
  return(list(beta = beta, R0 = R0, gamma = gamma1, A = A, Alpha = Alpha,
              latentTraj = LNATraj))
}

stan_res_convert0 = function(fit,N){

  R0 = As.mcmc.list(fit, pars = "R0")[[1]]
  nrep = length(R0)
  SI_matrix = matrix(As.mcmc.list(fit,pars = "SI")[[1]], byrow = T, ncol = nrep)
  gamma1 = As.mcmc.list(fit, pars = "gamma")[[1]]
  A = As.mcmc.list(fit, pars = "A")[[1]]
  Alpha = as.numeric(As.mcmc.list(fit, pars = "Alpha")[[1]])
  ngrid = dim(SI_matrix)[1] / 2 - 1
  nrep = length(R0)
  LNATraj = array(dim = c(ngrid+1,2,nrep))

  for(i in 1:nrep){

   # LNATraj[1,1,i] = N * Alpha[i] /(1 + Alpha[i])
  #  LNATraj[1,2,i] = N - LNATraj[1,1,i]

    for(j in 1:(ngrid +1)){
      LNATraj[j,,i] = SI_matrix[(2*j-1):(2*j),i]
    }
  }
  beta = R0 * gamma1 / N
  return(list(beta = beta, R0 = R0, gamma = gamma1, A = A, Alpha = Alpha))
}

stan_res_convertC0 = function(fit,N){

  R0 = As.mcmc.list(fit, pars = "R0")[[1]]
  nrep = length(R0)
  SI_matrix = matrix(As.mcmc.list(fit,pars = "SIs")[[1]], byrow = T, ncol = nrep)
  gamma1 = As.mcmc.list(fit, pars = "gamma")[[1]]
 # A = As.mcmc.list(fit, pars = "A")[[1]]
  ch = As.mcmc.list(fit, pars = "ch")[[1]]
  lambda = As.mcmc.list(fit, pars = "lambda")[[1]]
  Alpha = as.numeric(As.mcmc.list(fit, pars = "Alpha")[[1]])
  ngrid = dim(SI_matrix)[1] / 2 - 1
  nrep = length(R0)
  LNATraj = array(dim = c(ngrid+1,2,nrep))

  for(i in 1:nrep){

    # LNATraj[1,1,i] = N * Alpha[i] /(1 + Alpha[i])
    #  LNATraj[1,2,i] = N - LNATraj[1,1,i]

    for(j in 1:(ngrid +1)){
      LNATraj[j,,i] = SI_matrix[(2*j-1):(2*j),i]
    }
  }
  beta = R0 * gamma1 / N
  return(list(R0 = R0, gamma = gamma1, ch = ch, lambda= lambda, Alpha = Alpha,LatentTraj = LNATraj))
}



Prevalence_Stan_data = function(y,tt,gridsize,N,mu,period=1){
  if(length(y) != length(tt)){
    stop('the size of observation not equal to size of time')
  }
  ngrid = length(y) - 1
  n = ngrid * gridsize + 1
  ts = numeric(0)
  for(i in 1:ngrid){
    ts = c(ts, seq(tt[i],tt[i+1],length.out = gridsize + 1)[-(gridsize + 1)])
  }
  ts = c(ts, tt[ngrid + 1])
  if(n != length(ts)){
    stop("grid not equal")
  }

  return(list(n=n, ngrid = ngrid, gridsize = gridsize, period = period, N = N,
              ts = ts, tt = tt, mu = mu, y = y))
}

Prevalence_Stan_Initialize = function(input, nchain = 1){
  res = list()
  for(i in 1:nchain){
    R0 = 5
    gamma = exp(rnorm(1,-2,0.2))
    A = runif(1,0,1)
    rho = rbeta(1,2,2)
    Alpha = exp(rnorm(1,0,2))
    trajlist = SI_traj_initialize_stan(input, R0, gamma, input$mu, A, Alpha)
    res[[i]] = list(R0=R0,gamma=gamma,Alpha=Alpha,A = A,SI = trajlist, rho = rho)
  }
  return(res)
}


