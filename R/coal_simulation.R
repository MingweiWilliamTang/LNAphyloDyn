#'
#'
#'@param eff_pop a two dimensional vector of time and the corresponding trajectory of infectious population
#'@param t_correct the actual time point of sampled_time=0
#'
#'
#'
#'
coalsim_thin_sir <- function(eff_pop,t_correct,samp_times, n_sampled, lambda=1, ...)
{
  coal_times = NULL
  lineages = NULL

  curr = 1
  active_lineages = n_sampled[curr]
  time = samp_times[curr]

  traj = eff_pop[,2]/lambda
  ts = (t_correct-eff_pop[,1])
  active_lineages = n_sampled[curr]
  # time = 0
  i = which.min(ts>=0) - 1

  while (time <= max(samp_times) || active_lineages > 1)
  {
    if (active_lineages == 1)
    {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }

    time = time + rexp(1, 0.5*active_lineages*(active_lineages-1))

    if(i==0) break
    while(time>ts[i]){
      i = i-1
      if(i==0) {
        i = 1
        break
      }
    }
    thed = (traj[i] + traj[i+1])/2
    if (curr < length(samp_times) && time >= samp_times[curr + 1])
    {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }
    else if (runif(1) <= 1/thed)
    {
      coal_times = c(coal_times, time)
      lineages = c(lineages, active_lineages)
      active_lineages = active_lineages - 1
    }
  }

  return(list(coal_times = coal_times, lineages = lineages,
              intercoal_times = c(coal_times[1], diff(coal_times)),
              samp_times = samp_times, n_sampled = n_sampled))
}


coalsim_thin_sir_Hom = function(eff_pop,t_correct,n_sampled,lambda){
  coal_times = NULL
  lineages = NULL

  traj = eff_pop[,2]/lambda
  ts = (t_correct-eff_pop[,1])
  active_lineages = n_sampled
  time = 0
  i = which.min(ts>=0) - 1
  while (active_lineages > 1)
  {
    if(i==0) break
    time = time + rexp(1, 0.5*active_lineages*(active_lineages-1))

    while(time>ts[i]){
      i = i-1
      if(i==0) {
        i = 1
        break
      }
    }
    thred = (traj[i] + traj[i+1]) / 2
    if (runif(1) <= (1/ ( thred)) )
    {
      coal_times = c(coal_times, time)
      lineages = c(lineages, active_lineages)
      active_lineages = active_lineages - 1
    }
  }
  return(list(coal_times = coal_times, lineages = lineages,
              intercoal_times = c(coal_times[1], diff(coal_times)),
              samp_times = 0, n_sampled = n_sampled))

}



coal_loglik_hom = function(eff_pop,t_correct,gene,lambda){
  dt = eff_pop[2,1] - eff_pop[1,1]
  n = gene$lineages[1]
  traj = eff_pop[,2]
  ts = t_correct - eff_pop[,1]
  time = 0
  i = which.min(ts>=0) - 1
  j = 1
  loglik = 0
  while(n>1){
    # determine the time on the grid
    Ck = 0.5 * n * (n-1)
    time = gene$coal_times[j]
    integr = dt / (traj[i] * lambda)
    while(time > ts[i] ){
      i = i - 1
      integr = integr + dt/ (traj[i] * lambda)
    }
    loglik = loglik + log(Ck) - log(lambda) - log(traj[i]) - Ck * integr
    j = j+1
    n=n-1
  }
  return(loglik)
}



gene2 = coalsim_thin_sir(EP,45,samp_times = seq(0,30,length.out=50),rep(2,50),100)

EP2=out2[,c(1,3)]
gene2 = coalsim_thin_sir_Hom(EP2,80,200,1)

plot(EP[,1],EP[,2],type="l",ylab="infected population",xlab = "coalesent time")
points(45-gene2$coal_times,rep(100,length(gene2$coal_times)),pch=3,cex=0.5)

lines(90-res_BNPR$data$time,exp(res_BNPR$data$E_log))


EP = out[,c(1,3)]
gene = coalsim_thin_sir_Hom(EP,80,200,1)
plot(EP[,1],EP[,2],type="l",ylab="infected population",xlab = "coalesent time")
points(80-gene$coal_times,rep(100,length(gene$coal_times)),pch=3,cex=0.5)

# scaled by the population size
EP3 = out[,c(1,3)]
gene3 = coalsim_thin_sir(EP3,90,samp_times = seq(0,70,length.out=50),rep(1,50),10)
plot(EP3[,1],EP3[,2],type="l",ylab="infected population",xlab = "coalesent time")
points(90-gene3$coal_times,rep(50,length(gene3$coal_times)),pch=3,cex=0.5)

s2 = coal_lik_init(seq(0,70,length.out=50),rep(1,50),gene3$coal_times,seq(0,100,by=0.5),t_correct = 90)


EP = out2[,c(1,3)]

gene2 = coalsim_thin_sir(EP,45,samp_times = seq(0,30,length.out=50),rep(2,50),100)
plot(EP[,1],EP[,2],type="l",ylab="infected population",xlab = "coalesent time")
points(45-gene2$coal_times,rep(1000,length(gene2$coal_times)),pch=3,cex=0.5)
init2 = coal_lik_init(seq(0,30,length.out = 50),rep(1,50),gene2$coal_times,
                      seq(0,50,by=0.2),t_correct = 45)


coal_loglik(s2,out.log.coarse,90,10)
coal_loglik(s2,tjj,90,10)

res_BNPR = BNPR(data = gene, lengthout = 100)
plot_BNPR(res_BNPR, traj = exp, main="BNPR", ylim = c(1, 700))
nsamp = 500
nburnin = 100
res_MCMC = mcmc_sampling(data = gene, alg = 'splitHMC', nsamp = nsamp,
                         nburnin = nburnin)
plot_BNPR(BNPR_out = res_BNPR, traj = exp_traj, col = "blue", traj_col = "black")
lines(res_MCMC$med_fun, pch="", col='red', lwd=1)
lines(res_MCMC$low_fun, pch="", col='red', lwd=1, lty=1)
lines(res_MCMC$up_fun,  pch="", col='red', lwd=1, lty=1)
legend('topleft',c('Truth','BNPR',"splitHMC"),
       col=c('black','blue','red'),lwd=c(1,2.5,2.5),bty='n', lty=c(2,1,1))

traj = exp_traj
samp_times = 0
n_sampled  = 100
gene = coalsim(samp_times = samp_times, n_sampled = n_sampled, traj = traj, lower_bound = 1/20)


############################



# f is log effective samplesize


coal_loglik = function(init, f, t_correct,lambda,grad=FALSE)
{
  n0 = which.min(f[,1] < t_correct)
  f = f[(n0-1):1,3]

  if (init$ng != length(f))
    stop(paste("Incorrect length for f; should be", init$ng))



  f = rep(f, init$gridrep)
  llnocoal = init$D * init$C * exp(-f)*lambda
  if (!grad)
  {
    lls = init$y * (-f+log(lambda)) - llnocoal
    #print(lls)

    ll = sum(lls[!is.nan(lls)])

    return(ll)
  }
  else
  {
    dll = apply(init$rep_idx,1,function(idx)sum(-init$y[idx[1]:idx[2]]+llnocoal[idx[1]:idx[2]])) # gradient of log-likelihood wrt f_midpts

    return(dll)
  }
}

coal_lik_init = function(samp_times, n_sampled, coal_times, grid, t_correct)
{
  ns = length(samp_times)
  nc = length(coal_times)

  n0 = which.min(grid < t_correct)
  grid = grid[1:n0]
  ng = length(grid)-1

  if (length(samp_times) != length(n_sampled))
    stop("samp_times vector of differing length than n_sampled vector.")

  if (length(coal_times) != sum(n_sampled) - 1)
    stop("Incorrect length of coal_times: should be sum(n_sampled) - 1.")

  if (max(samp_times, coal_times) > max(grid))
    stop("Grid does not envelop all sampling and/or coalescent times.")

  t = sort(unique(c(samp_times, coal_times, grid)))
  l = rep(0, length(t))

  for (i in 1:ns)
    l[t >= samp_times[i]] = l[t >= samp_times[i]] + n_sampled[i]

  for (i in 1:nc)
    l[t >= coal_times[i]] = l[t >= coal_times[i]] - 1

  #print(l)

  if (sum((l < 1) & (t >= min(samp_times))) > 0)
    stop("Number of active lineages falls below 1 after the first sampling point.")

  mask = l > 0
  t = t[mask]
  l = head(l[mask], -1)

  gridrep = rep(0, ng)
  for (i in 1:ng)
    gridrep[i] = sum(t > grid[i] & t <= grid[i+1])

  C = 0.5 * l * (l-1)
  D = diff(t)

  y = rep(0, length(D))
  y[t[-1] %in% coal_times] = 1

  buckets = cut(x = samp_times, breaks = t,
                include.lowest = TRUE)
  tab <- aggregate(n_sampled ~ buckets, FUN = sum, labels = FALSE)
  count <- rep(0, length(D))
  count[as.numeric(tab$buckets)] <- tab$n_sampled
  count[head(t, -1) >= max(samp_times)] <- NA

  rep_idx = cumsum(gridrep)
  rep_idx = cbind(rep_idx-gridrep+1,rep_idx)

  return(list(t=t, l=l, C=C, D=D, y=y, count=count, gridrep=gridrep, ns=sum(n_sampled), nc=nc, ng=ng, rep_idx=rep_idx, args=list(samp_times=samp_times, n_sampled=n_sampled, coal_times=coal_times, grid=grid)))
}
