#' Simulate a coalescent time based on volz likelihood
#'
#' @param eff_pop Three column matrix for population on a time grid. The first column is time. The second colunmn is number of susceptible and the third is number of infected, with both in log-scale.
#' @param t_correct The time that we observe the first sampling time
#' @param samp_time The time for the sampling event, reverse time scale
#' @param n_sampled Number of sampled event corresponse to each sampling time
#' @param betaN A vector of infection rate beta at the same time grid as eff_pop
#'
#' @return A list of sampling time, number sampled, interval of time points and coalsecent time inverse to the time scale

volz_thin_sir <- function(eff_pop,t_correct,samp_times, n_sampled, betaN, ...)
{
  coal_times = NULL
  lineages = NULL

  curr = 1
  active_lineages = n_sampled[curr]
  time = samp_times[curr]

  traj = 1 / (betaN*eff_pop[,2] /((eff_pop[,3])/2))
  #traj[1]=0.00000000001
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

    time = time + rexp(1, 0.5*active_lineages*(active_lineages-1)*10)

    i = which.min(ts>=0) - 1
    while(time>ts[i]){
      i = i-1
      if(i==0) {
        print("i=0")
        print(length(coal_times))
        i = 1
        break
      }
    }
    thed = traj[i]
    if (curr < length(samp_times) && time >= samp_times[curr + 1])
    {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }
    else if (runif(1) <= 1/(10 * thed) )
    {
      time = min(c(time,t_correct))
      coal_times = c(coal_times, time)
      lineages = c(lineages, active_lineages)
      active_lineages = active_lineages - 1
    }
  }

  return(list(coal_times = coal_times, lineages = lineages,
              intercoal_times = c(coal_times[1], diff(coal_times)),
              samp_times = samp_times, n_sampled = n_sampled))
}



coalsim_thin_sir <- function(eff_pop,t_correct,samp_times, n_sampled, lambda=1, ...)
{
  coal_times = NULL
  lineages = NULL

  curr = 1
  active_lineages = n_sampled[curr]
  time = samp_times[curr]

  traj = eff_pop[,2]/lambda
  #traj[1]=0.00000000001
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

    time = time + rexp(1, 0.5*active_lineages*(active_lineages-1)*10)

    i = which.min(ts>=0) - 1
    while(time>ts[i]){
      i = i-1
      if(i==0) {
        print("i=0")
        i = 1
        break
      }
    }
    thed = traj[i]
    if (curr < length(samp_times) && time >= samp_times[curr + 1])
    {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }
    else if (runif(1) <= 1/(10 * thed) )
    {
      time = min(c(time,t_correct))
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

#' Generate the suffcient statistics for coalescent likelihood calculation
#'
#' @param Grid A vector of time points that we use to infer the trajectory. grid is in reverse time scale and should match the number the grid of trajectory. The first grid point corresponds to the first sampling event
#'
#' @return A list sufficient statistics for coalescent likelihood calculation:
#'
#'


coal_lik_init = function(samp_times, n_sampled, coal_times, grid, t_correct)
{
  ns = length(samp_times)
  nc = length(coal_times)
  samp_times = samp_times[ns:1]
  coal_times = coal_times[nc:1]

  #n0 = which.min(grid < t_correct)
  Tcoal = max(coal_times)

  n0 = which.min(grid < Tcoal)
  grid = grid[1:n0]


  #ng = length(grid)-1
  ng = n0
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
  grid_idx = c(0,cumsum(gridrep))
  if(gridrep[ng] == 0){
    grid_idx = grid_idx[1:(length(grid_idx)-1)]
  }
  return(list(t=t, l=l, C=C, D=D, y=y, count=count, gridrep=gridrep,
              ns=sum(n_sampled), nc=nc, ng=ng, rep_idx=rep_idx,
              args=list(samp_times=samp_times, n_sampled=n_sampled, coal_times=coal_times,
                        grid=grid),
              grid_idx = grid_idx))
}
