#' Simulate a coalescent time based on volz likelihood
#'
#' @param eff_pop Three column matrix for population on a time grid. The first column is time. The second colunmn is number of susceptible and the third is number of infected, with both in log-scale.
#' @param t_correct The time that we observe the first sampling time
#' @param samp_time The time for the sampling event, reverse time scale
#' @param n_sampled Number of sampled event corresponse to each sampling time
#' @param betaN A vector of infection rate beta at the same time grid as eff_pop
#'
#' @return A list of sampling time, number sampled, interval of time points and coalsecent time inverse to the time scale

volz_thin_sir <- function(eff_pop,t_correct,samp_times, n_sampled, betaN,lower = 1, ...)
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
      #print(curr)
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }

    time = time + rexp(1, 0.5*active_lineages*(active_lineages-1)/lower)

    #i = which.min(ts>=0) - 1
    while(time>ts[i]){
      i = i-1
      if(i==0) {
        print("i=0")
        print(length(coal_times))
        i = 1
        break
      }
    }
   # print(ts[i])
    thed = traj[i]
    if (curr < length(samp_times) && time >= samp_times[curr + 1])
    {
      curr = curr + 1
      #print(curr)
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }
    else if (runif(1) <= lower/(thed) )
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



coalsim_thin_sir <- function(eff_pop,t_correct,samp_times, n_sampled, lambda=1,lower = 0.1, ...)
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

    time = time + rexp(1, 0.5*active_lineages*(active_lineages-1)/lower)

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
    else if (runif(1) <= lower/( thed) )
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


coal_lik_init = function(samp_times, n_sampled, coal_times, grid)
{
  ns = length(samp_times)
  nc = length(coal_times)
  #samp_times = samp_times[ns:1]
  #coal_times = coal_times[nc:1]

  #n0 = which.min(grid < t_correct)
  Tcoal = max(coal_times)
  if(Tcoal > max(grid)){
    stop("coalescent time out of grid")
  }
  n0 = which.min(grid < Tcoal)
  grid = grid[1:n0]


  #ng = length(grid)-1
  ng = n0
  if (length(samp_times) != length(n_sampled))
    stop("samp_times vector of differing length than n_sampled vector.")

#  if (length(coal_times) != sum(n_sampled) - 1)
#    stop("Incorrect length of coal_times: should be sum(n_sampled) - 1.")

#  if (max(samp_times, coal_times) > max(grid))
 #   stop("Grid does not envelop all sampling and/or coalescent times.")

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
  print(ng)
  if(gridrep[ng] == 0){
    grid_idx = grid_idx[1:(length(grid_idx)-1)]
  }
  return(list(t=t, l=l, C=C, D=D, y=y, count=count, gridrep=gridrep,
              ns=sum(n_sampled), nc=nc, ng=ng, rep_idx=rep_idx,
              args=list(samp_times=samp_times, n_sampled=n_sampled, coal_times=coal_times,
                        grid=grid),
              grid_idx = grid_idx))
}




as.DatePhylo = function(phy, endTime, States){
  n = phy$Nnode + 1
  Sampling_Times = branching_sampling_times2(phy)
  Sampling_Times = endTime - Sampling_Times[n:(2 * n - 1)]
  names(Sampling_Times) = phy$tip.label
  return(DatedTree(phy,sampleTimes = Sampling_Times,sampleStates = States))
}





branching_sampling_times2 <- function(phy)
{
  phy = ape::new2old.phylo(phy)

 # if (class(phy) != "phylo")
 #   stop("object \"phy\" is not of class \"phylo\"")

  tmp <- as.numeric(phy$edge)
  nb.tip <- max(tmp)
  nb.node <- -min(tmp)
  xx <- as.numeric(rep(NA, nb.tip + nb.node))
  names(xx) <- as.character(c(-(1:nb.node), 1:nb.tip))
  xx["-1"] <- 0

  for (i in 2:length(xx))
  {
    nod <- names(xx[i])
    ind <- which(phy$edge[, 2] == nod)
    base <- phy$edge[ind, 1]
    xx[i] <- xx[base] + phy$edge.length[ind]
  }

  depth <- max(xx)
  branching_sampling_times <- depth - xx

  return(branching_sampling_times)
}


colikcpp_prepare <- function (bdt, times,DEMES, timeOfOriginBoundaryCondition = TRUE,
                       AgtYboundaryCondition = TRUE, maxHeight = Inf)
{

  bdt <- reorder.phylo(bdt, "postorder")
  bdt$heights <- signif(bdt$heights, digits = floor(1/bdt$maxHeight/10) +
                          6)
  #times <- tfgy[[1]]
  #Fs <- tfgy[[2]]
  #Gs <- tfgy[[3]]
  #Ys <- tfgy[[4]]
  m <- length(DEMES)
  if (m < 2)
    stop("Currently only models with at least two demes are supported")

  if (m == 2 & DEMES[2] == "V2" & ncol(bdt$sampleStates) ==
      1) {
    bdt$sampleStates <- cbind(bdt$sampleStates, rep(0, bdt$n))
    bdt$sortedSampleStates <- cbind(bdt$sortedSampleStates,
                                    rep(0, bdt$n))
  }
  fgyi <- 1:length(times)
  if (times[2] > times[1])
    fgyi <- length(times):1

  t0 <- times[fgyi[length(fgyi)]]
  if (timeOfOriginBoundaryCondition) {
    if ((bdt$maxSampleTime - t0) < bdt$maxHeight)
      return(-Inf)
  }

  heights <- times[fgyi[1]] - times[fgyi]
  ndesc <- rep(0, bdt$n + bdt$Nnode)
  for (iedge in 1:nrow(bdt$edge)) {
    a <- bdt$edge[iedge, 1]
    u <- bdt$edge[iedge, 2]
    ndesc[a] <- ndesc[a] + ndesc[u] + 1
  }
  uniqhgts <- sort(unique(bdt$heights))
  eventHeights <- rep(NA, (bdt$n + bdt$Nnode))
  eventIndicatorNode <- rep(NA, bdt$n + bdt$Nnode)
  events <- rep(NA, bdt$n + bdt$Nnode)
  hgts2node <- lapply(uniqhgts, function(h) which(bdt$heights ==
                                                    h))
  k <- 1
  l <- 1
  for (h in uniqhgts) {
    us <- hgts2node[[k]]
    if (k < length(hgts2node) | length(us) == 1) {
      if (length(us) > 1) {
        i_us <- sort(index.return = TRUE, decreasing = FALSE,
                     ndesc[us])$ix
        for (u in us[i_us]) {
          eventHeights[l] <- h
          events[l] <- ifelse(u <= bdt$n, 0, 1)
          eventIndicatorNode[l] <- u
          l <- l + 1
        }
      }
      else {
        eventHeights[l] <- h
        events[l] <- ifelse(us <= bdt$n, 0, 1)
        eventIndicatorNode[l] <- us
        l <- l + 1
      }
    }
    k <- k + 1
  }
  excl <- is.na(eventHeights) | is.na(events) | is.na(eventIndicatorNode) |
    (eventHeights > maxHeight)
  events <- events[!excl]
  eventIndicatorNode <- eventIndicatorNode[!excl]
  eventHeights <- eventHeights[!excl]
 # ll <- colik2cpp(heights, Fs[fgyi], Gs[fgyi], Ys[fgyi],
  #                events, eventIndicatorNode, eventHeights, t(bdt$sortedSampleStates),
   #               bdt$daughters, bdt$n, bdt$Nnode, m, AgtYboundaryCondition)
  return(list(heights = heights, events = events, eventIndicatorNode = eventIndicatorNode,
              eventHeights = eventHeights, sortedSampleStates = t(bdt$sortedSampleStates),
              daughters = bdt$daughters, n = bdt$n, Nnode = bdt$Nnode, m = m))
}


colik.fgy <- function (bdt, tfgy, timeOfOriginBoundaryCondition = TRUE,
                       AgtYboundaryCondition = TRUE, maxHeight = Inf, expmat = FALSE)
{

  bdt <- reorder.phylo(bdt, "postorder")
  bdt$heights <- signif(bdt$heights, digits = floor(1/bdt$maxHeight/10) +
                          6)
  times <- tfgy[[1]]
  Fs <- tfgy[[2]]
  Gs <- tfgy[[3]]
  Ys <- tfgy[[4]]
  m <- nrow(Fs[[1]])
  if (m < 2)
    stop("Currently only models with at least two demes are supported")
  DEMES <- names(Ys[[1]])
  if (is.null(DEMES)){
    DEMES <- rownames(Fs[[1]])
  }
  if (is.null(DEMES)){
    DEMES <- rownames(Gs[[1]])
  }
  if (is.null(DEMES)){
    stop('demographic model returns trajectory without deme names')
  }
  if (m == 2 & DEMES[2] == "V2" & ncol(bdt$sampleStates) ==
      1) {
    bdt$sampleStates <- cbind(bdt$sampleStates, rep(0, bdt$n))
    bdt$sortedSampleStates <- cbind(bdt$sortedSampleStates,
                                    rep(0, bdt$n))
  }
  fgyi <- 1:length(Fs)
  if (times[2] > times[1])
    fgyi <- length(Fs):1

  t0 <- times[fgyi[length(fgyi)]]
  if (timeOfOriginBoundaryCondition) {
    if ((bdt$maxSampleTime - t0) < bdt$maxHeight)
      return(-Inf)
  }

  heights <- times[fgyi[1]] - times[fgyi]
  ndesc <- rep(0, bdt$n + bdt$Nnode)
  for (iedge in 1:nrow(bdt$edge)) {
    a <- bdt$edge[iedge, 1]
    u <- bdt$edge[iedge, 2]
    ndesc[a] <- ndesc[a] + ndesc[u] + 1
  }
  uniqhgts <- sort(unique(bdt$heights))
  eventHeights <- rep(NA, (bdt$n + bdt$Nnode))
  eventIndicatorNode <- rep(NA, bdt$n + bdt$Nnode)
  events <- rep(NA, bdt$n + bdt$Nnode)
  hgts2node <- lapply(uniqhgts, function(h) which(bdt$heights ==
                                                    h))
  k <- 1
  l <- 1
  for (h in uniqhgts) {
    us <- hgts2node[[k]]
    if (k < length(hgts2node) | length(us) == 1) {
      if (length(us) > 1) {
        i_us <- sort(index.return = TRUE, decreasing = FALSE,
                     ndesc[us])$ix
        for (u in us[i_us]) {
          eventHeights[l] <- h
          events[l] <- ifelse(u <= bdt$n, 0, 1)
          eventIndicatorNode[l] <- u
          l <- l + 1
        }
      }
      else {
        eventHeights[l] <- h
        events[l] <- ifelse(us <= bdt$n, 0, 1)
        eventIndicatorNode[l] <- us
        l <- l + 1
      }
    }
    k <- k + 1
  }
  excl <- is.na(eventHeights) | is.na(events) | is.na(eventIndicatorNode) |
    (eventHeights > maxHeight)
  events <- events[!excl]
  eventIndicatorNode <- eventIndicatorNode[!excl]
  eventHeights <- eventHeights[!excl]
  ll <- colik2cpp(heights, Fs[fgyi], Gs[fgyi], Ys[fgyi],
                  events, eventIndicatorNode, eventHeights, t(bdt$sortedSampleStates),
                  bdt$daughters, bdt$n, bdt$Nnode, m, AgtYboundaryCondition)

  return(ll)
}




TruncTree_Fun = function(bigTree, date){
  cutID = sapply(bigTree$tip.label, function(x){
    a = strsplit(x,"[|\']")[[1]][7]
    return(as.Date(a) < as.Date(date))
  })
  return(drop.tip(bigTree, bigTree$tip.label[!cutID],subtree = F))
}


forwardSimulateODE_RST_core = function(param, LatentTraj, x_r, x_i, p, dt, Ngrid, beginning = F, tol=5){
  R_tstart = param[p + x_i[3] + 1] * prod(param[(p + x_i[2] + 1):(p+x_i[2]+x_i[1])])
  if(R_tstart > 0.999){
    return(Inf)
  }

  if(is.matrix(LatentTraj)){
    X1 = tail.matrix(LatentTraj,n=1)
  }else{
  X1 = LatentTraj
  }

  if(beginning == T){
    t = seq(LatentTraj[1,1], X1[1] + dt * Ngrid, by=dt)
    return(ODE_rk45_stop(X1[-1], t, param[-(1:p)], x_r, x_i, tol = tol));
  }else{
    t = seq(X1[1],X1[1] + dt * Ngrid, by=dt)
    return(ODE_rk45_stop(X1[-1], t, param[-(1:p)], x_r, x_i, tol = tol));
  }
}


############
SigmaFirst2 = function(n,mu = 0.0001){
  return(exp(mean(log(c(1:(n))))))
}

tauDist = function(a,b,dt){
  return(qgamma(a,b))
}


SEIR_simu_tree = function(tfgys, sampleTimes, sampleStates, tree = T){
  Tree = sim.co.tree.fgy(tfgys, sampleTimes, sampleStates, substitutionRates = NULL, sequenceLength = 0)
  class(Tree) = "phylo"
  if(tree == T){
    return(Tree)
  }else{
    return(summarize_phylo(Tree))
  }
}
