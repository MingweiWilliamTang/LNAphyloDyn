#Rcpp::sourceCpp('src/SIRS_period.cpp')
#Rcpp::sourceCpp("src/SIR_BD.cpp")
#source('~/Dropbox/Research/GP/code/LNAPhyloDyn/R/coal_simulation.R', echo=TRUE)
#source('~/Dropbox/Research/GP/code/LNAPhyloDyn/R/SIR_LNA_standard.R', echo=TRUE)

plot3 = function(traj){
  dev.off()
  par(mfrow = c(2,2),mar=c(4,4,1,1))
  plot(traj[,1],traj[,2],type = "l", xlab = "time", ylab = "suscepitible1")
  plot(traj[,1],traj[,3],type = "l", xlab = "time", ylab = "Infected")
  plot(traj[,1],traj[,4],type = "l", xlab = "time", ylab = "Recovered")
  #plot(traj[,1],traj[,5],type = "l", xlab = "time", ylab = "reconvered")
}



effpopfun = function(Traj,beta=0,lambda=1, volz = FALSE){
  if(volz){
    return(1 /(2 * Traj[,3] * beta / Traj[,2]))
  }else{
    return(Traj[,3] / lambda)
  }

}

#' Plot the trajectory of effective population size
#'
#' @param MCMC_res returned MCMC list
#' @param idxs The index choosed to make plot
#' @param volz If use volz likelihood the coalescent modeling
#'
#'
#'
random_effpoptraj_line = function(MCMC_res,idxs,volz = FALSE){
  N = MCMC_res$Trajectory[1,2,1] + MCMC_res$Trajectory[1,3,1]
  #plot(N)
  if(volz){
    for(i in idxs){
      lines(MCMC_res$Trajectory[,1,i],
            effpopfun(MCMC_res$Trajectory[,,i],beta = MCMC_res$par[i,3] * N,volz = T),
            lwd = 0.5, col = "grey")
    }
  }else{
    for(i in idxs){
      lines(MCMC_res$Trajectory[,1,i],
            MCMC_res$Trajectory[,3,i] / MCMC_res$par[i,5],
            lwd = 0.5, col = "yellow")
    }
  }
}



medianCur = function(MCMC_obj,ids,scale=1,col="red",row=3,med = T,volz = F,lwd = 2,lid=7){
  if(volz == F){
    scale = MCMC_obj$par[ids,lid]
  }
  if(med){
  lines(MCMC_obj$Trajectory[,1,1],apply(MCMC_obj$Trajectory[,row,ids]/scale,1,median),lty=2,
          col=col,lwd=lwd)
  }else{
  lines(MCMC_obj$Trajectory[,1,1],apply(MCMC_obj$Trajectory[,row,ids]/scale,1,mean),lty=2,
        col=col,lwd=lwd)
  }
}


effpopfun = function(Traj,beta=0,lambda=1, volz = FALSE){
  if(volz){
    return(1/(2 * Traj[,2] * beta / Traj[,3]))
  }else{
    return(Traj[,3] / lambda)
  }
}




CI_Curve = function(MCMC_obj,ids,scale = 1, col = "black", fill_col = "grey", row = 3,method = "qtile",alpha=0.05,fill = T){

  if(method == "NormApp"){
    midcur = apply(MCMC_obj$Trajectory[,row,ids]/scale,1,mean)
    midSd = apply(MCMC_obj$Trajectory[,row,ids]/scale,1,sd)
    qt = qnorm(1-alpha/2)
    upper = midcur + qt * midSd
    lower = midcur - qt * midSd
    #lines(MCMC_obj$Trajectory[,1,1],upper,lty=2,
    #      col=col,lwd=2)

    #lines(MCMC_obj$Trajectory[,1,1],lower,lty=2,
    #      col=col,lwd=2)
  }else if(method == "qtile"){
    qt1 = 1 - alpha/2
    qt2 = alpha/2
    upper = apply(MCMC_obj$Trajectory[,row,ids]/scale,1,function(x){
      return(quantile(x,qt1))
    }
      )
    lower = apply(MCMC_obj$Trajectory[,row,ids]/scale,1,function(x){
      return(quantile(x,qt2))
    })
   # lines(MCMC_obj$Trajectory[,1,1],midcur + qt * midSd,lty=2,
    #      col=col,lwd=2)

    #lines(MCMC_obj$Trajectory[,1,1],midcur - qt * midSd,lty=2,
     #     col=col,lwd=2)
  }
  if(fill == T){
    polygon(x = c(MCMC_obj$Trajectory[,1,1],rev(MCMC_obj$Trajectory[,1,1])),
            y = c(upper,rev(lower)),col = fill_col,border = NA)
  }
}


CI_Curve_eff = function(MCMC_obj,ids,scale = 1, col = "black", fill_col = "grey", row = 3,method = "qtile",alpha=0.05,fill = T,likelihood ="volz"){
  if(likelihood == "volz"){
    Mx = 1/(2 * MCMC_obj$Trajectory[,2,ids]/MCMC_obj$Trajectory[,3,ids])
  }else{
    Mx = MCMC_obj$Trajectory[,row,ids]
  }
  midcur = apply(Mx/scale,1,mean)
  if(method == "NormApp"){
    midSd = apply(Mx/scale,1,sd)
    qt = qnorm(1-alpha/2)
    upper = midcur + qt * midSd
    lower = midcur - qt * midSd
    #lines(MCMC_obj$Trajectory[,1,1],upper,lty=2,
    #      col=col,lwd=2)

    #lines(MCMC_obj$Trajectory[,1,1],lower,lty=2,
    #      col=col,lwd=2)
  }else if(method == "qtile"){
    qt1 = 1 - alpha/2
    qt2 = alpha/2
    upper = apply(Mx/scale,1,function(x){
      return(quantile(x,qt1))
    }
    )
    lower = apply(Mx/scale,1,function(x){
      return(quantile(x,qt2))
    })
    # lines(MCMC_obj$Trajectory[,1,1],midcur + qt * midSd,lty=2,
    #      col=col,lwd=2)

    #lines(MCMC_obj$Trajectory[,1,1],midcur - qt * midSd,lty=2,
    #     col=col,lwd=2)
  }
  if(fill == T){
    polygon(x = c(MCMC_obj$Trajectory[,1,1],rev(MCMC_obj$Trajectory[,1,1])),
            y = c(upper,rev(lower)),col = fill_col,border = NA)
  }
  lines(MCMC_obj$Trajectory[,1,1],midcur,col = col,lwd=2,lty=2)
}




vlineCI = function(data,cred = T){
  if(cred){
    m = median(data)
    s1 = quantile(data,0.025)
    s2 = quantile(data,1 - 0.025)
    abline(v = s1,col="blue",lwd=2,lty=2)
    abline(v = s2,col="blue",lwd=2,lty=2)
  }
  else{
  s = sd(data)
  m = meam(data)
  abline(v = m + 1.96*s,col="blue",lwd=2,lty=2)
  abline(v = m - 1.96*s,col="blue",lwd=2,lty=2)
  }
  abline(v = m ,col="blue",lwd=2)
}

randomR0_traj = function(times,MCMC_obj,R0_id,col_id,idx,ylim=c(0,2),main = ""){
  R0 = MCMC_obj$par[idx,R0_id]
  R0_traj = matrix(ncol= length(col_id)+2, nrow = length(idx))
  R0_traj[,1] = R0
  for(i in 1:length(col_id)){
    R0_traj[,i+1] = R0_traj[,i] * MCMC_obj$par[idx,col_id[i]]
  }
  i = length(col_id)
  R0_traj[,length(col_id)+2] = R0_traj[,i] * MCMC_obj$par[idx,col_id[i]]
  CIup = apply(R0_traj,2,function(x){
    return(quantile(x,0.975))
  })
  CIlow = apply(R0_traj,2,function(x){
    return(quantile(x,0.025))
  })
  m = apply(R0_traj,2,median)
  plot(times,m,type="l",ylab = "R0",col = "red",lwd = 2,ylim=ylim,main=main)
  polygon(x = c(times,rev(times)),
          y = c(CIup,rev(CIlow)),col = "grey",border = NA)
  lines(times,m,type="l",col = "red",lwd = 2)
}


randomR0_traj_V = function(times,MCMC_obj,R0_id,col_id,idx,xlim,ylim=c(0,2),main = ""){
  R0 = MCMC_obj$par[idx,R0_id]
  R0_traj = matrix(ncol= length(col_id)+1, nrow = length(idx))
  R0_traj[,1] = R0
  for(i in 1:length(col_id)){
    R0_traj[,i+1] = R0_traj[,i] * MCMC_obj$par[idx,col_id[i]]
  }
  i = length(col_id)
  #R0_traj[,length(col_id)+2] = R0_traj[,i] * MCMC_obj$par[idx,col_id[i]]
  CIup = apply(R0_traj,2,function(x){
    return(quantile(x,0.975))
    })
  CIlow = apply(R0_traj,2,function(x){
      return(quantile(x,0.025))
    })
  m = apply(R0_traj,2,median)
  times = times[2:(length(m))]
  step1 = stepfun(times, m, f = 0)
  step2 = stepfun(times, CIup, f = 0)
  step3 = stepfun(times, CIlow, f = 0)
  plot(step1,ylab = "R0",col = "red",lwd = 2.5,ylim=ylim,
       main=main,verticals = F,xlim=xlim,xlab = "time")
  #polygon(x = c(times,rev(times)),
   #       y = c(CIup,rev(CIlow)),col = "grey",border = NA)
  lines(step2, lty=2,lwd = 1,verticals = F, col = "blue",xlim=xlim)
  lines(step3, lty=2,lwd = 1,verticals = F, col = "blue",xlim=xlim)
  #lines(times,m,type="l",col = "red",lwd = 2)
}


coal_like_fast = function(traj, lambda = 1, coal_obj,t_correct,col=3){
  init = coal_lik_init(coal_obj$samp_times,coal_obj$n_sampled,coal_obj$coal_times,grid = traj[,1])
  return(coal_loglik(init,LogTraj(traj),t_correct = t_correct, lambda))
}

histwithCI = function(data,cred = T){
  hist(data);
  vlineCI(data,cred)
}

effpopPlot = function(NP_res,LNA_res,t_correct,idx,row=id,lid=5,ylab="",xlab="",ylim=c(0,10)){
  plot(t_correct - NP_res$x,NP_res$effpopmean,type="l",ylim = ylim,xlab = xlab,ylab=ylab)
  polygon(c(t_correct - NP_res$x, rev(t_correct - NP_res$x)), c(NP_res$effpop975,rev(NP_res$effpop025)),
          col=rgb(0,0,1,0.3),border = F)
  CI_Curve(LNA_res,idx,scale = LNA_res$par[idx,5],fill_col = rgb(1,0,0,0.3))
  medianCur(LNA_res,idx,col="red",lid = 5)
}




CI_Curve_eff2 = function(MCMC_obj,ids, col = "black", fill_col = "grey", Irow = 3,method = "qtile",alpha=0.05,fill = T,likelihood ="volz",
                         x_r,x_i,p=3){
  if(likelihood == "volz"){
    Mx = 1/(2 * MCMC_obj$Trajectory[,2,ids]/MCMC_obj$Trajectory[,Irow,ids])
    scale = sapply(ids,function(x){
      return(betaTs(MCMC_obj$par[x,(p+1):(p + x_i[1] + x_i[2])],
                   MCMC_obj$Trajectory[,1,1],x_r, x_i))})
    Mx = Mx/scale
    }else{
    Mx = MCMC_obj$Trajectory[,Irow,ids]
    scale = 1
  }
  midcur = apply(Mx,1,mean)
  if(method == "NormApp"){
    midSd = apply(Mx,1,sd)
    qt = qnorm(1-alpha/2)
    upper = midcur + qt * midSd
    lower = midcur - qt * midSd
    #lines(MCMC_obj$Trajectory[,1,1],upper,lty=2,
    #      col=col,lwd=2)

    #lines(MCMC_obj$Trajectory[,1,1],lower,lty=2,
    #      col=col,lwd=2)
  }else if(method == "qtile"){
    qt1 = 1 - alpha/2
    qt2 = alpha/2
    upper = apply(Mx,1,function(x){
      return(quantile(x,qt1))
    }
    )
    lower = apply(Mx,1,function(x){
      return(quantile(x,qt2))
    })
    # lines(MCMC_obj$Trajectory[,1,1],midcur + qt * midSd,lty=2,
    #      col=col,lwd=2)

    #lines(MCMC_obj$Trajectory[,1,1],midcur - qt * midSd,lty=2,
    #     col=col,lwd=2)
  }
  if(fill == T){
    polygon(x = c(MCMC_obj$Trajectory[,1,1],rev(MCMC_obj$Trajectory[,1,1])),
            y = c(upper,rev(lower)),col = fill_col,border = NA)
  }
  lines(MCMC_obj$Trajectory[,1,1],midcur,col = col,lwd=2,lty=2)
}
