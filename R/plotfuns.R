#Rcpp::sourceCpp('src/SIRS_period.cpp')
#Rcpp::sourceCpp("src/SIR_BD.cpp")
#source('~/Dropbox/Research/GP/code/LNAPhyloDyn/R/coal_simulation.R', echo=TRUE)
#source('~/Dropbox/Research/GP/code/LNAPhyloDyn/R/SIR_LNA_standard.R', echo=TRUE)

random_trajectory = function(t,gridsize,trMatrix,color="grey",ylab="infected population",lwd=0.5,margin=500){
  timegrid = seq(0,t,by=gridsize)
  k = dim(trMatrix)[2]
  plot(timegrid,trMatrix[,1],lwd=lwd,col=color,ylab=ylab,
       ylim=c(min(trMatrix[,1])-margin,max(trMatrix[,1])+margin),
       main="infectious population vs time",type="l")
  for(i in 2:k){
    lines(timegrid,trMatrix[,i],lwd=lwd,col=color)
  }
}


random_trajectory_line = function(t,gridsize,trMatrix,color="yellow",lwd=0.5){
  timegrid = seq(0,t,by=gridsize)
  k = dim(trMatrix)[2]
  lines(timegrid,trMatrix[,1],lwd=lwd,col=color,ylim=c(min(trMatrix[,1])-500,max(trMatrix[,1])+500),
        main="infectious population vs time",type="l")
  for(i in 2:k){
    lines(timegrid,trMatrix[,i],lwd=lwd,col=color)
  }
}
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
    return(1 /(2 * Traj[,2] * beta / Traj[,3]))
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
  N = MCMC_res$Trajectory[,2,1] + MCMC_res$Trajectory[,3,1]
  plot(N)
  if(volz){
    for(i in idxs){
      lines(MCMC_res$Trajectory[,1,i],
            effpopfun(MCMC_res$Trajectory[,,i],beta = MCMC_res$par[i,3] * N,volz = T),
            lwd = 0.5, col = "yellow")
    }
  }else{
    for(i in idxs){
      lines(MCMC_res$Trajectory[,1,i],
            MCMC_res$Trajectory[,3,i] / MCMC_res$par[i,5],
            lwd = 0.5, col = "yellow")
    }
  }
}



medianCur = function(MCMC_obj,ids,scale=1,col="red",row=3){
  lines(MCMC_obj$Trajectory[,1,1],apply(MCMC_obj$Trajectory[,row,ids]/scale,1,mean),lty=2,
        col=col,lwd=2)
}






CI_Curve = function(MCMC_obj,ids,scale = 1, col = "black", fill_col = "grey", row = 3,method = "NormApp",alpha=0.05,fill = T){

  if(method == "NormApp"){
    midcur = apply(MCMC_obj$Trajectory[,row,ids]/scale,1,mean)
    midSd = apply(MCMC_obj$Trajectory[,row,ids]/scale,1,sd)
    qt = qnorm(1-alpha/2)
    upper = midcur + qt * midSd
    lower = midcur - qt * midSd
    lines(MCMC_obj$Trajectory[,1,1],upper,lty=2,
          col=col,lwd=2)

    lines(MCMC_obj$Trajectory[,1,1],lower,lty=2,
          col=col,lwd=2)
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
    lines(MCMC_obj$Trajectory[,1,1],midcur + qt * midSd,lty=2,
          col=col,lwd=2)

    lines(MCMC_obj$Trajectory[,1,1],midcur - qt * midSd,lty=2,
          col=col,lwd=2)
  }
  if(fill == T){
    polygon(x = c(MCMC_obj$Trajectory[,1,1],rev(MCMC_obj$Trajectory[,1,1])),
            y = c(upper,rev(lower)),col = fill_col)
  }
}
