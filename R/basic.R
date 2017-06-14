loadModule("mod_Foo", TRUE)
#SIR_exact=list()
#SIR_exact$Pre=matrix(c(1,0,1,1,0,0),ncol=3)
#SIR_exact$Post = matrix(c(0,0,2,0,0,1),ncol=3)
#SIR_exact$h = function(x,t,th=c(theta1=0.00002,theta2=0.1)){
#  with(as.list(c(x,th)),{
#    return(c(theta1*X*Y, theta2*Y))
#  })
#}
#SIR_FRM = StepFRM(SIR_exact)

simuSIRS = function(theta1,theta2,theta3,S,I,R,times){
  SIRS_exact = list()
  SIRS_exact$Pre = matrix(c(1,1,0,0,1,0,0,0,1),ncol=3,byrow = T)
  SIRS_exact$Post = matrix(c(0,2,0,0,0,1,1,0,0),ncol=3,byrow = T)
  SIRS_exact$h = function(x,t,th=c(theta1=theta1,theta2=theta2,theta3 = theta3)){
    with(as.list(c(x,th)),{
      return(c(theta1*X*Y, theta2*Y,theta3 * Z))
    })
  }
  SIRS_FRM = StepFRM(SIRS_exact)
  simu_Traj = simTs(c(X = S, Y = I, Z = R), min(times), max(times), times[2] - times[1], SIRS_FRM)
  plot(times,simu_Traj[,2],type="l")
  return(cbind(times,simu_Traj))
}

#SIRS_F = simTs(c(X = 9500,Y = 500,Z = 1500),0,200,0.1,SIRS_FRM)
#plot(seq(0,200,by=0.1),SIRS_F[,2])

simuSIR = function(theta1,theta2,S,I,time){
  R = 0
  SIR = list()
  SIR$Pre=matrix(c(1,0,1,1,0,0),ncol=3)
  SIR$Post = matrix(c(0,0,2,0,0,1),ncol=3)
  SIR$h = function(x,t,th=c(theta1= theta1, theta2 = theta2)){
    with(as.list(c(x,th)),{
      return(c(theta1*X*Y, theta2*Y))
    })
  }
  SIR_F = StepFRM(SIR)
  simu_Traj = simTs(c(X = S, Y = I, Z = R), min(time), max(time), time[2] - time[1], SIR_F)
  plot(time,simu_Traj[,2],type="l")
  return(cbind(time,simu_Traj))
}



simuSEIR = function(theta1,theta2,theta3,S,E,I,time){
  R = 0
  SIR = list()
  SIR$Pre=matrix(c(1,0,1,0,1,0,0,0,1),byrow = T,ncol=3)
  SIR$Post = matrix(c(0,1,1,0,0,1,0,0,0),byrow = T,ncol=3)
  SIR$h = function(x,t,th=c(theta1= theta1, theta2 = theta2, theta3 = theta3)){
    with(as.list(c(x,th)),{
      return(c(theta1*X*Z, theta2*Y,theta3 * Z))
    })
  }
  SIR_F = StepFRM(SIR)
  simu_Traj = simTs(c(X = S, Y = E, Z = I), min(time), max(time), time[2] - time[1], SIR_F)
  #plot(time,simu_Traj[,2],type="l")
  return(cbind(time,simu_Traj))
}



#' Plot the posterior
#'
#' @param timegrid the time grid for the trajectory
#' @param trMatrix A matrix of trjactory for one population. The column corresponds to each posterior sample, the row corresponds to each time point
#' @param color The
#'
#'
random_trajectory = function(timegrid,trMatrix,color="grey",ylab="infected population",main="",lwd=0.5,ylim = NULL){
  k = dim(trMatrix)[2]
  if(is.null(ylim)){
  plot(timegrid,trMatrix[,1],lwd=lwd,col=color,ylab=ylab,
       ylim=c(min(trMatrix[,1],na.rm = T)-500,max(trMatrix[,1],na.rm = T)+500),
       main=main,type="l")
  }else{
    plot(timegrid,trMatrix[,1],lwd=lwd,col=color,ylab=ylab,
         ylim = ylim,
         main=main,type="l")
  }
  for(i in 2:k){
    lines(timegrid,trMatrix[,i],lwd=lwd,col=color)
  }
}



random_trajectory_line = function(timegrid,trMatrix,color="yellow",lwd=0.5, ylim = NULL){
 # timegrid = seq(t0,t,by=gridsize)
  k = dim(trMatrix)[2]
  if(is.null(ylim)){

    lines(timegrid,trMatrix[,1],lwd=lwd,col=color,
          ylim=c(min(trMatrix[,1])-500,max(trMatrix[,1])+500),
          main="population vs time",type="l")
  }else{
    lines(timegrid,trMatrix[,1],lwd=lwd,col=color,
          ylim=c(min(trMatrix[,1])-500,max(trMatrix[,1])+500),
          main="population vs time",type="l")

  }
  for(i in 2:k){
    lines(timegrid,trMatrix[,i],lwd=lwd,col=color)
  }
}


random_trajectory3D = function(t,gridsize,trCube){
  dev.off()
  par(mfrom(2,2))
  timegrid = seq(0,t,by=gridsize)
  k = dim(trCube)[3]
  picName = c("Susceptible","infectious","removed")
  for(j in 1:3){
    random_trajectory(t,gridsize,trCube[,i,])
    lines(timegrid,apply(trCube[,i,],1,mean),col="blue",lty=2,lwd=2)
    legend("topright",legend=c("samplepath",paste("mean ",picName[i])),col=c("grey","blue"))
  }
}

Convert_R0 = function(MCMC_obj,beta_row_id = 3, gamma_row_id=4,sample_id,S0){
  return(MCMC_obj$par[sample_id,beta_row_id]/MCMC_obj$par[sample_id,gamma_row_id]*S0)
}
