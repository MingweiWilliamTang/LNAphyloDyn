#' @useDynLib LNAPhyloDyn
#' @importFrom Rcpp sourceCpp
#'
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




random_trajectory = function(t,gridsize,trMatrix,color="grey",ylab="infected population",lwd=0.5){
  timegrid = seq(0,t,by=gridsize)
  k = dim(trMatrix)[2]
  plot(timegrid,trMatrix[,1],lwd=lwd,col=color,ylab=ylab,
       ylim=c(min(trMatrix[,1])-500,max(trMatrix[,1])+500),
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


random_trajectory3D = function(t,gridsize,trCube){
  dev.off()
  par(mfrom(2,2))
  timegrid = seq(0,t,by=gridsize)
  k = dim(trCube)[3]
  picName = c("Susceptible","infectious","removed")
  for(j in 1:3){
    random_trajectory(t,gridsize,trCube[,i,])
    lines(timegrid,apply(trCube[,i,],1,mean),col="blue",lty=2,lwd=2)
    length("topright",legend=c("samplepath",paste("mean ",picName[i])),col=c("grey","blue"))
  }
}
