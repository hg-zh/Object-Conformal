library(parallel)
Get_R<-function(r){
  set.seed(r)
  #library(rje)
  library(quantreg)
  library(gridExtra)
  #####
  #source('~/Documents/Hans/object conformal/simulation_XY_R/utils.R')
  source('utlis.R')
  alpha.sig<-0.1
  Coverage<-list();Length<-list();
  for (setting in c(1,2,3)) {
    for (n in c(500,1000,2000)) {
      #### setting 1: unimodal ####
      if(setting==1){
        work_grid_x<-seq(-1,1,0.02)
        workgrid_p<-seq(0,7,0.03)
        Training_x<-list();Training_y<-list()
        for (i in 1:n) {
          Training_x[[i]]<-runif(1,-1,1)
          Training_y[[i]]<-f(Training_x[[i]])+rnorm(1,0,1)
        }
        Training_set<-list(x=Training_x,y=Training_y)
        Training_set_L<-list()
        for (i in 1:n) {
          Training_set_L[[i]]<-list(x=Training_x[[i]],y=Training_y[[i]])
        }
        Calibration_x<-list();Calibration_y<-list()
        for (i in 1:n) {
          Calibration_x[[i]]<-runif(1,-1,1)
          Calibration_y[[i]]<-f(Calibration_x[[i]])+rnorm(1,0,1)
        }
        Calibration_set<-list(x=Calibration_x,y=Calibration_y)
        Calibration_L<-list()
        for (i in 1:n) {
          Calibration_L[[i]]<-list(x=Calibration_x[[i]],y=Calibration_y[[i]])
        }
        testing_x<-runif(101000,-1,1)
        # testing_y<-c()
        # for (i in 1:length(testing_x)) {
        #   x_tmp<-runif(1,max(-1,testing_x[i]-0.05),min(1,testing_x[i]+0.05))
        #   testing_y[i]<-f(x_tmp)+rnorm(1,0,0.1)
        # }
        Testing_L<-list(x=testing_x,y=f(testing_x)+rnorm(length(testing_x),0,1 ) )
        Testing_set<-list()
        for (i in 1:length(Testing_L$x)) {
          Testing_set[[i]]<-list(x=Testing_L$x[i],y=Testing_L$y[i])
        }
      }
      if(setting==2){
        work_grid_x<-seq(-1,1,0.02)
        workgrid_p<-seq(0,7,0.03)
        Training_x<-list();Training_y<-list()
        for (i in 1:n) {
          Training_x[[i]]<-runif(1,-1,1)
          Training_y[[i]]<-f(Training_x[[i]])+rnorm(1,0,0.5)
        }
        Training_set<-list(x=Training_x,y=Training_y)
        Training_set_L<-list()
        for (i in 1:n) {
          Training_set_L[[i]]<-list(x=Training_x[[i]],y=Training_y[[i]])
        }
        Calibration_x<-list();Calibration_y<-list()
        for (i in 1:n) {
          Calibration_x[[i]]<-runif(1,-1,1)
          Calibration_y[[i]]<-f(Calibration_x[[i]])+rnorm(1,0,0.5)
        }
        Calibration_set<-list(x=Calibration_x,y=Calibration_y)
        Calibration_L<-list()
        for (i in 1:n) {
          Calibration_L[[i]]<-list(x=Calibration_x[[i]],y=Calibration_y[[i]])
        }
        testing_x<-runif(101000,-1,1)
        # testing_y<-c()
        # for (i in 1:length(testing_x)) {
        #   x_tmp<-runif(1,max(-1,testing_x[i]-0.05),min(1,testing_x[i]+0.05))
        #   testing_y[i]<-f(x_tmp)+rnorm(1,0,0.1)
        # }
        Testing_L<-list(x=testing_x,y=f(testing_x)+rnorm(length(testing_x),0,0.5 ) )
        Testing_set<-list()
        for (i in 1:length(Testing_L$x)) {
          Testing_set[[i]]<-list(x=Testing_L$x[i],y=Testing_L$y[i])
        }
      }
      if(setting==3){
        work_grid_x<-seq(-1,1,0.02)
        workgrid_p<-seq(0,7,0.03)
        Training_x<-list();Training_y<-list()
        for (i in 1:n) {
          Training_x[[i]]<-runif(1,-1,1)
          Training_y[[i]]<-f(Training_x[[i]])+rnorm(1,0,0.1)
        }
        Training_set<-list(x=Training_x,y=Training_y)
        Training_set_L<-list()
        for (i in 1:n) {
          Training_set_L[[i]]<-list(x=Training_x[[i]],y=Training_y[[i]])
        }
        Calibration_x<-list();Calibration_y<-list()
        for (i in 1:n) {
          Calibration_x[[i]]<-runif(1,-1,1)
          Calibration_y[[i]]<-f(Calibration_x[[i]])+rnorm(1,0,0.1)
        }
        Calibration_set<-list(x=Calibration_x,y=Calibration_y)
        Calibration_L<-list()
        for (i in 1:n) {
          Calibration_L[[i]]<-list(x=Calibration_x[[i]],y=Calibration_y[[i]])
        }
        testing_x<-runif(101000,-1,1)
        # testing_y<-c()
        # for (i in 1:length(testing_x)) {
        #   x_tmp<-runif(1,max(-1,testing_x[i]-0.05),min(1,testing_x[i]+0.05))
        #   testing_y[i]<-f(x_tmp)+rnorm(1,0,0.1)
        # }
        Testing_L<-list(x=testing_x,y=f(testing_x)+rnorm(length(testing_x),0,0.1 ) )
        Testing_set<-list()
        for (i in 1:length(Testing_L$x)) {
          Testing_set[[i]]<-list(x=Testing_L$x[i],y=Testing_L$y[i])
        }
      }
      #### setting 2: heter case ####
      #### search grid ####
      ##### Proposed scores #####
      #####use cross validation to choose tuning parameter
      X0 <-as.matrix(unlist(Training_x)) 
      Y0 <- unlist(Training_y)
      X1 <-as.matrix( unlist(Calibration_x))
      Y1 <- unlist(Calibration_y)
      taus    <- seq(0.001,0.999,length=200)
      ys      <- quantile(unique(c(Y0,Y1)),seq(0.001,0.999,length=200))
      Y.test  <- unlist(Testing_L$y)
      X.test  <- as.matrix( unlist(Testing_L$x))
      res.qr    <- dcp.qr(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
      res.opt   <- dcp.opt(Y0,X0,Y1,X1,Y.test,X.test,taus,alpha.sig)
      res.dr    <- dcp.dr(Y0,X0,Y1,X1,Y.test,X.test,ys,taus,alpha.sig)
      ####
      Coverage[[as.character(setting) ]][[as.character(n)]]  <- binning(as.vector(X.test) ,cbind(res.qr$cov.qr,res.opt$cov.opt,res.dr$cov.dr),40)$cond
      Length[[as.character(setting)]][[as.character(n)]] <- binning(as.vector(X.test) ,cbind(res.qr$leng.qr,res.opt$leng.opt,res.dr$leng.dr),40)$cond
    }
  }
  return(list(Coverage=Coverage,Length=Length))
}
system.time({
  cl <- makeCluster(getOption("cl.cores", 50));
  res<-parLapply(cl, 1:200,Get_R)
  stopCluster(cl)
  rm(cl)
})

save( res, file = paste0("DCP_R1_new",format(Sys.time(),"%Y-%m-%d_%T"),".RData",sep="") ) 

