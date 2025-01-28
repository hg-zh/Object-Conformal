library(parallel)
Get_R<-function(r){
  set.seed(r)
  #library(rje)
  library(gridExtra)
  #####
  #source('~/Documents/Hans/object conformal/simulation_XY_R/utils.R')
  source('utils.R')
  Coverage<-list();Length<-list();Conf_Set<-list()
  for (setting in c(1,2,3)) {
    for (n in c(500,1000,2000)) {
      #### setting 1: unimodal ####
      if(setting==1){
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
        testing_x<-rep(work_grid_x,1000)
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
      if(setting==2){
        work_grid_x<-seq(-1,1,0.02)
        workgrid_p<-seq(0,7,0.02)
        Training_x<-list();Training_y<-list()
        for (i in 1:n) {
          Training_x[[i]]<-runif(1,-1,1)
          if(Training_x[[i]]<=0){
            Training_y[[i]]<-f(Training_x[[i]])+rnorm(1,0,0.5)
          }else{
            Training_y[[i]]<-f(Training_x[[i]])+rnorm(1,0,0.1)
          }
        }
        Training_set<-list(x=Training_x,y=Training_y)
        Training_set_L<-list()
        for (i in 1:n) {
          Training_set_L[[i]]<-list(x=Training_x[[i]],y=Training_y[[i]])
        }
        Calibration_x<-list();Calibration_y<-list()
        for (i in 1:n) {
          Calibration_x[[i]]<-runif(1,-1,1)
          if(Calibration_x[[i]]<=0){
            Calibration_y[[i]]<-f(Calibration_x[[i]])+rnorm(1,0,0.5)
          }else{
            Calibration_y[[i]]<-f(Calibration_x[[i]])+rnorm(1,0,0.1)
          }
        }
        Calibration_set<-list(x=Calibration_x,y=Calibration_y)
        Calibration_L<-list()
        for (i in 1:n) {
          Calibration_L[[i]]<-list(x=Calibration_x[[i]],y=Calibration_y[[i]])
        }
        testing_x<-rep(work_grid_x,1000)
        Testing_L<-list(x=testing_x,y=f(testing_x)+(testing_x>=0)*rnorm(length(testing_x),0,0.1)
                        +(testing_x<0)*rnorm(length(testing_x),0,0.5))
        Testing_set<-list()
        for (i in 1:length(Testing_L$x)) {
          Testing_set[[i]]<-list(x=Testing_L$x[i],y=Testing_L$y[i])
        }
      }
      #### setting 3: aleo ####
      if(setting==3){
        work_grid_x<-seq(-1,1,0.02)
        workgrid_p<-seq(0,5,0.02)
        Training_x<-list();Training_y<-list()
        for (i in 1:n) {
          Training_x[[i]]<-runif(1,-1,1)
          a<-rbinom(1,1,0.5)
          if(Training_x[[i]]<=0){
            Training_y[[i]]<-f(Training_x[[i]])+rnorm(1,0,0.1)
          }else if(a==0){
            Training_y[[i]]<-f(Training_x[[i]])+g(Training_x[[i]])+rnorm(1,0,0.1)
          }else{
            Training_y[[i]]<-f(Training_x[[i]])-0.2*g(Training_x[[i]])+rnorm(1,0,0.1)
          }
        }
        Training_set<-list(x=Training_x,y=Training_y)
        Training_set_L<-list()
        for (i in 1:n) {
          Training_set_L[[i]]<-list(x=Training_x[[i]],y=Training_y[[i]])
        }
        Calibration_x<-list();Calibration_y<-list()
        for (i in 1:n) {
          Calibration_x[[i]]<-runif(1,-1,1)
          a<-rbinom(1,1,0.5)
          if(Calibration_x[[i]]<=0){
            Calibration_y[[i]]<-f(Calibration_x[[i]])+rnorm(1,0,0.1)
          }else if(a==0){
            Calibration_y[[i]]<-f(Calibration_x[[i]])+g(Calibration_x[[i]])+rnorm(1,0,0.1)
          }else{
            Calibration_y[[i]]<-f(Calibration_x[[i]])-0.2*g(Calibration_x[[i]])+rnorm(1,0,0.1)
          }
        }
        Calibration_set<-list(x=Calibration_x,y=Calibration_y)
        Calibration_L<-list()
        for (i in 1:n) {
          Calibration_L[[i]]<-list(x=Calibration_x[[i]],y=Calibration_y[[i]])
        }
        testing_x<-rep(work_grid_x,1000)
        a<-rbinom(length(testing_x),1,0.5)
        Testing_L<-list(x=testing_x,y=f(testing_x)+(testing_x<=0)*(rnorm(length(testing_x),0,0.1))
                        +(testing_x>0)*( (a-0.2*(1-a))*g(testing_x)+ rnorm(length(testing_x),0,0.1)) )
        Testing_set<-list()
        for (i in 1:length(Testing_L$x)) {
          Testing_set[[i]]<-list(x=Testing_L$x[i],y=Testing_L$y[i])
        }
      }
      #### search grid ####
      gap=0.005
      grid <- expand.grid(x = work_grid_x, y = seq(-1,3,gap))
      Search_grid<-list()
      for (i in 1:length(grid$x)) {
        Search_grid[[i]]<-list(x=grid$x[i],y=grid$y[i])
      }
      Coverage[[as.character(setting) ]][[as.character(n)]]<-list()
      Length[[as.character(setting)]][[as.character(n)]]<-list()
      Conf_Set[[as.character(setting)]][[as.character(n)]]<-list()
      ##### Proposed scores #####
      #####use cross validation to choose tuning parameter
      if(setting==1){
        if(n==500){bd_cand<-c(seq(0.06,0.15,length=20))}
        if(n==1000){bd_cand<-c(seq(0.03,0.12,length=20))}
        if(n==2000){bd_cand<-c(seq(0.02,0.10,length=20))}
      }
      if(setting==2){
        if(n==500){bd_cand<-c(seq(0.06,0.2,length=20))}
        if(n==1000){bd_cand<-c(seq(0.03,0.15,length=20))}
        if(n==2000){bd_cand<-c(seq(0.02,0.15,length=20))}
      }
      if(setting==3){
        if(n==500){bd_cand<-c(seq(0.06,0.15,length=20))}
        if(n==1000){bd_cand<-c(seq(0.03,0.15,length=20))}
        if(n==2000){bd_cand<-c(seq(0.02,0.15,length=20))}
      }
      for (i in 1:length(bd_cand)) {
        Conf_Set[[as.character(setting)]][[as.character(n)]][[i]]<-Get_conf_set(h=bd_cand[i],metric=met_R1,
                                                          R_depth=R_n,Search_grid=Search_grid,Training_set_L,
                            Training_set,workgrid_p,Calibration_L,Calibration_set)
        C_set_proposed<-Conf_Set[[as.character(setting)]][[as.character(n)]][[i]]
        Data_C_x<-c();Data_C_y<-c()
        for (j in 1:length(C_set_proposed)) {
          Data_C_x[j]<-C_set_proposed[[j]]$x;
          Data_C_y[j]<-C_set_proposed[[j]]$y
        }
        Cov_temp<-c();Len_tmp<-c()
        for (j in 1:length(work_grid_x)) {
          Ind<-which(abs(testing_x-work_grid_x[j])<=0.01)
          if(length(Ind)==0){
            Cov_temp[j]<-NA
            Len_tmp[j]<-NA
          }else{
            C_set<-Data_C_y[which(Data_C_x==work_grid_x[j])]
            Res<-belongs_to_set(Testing_L$y[Ind],C_set,gap =gap )
            Cov_temp[j]<-sum(Res)/length(Res)
            Len_tmp[j]<-(length(C_set)-1)*gap
          }
        }
        Coverage[[as.character(setting)]][[as.character(n)]][[i]]<-Cov_temp
        Length[[as.character(setting)]][[as.character(n)]][[i]]<-Len_tmp
      }
    }
  }
  return(list(Coverage=Coverage,Length=Length))
}
system.time({
  cl <- makeCluster(getOption("cl.cores", 100));
  res<-parLapply(cl, 1:200,Get_R)
  stopCluster(cl)
  rm(cl)
})

save( res, file = paste0("May21th_R1_proposed_geq.RData") ) 


