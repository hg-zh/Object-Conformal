library(parallel)
Get_R<-function(r){
  set.seed(r)
  library(rje)
  library(gridExtra)
  library(predictionBands)
  library(FlexCoDE)
  file_path <- "Output.txt"
  file_connection <- file(file_path, open = "a")
  writeLines(paste0("The run ",as.character(r)," is running"), file_connection)
  close(file_connection)
  #####
  source('utils.R')
  #source('/home/hangzh/Object_conformal/May6th/utils.R')
  Coverage<-list();Length<-list()
  for (setting in c(1,2,3)) {
    for (n in c(500,1000,2000,4000)) {
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
        testing_x<-rep(work_grid_x,each=100)
        Testing_L<-list(x=testing_x,y=f(testing_x)+rnorm(length(testing_x),0,0.1))
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
        testing_x<-rep(work_grid_x,each=100)
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
        a<-rbinom(length(testing_x),1,0.5)
        testing_x<-rep(work_grid_x,each=100)
        Testing_L<-list(x=testing_x,y=f(testing_x)+(testing_x<=0)*(rnorm(length(testing_x),0,0.1))
                        +(testing_x>0)*( (a-0.2*(1-a))*g(testing_x)+ rnorm(length(testing_x),0,0.1)) )
        Testing_set<-list()
        for (i in 1:length(Testing_L$x)) {
          Testing_set[[i]]<-list(x=Testing_L$x[i],y=Testing_L$y[i])
        }
      }
      
      
      #### search grid ####
      grid <- expand.grid(x = seq(-1, 1, length.out = 101), y = seq(-2,3, length.out = 101))
      Search_grid<-list()
      for (i in 1:length(grid$x)) {
        Search_grid[[i]]<-list(x=grid$x[i],y=grid$y[i])
      }
      Coverage[[as.character(setting)]][[as.character(n)]]<-matrix(0,3,ncol =length(work_grid_x) );
      Length[[as.character(setting)]][[as.character(n)]]<-matrix(0,3,ncol =length(work_grid_x) )
      ##### Proposed scores #####
      #####use cross validation to choose tuning parameter
      # bd_cand<-c(seq(0.05,0.2,0.01),0.25,0.3)
      # len_of_C<-c()
      # for (i in 1:length(bd_cand)) {
      #   C_set<-Get_conf_set(h=bd_cand[i],metric=met_R1,R_depth=R_n,Search_grid=Search_grid,Training_set_L,
      #                       Training_set,workgrid_p,Calibration_L,Calibration_set)
      #   len_of_C[i]<-length(C_set)
      # }
      # C_set_proposed<-Get_conf_set(h=bd_cand[which.min(len_of_C)],metric=met_R1,R_depth=R_n,Search_grid=Search_grid,
      #                              Training_set_L,Training_set,workgrid_p,Calibration_L,Calibration_set)
      # Data_C_x<-c();Data_C_y<-c()
      # for (i in 1:length(C_set_proposed)) {
      #   Data_C_x[i]<-C_set_proposed[[i]]$x;
      #   Data_C_y[i]<-C_set_proposed[[i]]$y
      # }
      # for (i in 1:length(work_grid_x)) {
      #   gap=0.05
      #   Ind<-which(testing_x==work_grid_x[i])
      #   C_set<-Data_C_y[which(Data_C_x==work_grid_x[i])]
      #   Res<-belongs_to_set(Testing_L$y[Ind],C_set,gap =gap )
      #   Coverage[[as.character(setting)]][[as.character(n)]][1,i]<-sum(Res)/length(Res)
      #   Length[[as.character(setting)]][[as.character(n)]][1,i]<-(length(C_set)-1)*gap
      # }
      ##### HPD #####
      fit<-fit_predictionBands(x=c(unlist(Training_x),unlist(Calibration_x)),
                               y=c(unlist(Training_y),unlist(Calibration_y)),0.5,0.4,0.1)
      xnew =matrix(work_grid_x,length(work_grid_x),1)
      bands<-predict(fit,xnew,type="hpd")
      plot(bands)
      for (i in 1:length(work_grid_x)) {
        Ind<-which(testing_x==work_grid_x[i])
        gap=bands$y_grid[2]-bands$y_grid[1]
        C_set<-bands$y_grid[bands$prediction_bands_which_belong[[i]]]
        Res<-belongs_to_set(Testing_L$y[Ind],C_set,gap =gap )
        Coverage[[as.character(setting)]][[as.character(n)]][2,i]<-sum(Res)/length(Res)
        Length[[as.character(setting)]][[as.character(n)]][2,i]<-(length(C_set)-1)*gap
      }
      
      ##### residual #####
      # Training_df<-data.frame(x=unlist(Training_set$x),y=unlist(Training_set$y))
      # loess_fit <- loess(y ~ x, data = Training_df,degree = 1)
      # Res_on_Cal<-abs(predict(loess_fit, newdata = data.frame(x = unlist(Calibration_x)))-unlist(Calibration_y) )
      # q<-quantile(Res_on_Cal,0.9*(1+1/length(Calibration_L)),names = FALSE,na.rm = T)
      # Pred_grid<-predict(loess_fit, newdata = data.frame(x = unlist(work_grid_x)))
      # Pred_grid[1]<-Pred_grid[2];Pred_grid[length(work_grid_x)]<-Pred_grid[length(work_grid_x)-1]
      # for (i in 1:length(work_grid_x)) {
      #   Ind<-which(testing_x==work_grid_x[i])
      #   gap=2*q
      #   C_set<-c(Pred_grid[i]-q,Pred_grid[i]+q)
      #   Res<-belongs_to_set(Testing_L$y[Ind],C_set,gap =gap )
      #   Coverage[[as.character(setting)]][[as.character(n)]][3,i]<-sum(Res)/length(Res)
      #   Length[[as.character(setting)]][[as.character(n)]][3,i]<-(length(C_set)-1)*gap
      # }
    }
  }
  save(Coverage,Length,file =paste0(as.character(r),".RData")  )
  file_connection <- file(file_path, open = "a")
  writeLines(paste0("The run ",as.character(r)," at n ",as.character(n)," is finished"), file_connection)
  close(file_connection)
  return(list(Coverage=Coverage,Length=Length))
}
system.time({
  cl <- makeCluster(getOption("cl.cores", 100));
  res<-parLapply(cl, 101:200,Get_R)
  stopCluster(cl)
  rm(cl)
})

save( res, file = paste0("May6th_HPD3.RData") ) 


