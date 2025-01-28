library(parallel)
Get_R<-function(r){
  set.seed(r)
  #library(rje)
  library(gridExtra)
  #####
  #source('~/Documents/Hans/object conformal/simulation_XY_R/utils.R')
  source('utils.R')
  Coverage<-list();Length<-list();Conf_Set<-list();MSE_beta<-list();Beta<-list()
  beta0<-c(1,0,0,0)
  setting<-1
  for (n in c(200,500)) {
    #### setting 1: unimodal ####
    work_grid_x<-seq(-0.9,0.9,0.05)
    workgrid_p<-seq(0,3.5,length=101)
    Training_x<-list();Training_y<-list()
    for (i in 1:n) {
      Training_x[[i]]<-c(runif(1,-1,1),runif(1,-1,1),runif(1,-1,1),runif(1,-1,1))
      # mu<-c(1,0,0)
      # V<-c(0,rnorm(2,0,0.2))
      mu<-as.vector(f_sphere(as.numeric(Training_x[[i]]%*%beta0))) 
      V<-c(0,0,rnorm(1,0,0.5))
      Training_y[[i]]<-cos(norm(V,type = "2"))*mu+
        sin(norm(V,type = "2"))*V/norm(V,type = "2")
    }
    Training_set<-list(x=Training_x,y=Training_y)
    Training_set_L<-list()
    for (i in 1:n) {
      Training_set_L[[i]]<-list(x=Training_x[[i]],y=Training_y[[i]])
    }
    Calibration_x<-list();Calibration_y<-list()
    for (i in 1:n) {
      Calibration_x[[i]]<-c(runif(1,-1,1),runif(1,-1,1),runif(1,-1,1),runif(1,-1,1))
      # mu<-c(1,0,0)
      # V<-c(0,rnorm(2,0,0.2))
      mu<-as.vector(f_sphere(as.numeric(Calibration_x[[i]]%*%beta0))) 
      V<-c(0,0,rnorm(1,0,0.5))
      Calibration_y[[i]]<-cos(norm(V,type = "2"))*mu+
        sin(norm(V,type = "2"))*V/norm(V,type = "2")
    }
    Calibration_set<-list(x=Calibration_x,y=Calibration_y)
    Calibration_set_L<-list()
    for (i in 1:n) {
      Calibration_set_L[[i]]<-list(x=Calibration_x[[i]],y=Calibration_y[[i]])
    }
    testing_x<-runif(37000,-1,1);Testing_y<-list()
    for (i in 1:length(testing_x)) {
      
      mu<-as.vector(f_sphere(testing_x[i]))
      V<-c(0,0,rnorm(1,0,0.5))
      Testing_y[[i]]<-cos(norm(V,type = "2"))*mu+
        sin(norm(V,type = "2"))*V/norm(V,type = "2")
    }
    Testing_L<-list(x=testing_x,y=Testing_y)
    Testing_set<-list()
    for (i in 1:length(Testing_L$x)) {
      Testing_set[[i]]<-list(x=Testing_L$x[i],y=Testing_L$y[i])
    }
    #### solve beta ####
    if(n<=500){bd_can_beta<-c(0.1,0.15,0.2,0.25,0.3)}
    if(n==1000){bd_can_beta<-c(0.05,0.1,0.15,0.2,0.25)}
    if(n==2000){bd_can_beta<-c(0.05,0.1,0.15,0.2,0.25)}
    ID_loss<-function(x,Training_set,h=NULL,work_grid_x){
      theta<-c(cos(x[1]),sin(x[1])*cos(x[2]),sin(x[1])*sin(x[2])*cos(x[3]),sin(x[1])*sin(x[2])*sin(x[3]) )
      X_theta <- unlist(lapply(Training_set$x,  function(x) t(x)%*%theta )) 
      Y_theta<-ll_frechet_CV(X_theta,Training_set$y, 
                             h=NULL, X_theta, cv_folds = 2,bd_can_beta) 
      return( sum( mapply(function(x,y) met_sphere(x,y),Y_theta,Training_set$y )  ))
    }
    hat_theta<-optim(c(0,0,0),function(x) ID_loss(x,Training_set,h=NULL,work_grid_x),
                     lower=c(-pi,-pi,-pi),upper = c(pi,pi,pi),method="L-BFGS-B") 
    beta_hat<-c(cos(hat_theta$par[1]),sin(hat_theta$par[1])*cos(hat_theta$par[2]),
                sin(hat_theta$par[1])*sin(hat_theta$par[2])*cos(hat_theta$par[3]),
                sin(hat_theta$par[1])*sin(hat_theta$par[2])*sin(hat_theta$par[3]) )
    MSE_beta[[as.character(setting)]][[as.character(n)]]<-sqrt(sum( (beta_hat-beta0 )^2 ) )
    Beta[[as.character(setting)]][[as.character(n)]]<-beta_hat
    #### convert to 1d ####
    Training_x_h<-list();Calibration_x_h<-list()
    for (i in 1:n) {
      Training_x_h[[i]]<-as.numeric(Training_x[[i]]%*%beta_hat)
      Calibration_x_h[[i]]<-as.numeric(Calibration_x[[i]]%*%beta_hat)
    }
    Index_train<-which((Training_x_h)<=1&(Training_x_h)>=-1)
    Training_x_h<-Training_x_h[Index_train]
    Index_cal<-which((Calibration_x_h)<=1&(Calibration_x_h)>=-1)
    Calibration_x_h<-Calibration_x_h[Index_cal]
    Training_set_h<-list(x=Training_x_h,y=Training_y[Index_train])
    Calibration_set_h<-list(x=Calibration_x_h,y=Calibration_y[Index_cal])
    Training_set_L_h<-list();Calibration_set_L_h<-list()
    for (i in 1:length(Training_x_h)) {
      Training_set_L_h[[i]]<-list(x=Training_set_h$x[[i]] ,y=Training_set_h$y[[i]] )
    }
    for (i in 1:length(Calibration_x_h)) {
      Calibration_set_L_h[[i]]<-list(x=Calibration_set_h$x[[i]] ,y=Calibration_set_h$y[[i]] )
    }
    #### search grid ####
    u <- seq(0, pi, length.out = 100)
    v <- seq(0, 2 * pi, length.out = 100)
    x <- outer(cos(u), cos(v), FUN = "*")
    y <- outer(cos(u), sin(v), FUN = "*")
    z <- outer(sin(u), rep(1, length(v)), FUN = "*")
    search_grid <- data.frame(x = c(c(x),c(x)), y = c(c(y),c(y)) , z = c(c(z),-c(z)) )
    Search_grid <- vector("list", length(work_grid_x) * nrow(search_grid))
    idx <- 1
    # Nested lapply() functions to replace for-loops
    Search_grid <- do.call(c, lapply(1:length(work_grid_x), function(i) {
      lapply(1:nrow(search_grid), function(j) {
        temp <- list(x = work_grid_x[i], y =c(search_grid[j, 1],search_grid[j, 2],search_grid[j, 3]) )
        idx <<- idx + 1
        temp
      })
    }))
    Coverage[['1']][[as.character(n)]]<-list()
    Length[['1']][[as.character(n)]]<-list()
    Conf_Set[['1']][[as.character(n)]]<-list()
    ##### Proposed scores #####
    #####use cross validation to choose tuning parameter
    if(n<=500){bd_cand<-c(seq(0.1,0.25,length=10))}
    if(n==800){bd_cand<-c(seq(0.08,0.2,length=10))}
    if(n>800){bd_cand<-c(seq(0.04,0.1,length=10))}
    for (i in 1:length(bd_cand)) {
      Conf_Set[[as.character(setting)]][[as.character(n)]][[i]]<-Get_conf_set(h=bd_cand[i],metric=met_sphere,
                                                                              R_depth=R_n,Search_grid=Search_grid,Training_set_L=Training_set_L_h,
                                                                              Training_set=Training_set_h,workgrid_p,Calibration_L=Calibration_set_L_h,
                                                                              Calibration_set=Calibration_set_h)
      C_set_proposed<-Conf_Set[[as.character(setting)]][[as.character(n)]][[i]]
      Data_C_x<-c();Data_C_y<-list()
      for (j in 1:length(C_set_proposed)) {
        Data_C_x[j]<-C_set_proposed[[j]]$x;
        Data_C_y[[j]]<-C_set_proposed[[j]]$y
      }
      Cov_temp<-c();Len_tmp<-c()
      gap<-met_sphere(Search_grid[[1]]$y,Search_grid[[2]]$y)
      for (j in 1:length(work_grid_x)) {
        Ind_test<-which(abs(testing_x-work_grid_x[j])<=(work_grid_x[2]-work_grid_x[1] )/2 )
        C_set<-Data_C_y[which(Data_C_x==work_grid_x[j])]
        if(length(C_set)==0){
          Cov_temp[j]<-0
          Len_tmp[j]<-0
        }else{
          Res<-unlist(lapply(Testing_L$y[Ind_test], function(x) belongs_to_shpere(x,C_set,gap*sqrt(3)/(2*sqrt(2) ))))
          Cov_temp[j]<-sum(Res)/length(Res)
          Len_tmp[j]<-(length(C_set)/length(search_grid))*4*pi
        }
      }
      Coverage[[as.character(setting)]][[as.character(n)]][[i]]<-Cov_temp
      Length[[as.character(setting)]][[as.character(n)]][[i]]<-Len_tmp
    }
  }
  return(list(Coverage=Coverage,Length=Length,MSE=MSE_beta,Beta=Beta))
}
system.time({
  cl <- makeCluster(getOption("cl.cores", 50));
  res<-parLapply(cl, 1:200,Get_R)
  stopCluster(cl)
  rm(cl)
})

save(res, file = paste0("Muti_S2_a_siga05.RData")) 
