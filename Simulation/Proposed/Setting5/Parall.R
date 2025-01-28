library(parallel)
Get_R<-function(r){
  set.seed(r)
  library(gridExtra)
  library(truncnorm)
  #####
  #source('~/Documents/Hans/object conformal/simulation_XY_R/utils.R')
  source('utils.R')
  Coverage<-list();Length<-list();Conf_Set<-list();bd_L<-list()
  for (n in c(400,800,1600)) {
    work_grid_y<-seq(0,1,0.01)
    ep<-matrix(0,nrow = n,ncol = length(work_grid_y))#error
    for (i in 1:n) {
      ep[i,]<-Otime(runif(1,-0.47,0.47),qbeta(work_grid_y,2,2),work_grid_y)
    }
    work_grid_x<-seq(-0.9,0.9,0.02)
    workgrid_p<-seq(0,1,length=101)
    Training_x<-list();Training_y<-list()
    for (i in 1:n) {
      Training_x[[i]]<-runif(1,-1,1)
      Training_y[[i]]<-approx(work_grid_y,ep[i,],qtruncnorm(work_grid_y,a=0,b=1,mean = 0.8*f(Training_x[[i]]),sd=0.5) ,rule=2)$y
    }
    Training_set<-list(x=Training_x,y=Training_y)
    Training_set_L<-list()
    for (i in 1:n) {
      Training_set_L[[i]]<-list(x=Training_x[[i]],y=Training_y[[i]])
    }
    Calibration_x<-list();Calibration_y<-list()
    ep<-matrix(0,nrow = n,ncol = length(work_grid_y))#error
    for (i in 1:n) {
      ep[i,]<-Otime(runif(1,-0.47,0.47),qbeta(work_grid_y,2,2),work_grid_y)
    }
    for (i in 1:n) {
      Calibration_x[[i]]<-runif(1,-1,1)
      Calibration_y[[i]]<-approx(work_grid_y,ep[i,],qtruncnorm(work_grid_y,a=0,b=1,
                                                               mean = 0.8*f(Calibration_x[[i]]),sd=0.5) ,rule=2)$y
    }
    Calibration_set<-list(x=Calibration_x,y=Calibration_y)
    Calibration_L<-list()
    for (i in 1:n) {
      Calibration_L[[i]]<-list(x=Calibration_x[[i]],y=Calibration_y[[i]])
    }
    #testing_x<-runif(10000,-1,1);Testing_y<-list()
    testing_x<-rep(work_grid_x,each=100);Testing_y<-list()
    ep<-matrix(0,nrow = length(testing_x),ncol = length(work_grid_y))#error
    for (i in 1:length(testing_x)) {
      ep[i,]<-Otime(runif(1,-0.47,0.47),qbeta(work_grid_y,2,2),work_grid_y)
    }
    for (i in 1:length(testing_x)) {
      Testing_y[[i]]<-approx(work_grid_y,ep[i,],qtruncnorm(work_grid_y,a=0,b=1,
                                                           mean = 0.8*f(testing_x[i]),sd=0.5) ,rule=2)$y
    }
    Testing_L<-list(x=testing_x,y=Testing_y)
    Testing_set<-list()
    for (i in 1:length(Testing_L$x)) {
      Testing_set[[i]]<-list(x=Testing_L$x[i],y=Testing_L$y[i])
    }
    #### search grid ####
    ## set search grid
    mu_s<-seq(0,1,length=101)
    sigma_s<-seq(0.3,0.7,length=21)
    Search_grid<-list()
    {id<-1
      for (i in 1:length(work_grid_x)) {
        for (k in 1:length(mu_s)) {
          for (j in 1:length(sigma_s)) {
            Search_grid[[id]]<-list(x=work_grid_x[i],y=qtruncnorm(work_grid_y,a=0,b=1,
                                                                  mean = mu_s[k],sd=sigma_s[j] ))
            id<-id+1
          }
        }
      }}
    Coverage[['1']][[as.character(n)]]<-list()
    Length[['1']][[as.character(n)]]<-list()
    Conf_Set[['1']][[as.character(n)]]<-list()
    ##### Proposed scores #####
    #####use cross validation to choose tuning parameter
    if(n==200){bd_cand<-c(seq(0.06,0.1,length=5))}
    if(n==800){bd_cand<-c(seq(0.04,0.08,length=5))}
    if(n>800){bd_cand<-c(seq(0.02,0.4,length=5))}
    len_of_C<-c()
    for (i in 1:length(bd_cand)) {
      Conf_Set[['1']][[as.character(n)]][[i]]<-Get_conf_set(h=bd_cand[i],metric=met_W2,R_depth=R_n,
                                                            Search_grid=Search_grid,Training_set_L,
                                                            Training_set,workgrid_p,Calibration_L,Calibration_set)
      C_set_proposed<-Conf_Set[['1']][[as.character(n)]][[i]]
      Data_C_x<-c();Data_C_y<-list()
      for (j in 1:length(C_set_proposed)) {
        Data_C_x[j]<-C_set_proposed[[j]]$x;
        Data_C_y[[j]]<-C_set_proposed[[j]]$y
      }
      Cov_temp<-c();Len_tmp<-c()
      gap<-met_W2(Search_grid[[1]]$y,Search_grid[[2]]$y)
      for (j in 1:length(work_grid_x)) {
        Ind_test<-which(abs(testing_x-work_grid_x[j])<=(work_grid_x[2]-work_grid_x[1] )/2 )
        C_set<-Data_C_y[which(Data_C_x==work_grid_x[j])]
        if(length(C_set)==0){
          Cov_temp[j]<-0
          Len_tmp[j]<-0
        }else{
          Res<-unlist(lapply(Testing_L$y[Ind_test], function(x) belongs_to_W2(x,C_set,gap)))
          Cov_temp[j]<-sum(Res)/length(Res)
          Len_tmp[j]<-(length(C_set)/(length(Search_grid)/length(work_grid_x)))
        }
      }
      Coverage[['1']][[as.character(n)]][[i]]<-Cov_temp
      Length[['1']][[as.character(n)]][[i]]<-Len_tmp
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

save( res, file = paste0("W2_proposed_May28_2.RData") ) 






#### plot the conformal set ####


