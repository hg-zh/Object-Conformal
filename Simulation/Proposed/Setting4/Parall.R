library(parallel)
Get_R<-function(r){
  set.seed(r)
  library(gridExtra)
  #####
  #source('~/Documents/Hans/object conformal/simulation_XY_R/utils.R')
  source('utils.R')
  Coverage<-list();Length<-list();Conf_Set<-list();bd_L<-list()
  for (n in c(200,400,800,1600)) {
    #n=500
    #### setting 1: unimodal ####
    work_grid_x<-seq(-0.9,0.9,0.05)
    workgrid_p<-seq(0,3.5,length=101)
    Training_x<-list();Training_y<-list()
    for (i in 1:n) {
      Training_x[[i]]<-runif(1,-1,1)
      # mu<-c(1,0,0)
      # V<-c(0,rnorm(2,0,0.2))
      mu<-as.vector(f_sphere(Training_x[[i]])) 
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
      Calibration_x[[i]]<-runif(1,-1,1)
      # mu<-c(1,0,0)
      # V<-c(0,rnorm(2,0,0.2))
      mu<-as.vector(f_sphere(Calibration_x[[i]])) 
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
    #### search grid ####
    ## set search grid
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
    if(n==200){bd_cand<-c(seq(0.1,0.25,length=20))}
    if(n==400){bd_cand<-c(seq(0.06,0.21,length=20))}
    if(n==800){bd_cand<-c(seq(0.04,0.19,length=20))}
    if(n>800){bd_cand<-c(seq(0.02,0.17,length=20))}
    len_of_C<-c()
    for (i in 1:length(bd_cand)) {
      Conf_Set[['1']][[as.character(n)]][[i]]<-Get_conf_set(h=bd_cand[i],metric=met_sphere,R_depth=R_n,
                                    Search_grid=Search_grid,Training_set_L=Training_set_L,
                                     Training_set=Training_set,workgrid_p,Calibration_L=Calibration_set_L,
                                         Calibration_set=Calibration_set)
      C_set_proposed<-Conf_Set[['1']][[as.character(n)]][[i]]
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

save( res, file = paste0("Sphere_proposed_May28.RData") ) 






#### plot the conformal set ####
library(plotly)
C_set_proposed<-Conf_Set[['1']][[as.character(500)]][[1]]
C_set_list_x<-list();C_set_list_y<-list()
for(i in 1:length(C_set_proposed)){
  C_set_list_x[[i]]<-C_set_proposed[[i]]$x
  C_set_list_y[[i]]<-C_set_proposed[[i]]$y
}
SG_plot<-data.frame()

for (j in 1:length(work_grid_x)) {
  DF_tmp<-search_grid
  DF_tmp$time<-rep(work_grid_x[j],length(search_grid$x))
  SG_plot<-rbind(SG_plot,DF_tmp)
}

df_points<-data.frame(x=matrix(unlist(Training_set$y),3,500)[1,],y=matrix(unlist(Training_set$y),3,500)[2,],
                      z=matrix(unlist(Training_set$y),3,500)[3,],time=ceiling((unlist(Training_set$x)-0.1)*5)/5 )


df_C<-data.frame(x=matrix(unlist(C_set_list_y),3,length(C_set_list_y))[1,],
                 y=matrix(unlist(C_set_list_y),3,length(C_set_list_y))[2,],
                 z=matrix(unlist(C_set_list_y),3,length(C_set_list_y))[3,],time=unlist(C_set_list_x)) 


# Define the range of values for the slider
slider_values <- unique(SG_plot$time)

# Create the interactive 3D plot
p <- plot_ly() %>%
  add_trace(
    data = SG_plot,
    x = ~x,
    y = ~y,
    z = ~z,
    frame = ~time, # Define the variable that changes with the slider
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 5, color = "gray", opacity = 0.1)
  ) %>%
  add_trace(
    data = df_points,
    x = ~x,
    y = ~y,
    z = ~z,
    frame = ~time,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 5, color = "blue", opacity = 1)
  ) %>%
  add_trace(
    data = df_C,
    x = ~x,
    y = ~y,
    z = ~z,
    frame = ~time,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 5, color = "yellow", opacity = 0.5)
  ) %>%
  layout(showlegend = FALSE) %>%
  animation_slider( # Add a slider
    currentvalue = list(prefix = "Time: "),
    tickprefix = "Time: ",
    ticks = "outside",
    tickvals = slider_values,
    ticktext = slider_values,
    font = list(size = 14)
  )
p



for (i in 1:1100) {
  if(abs(Testing_L$y[[i]][1]-0.4584)<=0.01  ){
    print(i)
  }
}
Testing_L$y[[135]]


