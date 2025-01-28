source("utils.R")
load("~/Documents/Hans/object conformal/energy data/pre-processed.RData")
#### turn to the conformal dataset
set.seed(2023)
df_all<-data.frame()
years<-1990:2021;work_grid_x<-years
workgrid_p<-seq(0,2,length=101)
for (k in 1:length(years)) {
  Plot_df<-data.frame(x=numeric(length(Data_year[[k]])),y=numeric(length(Data_year[[k]])),
                      z=numeric(length(Data_year[[k]])),state=character(length(Data_year[[k]])))
  for (i in 1:length(Data_year[[k]])) {
    Plot_df$x[i]<-Data_year[[k]][[i]][1]
    Plot_df$y[i]<-Data_year[[k]][[i]][2]
    Plot_df$z[i]<-Data_year[[k]][[i]][3]
    Plot_df$state[i]<-names(Data_year[[k]])[i]
  }
  Plot_df$time<-rep(years[k],length(Plot_df$x))
  df_all<-rbind(df_all,Plot_df)
}
ind_train<-sort(sample.int(length(df_all$x),floor(length(df_all$x)/2)))
Data_train<-df_all[ind_train,]
Data_calibration<-df_all[-ind_train,]
Training_x<-list();Training_y<-list()
for (i in 1:length(Data_train$x)) {
  Training_x[[i]]<-Data_train$time[i]
  Training_y[[i]]<-c(Data_train$x[i],Data_train$y[i],Data_train$z[i])
}
Training_set<-list(x=Training_x,y=Training_y)
Training_set_L<-list()
for (i in 1:length(Data_train$x)) {
  Training_set_L[[i]]<-list(x=Training_x[[i]],y=Training_y[[i]])
}
Calibration_x<-list();Calibration_y<-list()
for (i in 1:length(Data_calibration$x)) {
  Calibration_x[[i]]<-Data_calibration$time[i]
  Calibration_y[[i]]<-c(Data_calibration$x[i],Data_calibration$y[i],Data_calibration$z[i])
}
Calibration_set<-list(x=Calibration_x,y=Calibration_y)
Calibration_L<-list()
for (i in 1:length(Data_calibration$x)) {
  Calibration_L[[i]]<-list(x=Calibration_x[[i]],y=Calibration_y[[i]])
}
#### search grid ####
## set search grid
u <- seq(0, pi, length.out = 100)
v <- seq(0, 2 * pi, length.out = 100)
x <- outer(cos(u), cos(v), FUN = "*")
y <- outer(cos(u), sin(v), FUN = "*")
z <- outer(sin(u), rep(1, length(v)), FUN = "*")
search_grid <- data.frame(x = c(c(x),c(x)), y = c(c(y),c(y)) , z = c(c(z),-c(z)) )
search_grid<-search_grid[(search_grid$x>=0)&(search_grid$y>=0)&(search_grid$z>=0), ]
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
##
bd<-seq(1.5,2.5,length=11)
Len<-rep(0,length(bd))
for (i in 1:length(bd)) {
  C_set_proposed<-Get_conf_set(h=bd[i],metric=met_sphere,R_depth=R_n,
                               Search_grid=Search_grid,Training_set_L,
                               Training_set,workgrid_p,Calibration_L,Calibration_set,alpha = 0.9)
  Len[i]<-length(C_set_proposed)
}
which.min(Len)
C_set_proposed<-Get_conf_set(h=bd[which.min(Len)],metric=met_sphere,R_depth=R_n,
                             Search_grid=Search_grid,Training_set_L,
                             Training_set,workgrid_p,Calibration_L,Calibration_set,alpha = 0.9)
## plot
SG_plot<-data.frame()
for (j in 1:length(years)) {
  DF_tmp<-search_grid
  DF_tmp$time<-rep(years[j],length(search_grid$x))
  SG_plot<-rbind(SG_plot,DF_tmp)
}
PT_plot<-data.frame()
for (k in 1:length(years)) {
  Plot_df<-data.frame(x=numeric(length(Data_year[[k]])),y=numeric(length(Data_year[[k]])),
                      z=numeric(length(Data_year[[k]])),state=character(length(Data_year[[k]])))
  for (i in 1:length(Data_year[[k]])) {
    Plot_df$x[i]<-Data_year[[k]][[i]][1]
    Plot_df$y[i]<-Data_year[[k]][[i]][2]
    Plot_df$z[i]<-Data_year[[k]][[i]][3]
    Plot_df$state[i]<-names(Data_year[[k]])[i]
  }
  Plot_df$time<-rep(years[k],length(Plot_df$x))
  PT_plot<-rbind(PT_plot,Plot_df)
}
C_set_list_x<-list();C_set_list_y<-list()
for(i in 1:length(C_set_proposed)){
  C_set_list_x[[i]]<-C_set_proposed[[i]]$x
  C_set_list_y[[i]]<-C_set_proposed[[i]]$y
}
df_C<-data.frame(x=matrix(unlist(C_set_list_y),3,length(C_set_list_y))[1,],
                 y=matrix(unlist(C_set_list_y),3,length(C_set_list_y))[2,],
                 z=matrix(unlist(C_set_list_y),3,length(C_set_list_y))[3,],time=unlist(C_set_list_x)) 
##
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
    data = df_C,
    x = ~x,
    y = ~y,
    z = ~z,
    frame = ~time,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 5, color = "red", opacity = 0.3)
  ) %>%
  add_trace(
    data = PT_plot,
    x = ~x,
    y = ~y,
    z = ~z,
    frame = ~time,
    type = "scatter3d",
    mode = "markers",
    # color=~state,
    marker = list(size = 5, color="blue", opacity = 1)
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
p <- layout( p, scene = list(
                    xaxis = list(title = "fossil fuels"),
                    yaxis = list(title = "natural gas"),
                    zaxis = list(title = "renewable energy")
                  ))
p






