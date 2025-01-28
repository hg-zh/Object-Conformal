library(ggplot2)
library(reshape2)
library(cowplot)
library(grid)
load("~/Documents/Hans/object conformal/taxi/Weekday_2023.RData")
source("utlis.R")
set.seed(2024)
#### time rage 4 - 20 
x_min=4;x_max<-20
Data_raw<-Data_wd
Data_wd<-list()
for (i in 1:length(Data_raw)) {
  if((!any(is.na(Data_raw[[i]]$y)))&(Data_raw[[i]]$x<=x_max)&(Data_raw[[i]]$x>=x_min)  ){
    Data_wd[[length(Data_wd)+1]]<-Data_raw[[i]]
  }
}
Wd_ind<-c()
for (i in 1:length(Data_wd)) {
  Wd_ind[length(Wd_ind)+1]<-Data_wd[[i]]$date
}
length(unique(Wd_ind))


n<-length(Data_wd)
n_train<-floor(4*n/10);n_cal<-floor(4*n/10);n_test<-n-n_train-n_cal
n_perm<-sample.int(n,n)
id_train<-sort(n_perm[1:n_train]);id_cal<-sort(n_perm[(n_train+1):(n_train+n_cal)])
id_test<-sort(n_perm[(n_train+n_cal+1):n])
Training_L<-Data_wd[id_train];Calibration_L<-Data_wd[id_cal];Testing_L<-Data_wd[id_test]
Training_x<-list();Training_y<-list()
for (i in 1:n_train) {
  Training_x[[i]]<-Training_L[[i]]$x
  Training_y[[i]]<-Training_L[[i]]$y
}
Training_set<-list(x=Training_x,y=Training_y)
Calibration_x<-list();Calibration_y<-list()
for (i in 1:n_cal) {
  Calibration_x[[i]]<-Calibration_L[[i]]$x
  Calibration_y[[i]]<-Calibration_L[[i]]$y
}
Calibration_set<-list(x=Calibration_x,y=Calibration_y)
Testing_x<-list();Testing_y<-list()
for (i in 1:n_test) {
  Testing_x[[i]]<-Testing_L[[i]]$x
  Testing_y[[i]]<-Testing_L[[i]]$y
}
Testing_set<-list(x=Testing_x,y=Testing_y)
#### load holiday data ####
load("~/Documents/Hans/object conformal/taxi/Holiday_2023.RData")
Data_raw_ho<-Data_ho
Data_ho<-list()
for (i in 1:length(Data_raw_ho)) {
  if((!any(is.na(Data_raw_ho[[i]]$y)))&(Data_raw_ho[[i]]$x>=x_min)&(Data_raw_ho[[i]]$x<=x_max)  ) {
    Data_ho[[length(Data_ho)+1]]<-Data_raw_ho[[i]]
  }
}
Testing_L_ho<-Data_ho

Testingho_x<-list();Testingho_y<-list()
for (i in 1:n_test) {
  Testingho_x[[i]]<-Testing_L_ho[[i]]$x
  Testingho_y[[i]]<-Testing_L_ho[[i]]$y
}
Testing_set_ho<-list(x=Testingho_x,y=Testingho_y)


#### conformal ####
h<-2
metric<-met_FB
workgrid_p<-seq(0,5,length=101)
dist_profile <- t(matrix(unlist(lapply(Training_L, function(x) Dist_profile(x$y,x$x, workgrid_p, Training_set,metric,h))),
                         length(workgrid_p) ,length(Training_L) ) )
inv_dist_profile <- apply(dist_profile, 1, function(x) Get_Inverse(workgrid_p, x)$y)
Q_train <- t(inv_dist_profile)
Ranks_on_train<-unlist(lapply(Training_L,function(x) R_n(x$y,x$x,Training_set,metric,F_mean=NULL,Q_train,h,workgrid_p)))
Ranks_on_cal<-unlist(lapply(Calibration_L,function(x) R_n(x$y,x$x,Training_set,metric,F_mean=NULL,Q_train,h,workgrid_p)))
Ranks_on_test<-unlist(lapply(Testing_L,function(x) R_n(x$y,x$x,Training_set,metric,F_mean=NULL,Q_train,h,workgrid_p)))
Ranks_on_test_ho<-unlist(lapply(Testing_L_ho,function(x) R_n(x$y,x$x,Training_set,metric,F_mean=NULL,Q_train,h,workgrid_p)))
S_cal_list<-list()
for (i in 1:length(Calibration_L)) {
  S_cal_list[[i]]<-list(z=Ranks_on_cal[i],x=Calibration_set$x[[i]])
}
S_test_list<-list()
for (i in 1:n_test) {
  S_test_list[[i]]<-list(z=Ranks_on_test[i],x=Testing_set$x[[i]])
}
S_test_list_ho<-list()
for (i in 1:n_test) {
  S_test_list_ho[[i]]<-list(z=Ranks_on_test_ho[i],x=Testing_set_ho$x[[i]])
}
Scores_on_cal<-unlist(lapply(S_cal_list, function(x) H(x$z,x$x,Training_set,Ranks_on_train,h) ))
Scores_on_test<-unlist(lapply(S_test_list, function(x) H(x$z,x$x,Training_set,Ranks_on_train,h) ))
Scores_on_test_ho<-unlist(lapply(S_test_list_ho, function(x) H(x$z,x$x,Training_set,Ranks_on_train,h) ))
q_H<-quantile(Scores_on_cal,0.95*(1+1/length(Calibration_L)),names = FALSE)
sum(Scores_on_test<=q_H)/length(Scores_on_test)
sum(Scores_on_test_ho<=q_H)/length(Scores_on_test_ho)
#### plot the conditional coverage
l=101
df_cd<-data.frame(x=seq(x_min,x_max,length=l),wd=rep(0,l),ho=rep(0,l))
for (j in 1:l) {
  x<-df_cd$x[j]
  Index_wd<-which((Testing_set$x>=x-h)&(Testing_set$x<=x+h) )
  df_cd$wd[j]<-sum(Scores_on_test[Index_wd]<=q_H)/length(Scores_on_test[Index_wd])
  Index_ho<-which((Testing_set_ho$x>=x-h)&(Testing_set_ho$x<=x+h) )
  df_cd$ho[j]<-sum(Scores_on_test_ho[Index_ho]<=q_H)/length(Scores_on_test_ho[Index_ho])
}
plot <- ggplot(data = df_cd) + 
  geom_line(aes(x = x, y = ho), color = "black", size = 1, linetype = "dashed") +
  geom_hline(yintercept = 0.95, color = "red", size = 0.6, linetype = "solid") + 
  geom_line(aes(x = x, y = wd), color = "black", size = 1, linetype = "solid") + 
  theme_minimal() + ylim(0.4,1)+
  labs(x = "hours in a day",
       y = "conditional coverage") +
  scale_color_manual(values = c("wd" = "blue", "ho" = "red")) +
  theme(
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)  # Adds a black border
  )
#### heatmap of frechet mean ####
time_of_Int<-15;h<-1
Index_time<-which((Testing_set$x>=time_of_Int-h)&(Testing_set$x<=time_of_Int+h))
F_mean<-matrix(0,13,13)
for(i in 1:length(Index_time)){
  F_mean<-F_mean+Training_L[[Index_time[i]]]$y/length(Index_time)
}
F_mean_long <- melt(F_mean)

# Naming the columns for clarity (optional)
colnames(F_mean_long) <- c("X", "Y", "Value")
F_mean_long$group<-rep("F-mean",length(F_mean_long$X))

# Create the heatmap
p_mean<-ggplot(F_mean_long, aes(X, Y, fill = Value)) + 
  geom_tile() + # Use geom_tile() for heatmap
  scale_fill_gradient(low = "blue", high = "red") + # Customize color gradient
  theme_minimal() + # Optional: use a minimal theme
  labs(x = "regions", y = "regions", fill = NULL)+
  theme(legend.position="none")

id_low<-which.min(Ranks_on_train[Index_time])

F_low<-melt(Training_L[[Index_time[id_low] ]]$y)
colnames(F_low) <- c("X", "Y", "Value")
F_low$group<-rep("Lowest",length(F_low$X))

p_low<-ggplot(F_low, aes(X, Y, fill = Value)) + 
  geom_tile() + # Use geom_tile() for heatmap
  scale_fill_gradient(low = "blue", high = "red") + # Customize color gradient
  theme_minimal() + # Optional: use a minimal theme
  labs(x = "regions", y = "regions", fill = NULL)+
  theme(legend.position="none")

id_high<-which.min(abs(Ranks_on_train[Index_time]-q_H))
F_high<-melt(Training_L[[Index_time[id_high] ]]$y)
colnames(F_high) <- c("X", "Y", "Value")
F_high$group<-rep("Highest",length(F_high$X))

p_high<-ggplot(F_high, aes(X, Y, fill = Value)) + 
  geom_tile() + # Use geom_tile() for heatmap
  scale_fill_gradient(low = "blue", high = "red") + # Customize color gradient
  theme_minimal() + # Optional: use a minimal theme
  labs(x = "regions", y = "regions", fill = NULL)+
  theme(legend.position="none")

plot_grid(p_mean,p_low,p_high)



df_plot <- rbind(F_mean_long, F_low,F_high)
df_plot$group <- factor(df_plot$group, levels = c("F-mean", "Lowest", "Highest"))

# Now, plot using ggplot2 and ensure facets are ordered as per factor levels
heatmap_plot <- ggplot(df_plot, aes(x = X, y = Y, fill = Value)) +
  geom_tile(width = 1.01, height = 1.01) +  # Creates the heatmap
  facet_wrap(~group, scales = "free") +  # Facets by group, ensuring order
  scale_fill_gradient(low = "blue", high = "red") +  # Sets the color gradient
  theme_minimal() +  # Applies a minimal theme
  labs(fill = NULL,x='regions',y='regions') +  # Labels the color scale
  theme(
    panel.background = element_blank(),  # Removes background
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank(),  # Removes minor grid lines
    axis.text.x = element_blank(),  # Removes X axis text
    axis.text.y = element_blank(),  # Removes Y axis text
    axis.ticks = element_blank(),  # Removes axis ticks
    panel.border = element_rect(colour = "black", fill=NA),  # Adds black border
    plot.margin = unit(c(2,2,2,2), "mm"),  # Removes margins around the plot
    panel.spacing = unit(0, "mm"),  # Removes spacing between panels
    strip.background = element_rect(colour = "black", fill=NA)  # Optional: Adds border around facet labels
  )



# Display the plot
print(heatmap_plot)




