Get_Inverse <- function(u, y) {
  min_y <- range(y)[1]
  max_y <- range(y)[2]
  min_u <- range(u)[1]
  max_u <- range(u)[2]
  temp_x <- seq(min_y, max_y, length = length(y))
  temp_y <- numeric(length = length(u))
  temp_y <- approx(y, u, temp_x, rule = 2,ties=min)$y
  return(list(x = temp_x, y = temp_y))
}
met_R1<-function(a,b){return(abs(a-b))}
met_sphere<-function(a,b){
  return(acos(min(t(a)%*%b,1) ) )
}
Dist_profile<-function(w,x,t,Training_set,metric,h){
  Index<-which((Training_set$x>=x-h)&(Training_set$x<=x+h) )
  if(length(Index)==0){
    Index<-which.min(abs(unlist(Training_set$x)-x))
    S <- unlist(lapply(Training_set$y[Index],  function(x) metric(x,w)  )) 
    W<-1
    Out <- sapply(t, function(x) sum((S <= x)*W)/sum(W))
  }else{
    S <- unlist(lapply(Training_set$y[Index],  function(x) metric(x,w)  )) 
    W<-epanechnikov_kernel((unlist(Training_set$x[Index]) -x)/h )
  }
  Out <- sapply(t, function(x) sum((S <= x)*W)/sum(W))
  return(Out)
}
R_n<-function(y,x,Training_set,metric,F_mean,Q_train,h,workgrid_p){
  Index<-which((Training_set$x>=x-h)&(Training_set$x<=x+h) )
  if(length(Index)==0){
    Index<-which.min(abs(unlist(Training_set$x)-x))
    W<-1
    Q_train_in<-Q_train[Index,]
    Q_i<-Get_Inverse(workgrid_p,Dist_profile(y,x,workgrid_p,Training_set,metric,h))$y
    R_i<-sum( abs(Q_train_in-rep(1,length(Index))%*%t(Q_i))*W )/(length(Q_i))
  }else{
    W<-epanechnikov_kernel((unlist(Training_set$x[Index]) -x)/h )
    Q_train_in<-Q_train[Index,]
    Q_i<-Get_Inverse(workgrid_p,Dist_profile(y,x,workgrid_p,Training_set,metric,h))$y
    R_i<-sum( abs(Q_train_in-rep(1,length(Index))%*%t(Q_i))*W )/(length(Q_i))
  }
  return(expit(-R_i/sum(W)) )
}
H<-function(z,x,Training_set,Ranks_on_train,h){
  Index<-which((Training_set$x>=x-h)&(Training_set$x<=x+h) )
  W<-epanechnikov_kernel((unlist(Training_set$x[Index]) -x)/h )
  return(sum((Ranks_on_train[Index]<=z)*W)/sum(W) )
}
Get_conf_set<-function(h=0.1,metric=met_sphere,R_depth,Search_grid,Training_set_L,
                       Training_set,workgrid_p,Calibration_L,Calibration_set){
  dist_profile <- t(matrix(unlist(lapply(Training_set_L, function(x) Dist_profile(x$y,x$x, workgrid_p, Training_set,metric,h))),
                           length(workgrid_p) ,length(Training_set_L) ) )
  inv_dist_profile <- apply(dist_profile, 1, function(x) Get_Inverse(workgrid_p, x)$y)
  Q_train <- t(inv_dist_profile)
  Ranks_on_train<-unlist(lapply(Training_set_L,function(x) R_depth(x$y,x$x,Training_set,metric,F_mean=NULL,Q_train,h,workgrid_p)))
  Ranks_on_cal<-unlist(lapply(Calibration_L,function(x) R_depth(x$y,x$x,Training_set,metric,F_mean=NULL,Q_train,h,workgrid_p)))
  S_cal_list<-list()
  for (i in 1:length(Calibration_L)) {
    S_cal_list[[i]]<-list(z=Ranks_on_cal[i],x=Calibration_set$x[[i]])
  }
  Scores_on_cal<-unlist(lapply(S_cal_list, function(x) H(x$z,x$x,Training_set,Ranks_on_train,h) ))
  q_H<-quantile(Scores_on_cal,0.1*(1+1/length(Calibration_L)),names = FALSE)
  ##serach grid
  Ranks_on_SG<-unlist(lapply(Search_grid,function(x) R_depth(x$y,x$x,Training_set,metric,F_mean=NULL,Q_train,h,workgrid_p)))
  S_SG_list<-list()
  for (i in 1:length(Search_grid)) {
    S_SG_list[[i]]<-list(z=Ranks_on_SG[i],x=Search_grid[[i]]$x)
  }
  Scores_on_SG<-unlist(lapply(S_SG_list, function(x) H(x$z,x$x,Training_set,Ranks_on_train,h) ))
  C_set<-Search_grid[which(Scores_on_SG>=q_H)] 
  return(C_set)
}
expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}
f<-function(x){return((x-1)^2*(x+1) )}
g<-function(x){
  out<-rep(0,length(x))
  id<-which(x>=0)
  out[id]<-sqrt(x[id])
  return(2*out)
}
g_sphere<-function(x){
  out<-matrix(0,length(x),3)
  for (i in which(x>0)) {
    if(rbinom(1,1,0.5)){
      out[i,]<-c(0,cos(pi*x[i]/2),sin(pi*x[i]/2))
    }else{
      out[i,]<-c(0,cos(pi*x[i]/2),-sin(pi*x[i]/2))
    }
  }
  return(out)
}




belongs_to_set <- function(datapoints, input_set, gap = 1) {
  # Calculate the differences between consecutive elements in the input_set
  diffs <- diff(input_set)
  
  # Find the indices where the difference is greater than the gap
  interval_indices <- c(0, which(diffs > gap+sqrt(.Machine$double.eps)), length(input_set))
  
  # Define a function to check if a single datapoint belongs to any interval
  check_datapoint <- function(datapoint, input_set, interval_indices) {
    start <- input_set[interval_indices[-length(interval_indices)] + 1]
    end <- input_set[interval_indices[-1]]
    any(datapoint >= start-sqrt(.Machine$double.eps) & datapoint <= end+sqrt(.Machine$double.eps))
  }
  
  # Apply the check_datapoint function to each element of datapoints
  belongs <- sapply(datapoints, check_datapoint, input_set, interval_indices)
  
  return(belongs)
}
met_SPD <- function(P,Q) {
  return(norm(logm(solve(sqrtm(P))%*%Q%*%solve(sqrtm(P))) ,"F") )
}
mu_SPD<-function(x){
  if(x<=0){
    if(rbinom(1,1,0.5)){
      return(matrix(c(1,(1+x)/2,(1+x)/2,1),2,2))
    }else{
      return(matrix(c(1,(-x-1)/2,(-x-1)/2,1),2,2))
    }
  }else{
    if(rbinom(1,1,0.5)){
      return(matrix(c(1+x,(1+x)/2,(1+x)/2,1),2,2))
    }else{
      return(matrix(c(1,(-x-1)/2,(-x-1)/2,1+x),2,2))
    }
  }
}


belongs_to_shpere<-function(datapoint,input_set,tolerance){
  for (i in 1:length(input_set)) {
    if(met_sphere(datapoint,input_set[[i]] )<=tolerance){
      return(TRUE)
    }
  }
  return(FALSE)
}
epanechnikov_kernel <- function(u) {
  ifelse(abs(u) <= 1, 3/4 * (1 - u^2), 0)
}

local_linear_smoother_CV <- function(X, Y, h = NULL, x_vec, lambda = 1e-6, cv_folds = 5, h_candidates = seq(0.1, 2, by=0.1)) {
  
  # Epanechnikov kernel
  epanechnikov_kernel <- function(u) {
    ifelse(abs(u) <= 1, (3/4)*(1 - u^2), 0)
  }
  
  # Internal function for performing cross-validation
  choose_bandwidth <- function(X, Y, h_candidates) {
    n <- length(X)
    fold_indices <- sample(rep(1:cv_folds, length.out=n))
    best_mse <- Inf
    best_h <- NULL
    
    # Loop over candidate bandwidths
    for (h in h_candidates) {
      total_mse <- 0
      for (i in 1:cv_folds) {
        # Split data into training and validation sets
        train_idx <- which(fold_indices != i)
        valid_idx <- which(fold_indices == i)
        X_train <- X[train_idx]
        Y_train <- Y[train_idx]
        X_valid <- X[valid_idx]
        Y_valid <- Y[valid_idx]
        
        # Predict using local_linear_smoother without CV
        Y_pred <- local_linear_smoother(X_train, Y_train, h, X_valid) 
        fold_mse <- mean((Y_pred - Y_valid)^2, na.rm=TRUE)
        total_mse <- total_mse + fold_mse
      }
      
      # Calculate average MSE
      avg_mse <- total_mse / cv_folds
      
      # Update best bandwidth if necessary
      if (avg_mse < best_mse) {
        best_mse <- avg_mse
        best_h <- h
      }
    }
    return(best_h)
  }
  
  # If h is NULL, perform cross-validation to choose bandwidth
  if (is.null(h) && cv_folds > 0) {
    h <- choose_bandwidth(X, Y, h_candidates)
  }
  
  return(local_linear_smoother(X, Y, h, x_vec, lambda = 1e-6))
}
local_linear_smoother <- function(X, Y, h, x_vec, lambda = 1e-6) {
  # Epanechnikov kernel
  epanechnikov_kernel <- function(u) {
    ifelse(abs(u) <= 1, (3/4)*(1 - u^2), 0)
  }
  
  smoothed_values <- sapply(x_vec, function(x) {
    # Calculate weights using kernel function
    weights <- sapply((x - X) / h, epanechnikov_kernel)
    
    # Check if weights are almost all zeros
    if (max(weights) < 1e-10) {
      return(Y[which.min(abs(x-X))])
    }
    
    # Build design matrix for local linear regression
    Z <- cbind(1, X - x)
    
    # Regularization matrix
    reg <- lambda * diag(2)
    
    # Solve for local coefficients with regularization
    b <- solve(t(Z) %*% (Z * weights) + reg, t(Z) %*% (Y * weights))
    
    # a(x) is the first coefficient
    estimate <- b[1]
    
    return(estimate)
  })
  
  return(smoothed_values)
}

ll_frechet <- function(X, Y ,h, x_vec) {
  # Epanechnikov kernel
  Ln<-function(params,Y,weights,mu,Z) {
    d_Y<- unlist(lapply(Y,function(x) met_sphere(x,c(sin(params[1])*cos(params[2]),
                                                     sin(params[1])*sin(params[2]),cos(params[1]) ))))
    return( as.numeric(t((weights)*(mu[3]-mu[2]*Z[,2]))%*%d_Y    ) )
  }
  smoothed_values <- lapply(x_vec, function(x) {
    # Calculate weights using kernel function
    weights <- sapply(( X-x) / h, epanechnikov_kernel)
    Z <- cbind(1, X-x,(X-x)^2 )
    mu<-t(Z)%*%weights
    
    # Check if weights are almost all zeros
    if (max(weights) < 1e-10) {
      return(unlist(Y[which.min(abs(x-X))]))
    }
    
    hat_phi <- optim(c(0, 0), function(x) Ln(x,Y,weights,mu,Z))
    phi<-hat_phi$par[1];psi<-hat_phi$par[2]
    return(c(sin(phi)*cos(psi),sin(phi)*sin(psi),cos(phi) ))
  })
  return(smoothed_values)
}

ll_frechet_CV <- function(X, Y, h = NULL, x_vec, cv_folds = 5, h_candidates = seq(0.1, 1, by=0.1)) {
  
  # Epanechnikov kernel
  
  
  # Internal function for performing cross-validation
  choose_bandwidth <- function(X, Y, h_candidates) {
    n <- length(X)
    fold_indices <- sample(rep(1:cv_folds, length.out=n))
    best_mse <- Inf
    best_h <- NULL
    
    # Loop over candidate bandwidths
    for (h in h_candidates) {
      total_mse <- 0
      for (i in 1:cv_folds) {
        # Split data into training and validation sets
        train_idx <- which(fold_indices != i)
        valid_idx <- which(fold_indices == i)
        X_train <- X[train_idx]
        Y_train <- Y[train_idx]
        X_valid <- X[valid_idx]
        Y_valid <- Y[valid_idx]
        
        # Predict using local_linear_smoother without CV
        Y_pred <- ll_frechet(X_train, Y_train, h, X_valid) 
        fold_mse <- sum(mapply(function(x,y) met_sphere(x,y),Y_pred,Y_valid ) ) 
        total_mse <- total_mse + fold_mse
      }
      
      # Calculate average MSE
      avg_mse <- total_mse / cv_folds
      
      # Update best bandwidth if necessary
      if (avg_mse < best_mse) {
        best_mse <- avg_mse
        best_h <- h
      }
    }
    return(best_h)
  }
  
  # If h is NULL, perform cross-validation to choose bandwidth
  if (is.null(h) && cv_folds > 0) {
    h <- choose_bandwidth(X, Y, h_candidates)
  }
  
  return(ll_frechet(X, Y, h, x_vec))
}

f_sphere<-function(x){
  out<-matrix(0,length(x),3)
  for (i in 1:length(x)) {
    out[i,]<-c(sin(abs(pi*x[i]/4)),cos(abs(pi*x[i]/4)),0)
  }
  return(out)
}
epanechnikov_kernel <- function(u) {
  ifelse(abs(u) <= 1, (3/4)*(1 - u^2), 0)
}


