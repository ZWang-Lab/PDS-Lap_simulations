library(MASS)
library(gtools)
library(glmnet)
library(nnet)
library(dplyr)
library(gridExtra)
library(pROC)

est_binary_laplacian_newton=function(X,Y,Lap,lambdal,
                                     beta0=rep(0,dim(X)[2]),
                                     maxiter=10000,tol=10^(-6),rate=0.1){
  fit_binary_laplacian=function(beta0,X,Y,Lap,lambdal){
    eta=X%*%beta0
    mu0 <- inv.logit(eta)
    dfdbeta=-t((Y-mu0))%*%X/(dim(X)[1])+c(2*lambdal*Lap%*%beta0)
    halfX=sqrt(c((1-mu0)*mu0))*X
    d2fdbeta2=t(halfX)%*%halfX/(dim(X)[1])
    d2fdbeta2=d2fdbeta2+2*lambdal*Lap
    # calculate inverse
    temp_inv = try(solve(d2fdbeta2))
    if(inherits(temp_inv, "try-error")){
      return("Fail to inverse")
    }
    beta=beta0-rate*c(temp_inv%*%t(dfdbeta))
    return(beta)
  }
  
  n=0
  while(n<maxiter){
    beta=fit_binary_laplacian(beta0,X,Y,Lap,lambdal)
    if(length(beta)==1) {
      print(beta)
      break
    }
    if(max(abs(beta0-beta))<=tol*max(abs(beta0))) break
    n=n+1;beta0=beta
  }
  return(beta)
}

est_binary_laplacian_cv=function(est_func,X,Y,Lap,lambda_list,k_folds=5,this.seed=123){
  # Create cross-validation folds
  set.seed(this.seed)  # For reproducibility
  folds <- sample(rep(1:k_folds, length.out = nrow(X)))
  Lap_full = matrix(0,ncol = ncol(X),nrow = ncol(X))
  Lap_full[(ncol(X)-nrow(Lap)+1):ncol(X),
           (ncol(X)-nrow(Lap)+1):ncol(X)]=Lap
  
  # Initialize a vector to store performance metrics for each lambda
  performance <- numeric(length(lambda_list))
  
  # Loop over each lambda value
  for(i in c(1:length(lambda_list))){
    lambdal <- lambda_list[i]
    
    # Initialize a vector to store the performance metric for each fold
    cv_performance <- numeric(k_folds)
    
    # Perform k-fold cross-validation
    for(fold in 1:k_folds){
      cat("Lambda ",lambdal,"Fold ",fold,".\n")
      # Split the data into training and validation sets
      test_indices <- which(folds == fold)
      train_indices <- setdiff(seq_len(nrow(X)), test_indices)
      
      X_train <- X[train_indices, ]
      Y_train <- Y[train_indices]
      
      X_test <- X[test_indices, ]
      Y_test <- Y[test_indices]
      
      remove_col = which(colSums(X_train!=0)==0)
      if(length(remove_col)>0){
        X_train = X_train[,-remove_col]
        X_test = X_test[,-remove_col]
        this.Lap = Lap_full[-remove_col,-remove_col]
      }else{
        this.Lap = Lap_full
      }
      
      # Fit the model on the training data with the current lambda
      this.fit <- est_func(X=X_train, Y=Y_train,Lap=this.Lap, 
                           lambdal=lambdal)
      
      if(length(this.fit)==1){
        return(this.fit)
      }else{
        model = this.fit
      }
      # Make predictions on the validation data
      predictions <- c(inv.logit(X_test%*%model))
      # criteria with AUC
      roc_obj <- roc(Y_test,predictions,levels = c(0, 1), direction = "<")
      auc_value_pROC <-  auc(roc_obj)
      cv_performance[fold] <- auc_value_pROC  # Calculate the misclassification error
    }
    # Store the average performance across folds for the current lambda
    performance[i] <- mean(cv_performance)
  }
  return(list(performance=performance,
              lambda_grid=lambda_list,
              penalty_matrix=Lap_full))
}