library(glmnet)
library(gtools)
library(MASS)
library(nnet)
library(dplyr)
library(gridExtra)
library(pROC)

source("PDS-Lap.R")
source("generate_data.R")

set.seed(1)
# fixed coefficients 
# covariate effect on treatment selection
alphas0 <- matrix(round(rnorm(3*22,0,1),2),
                  nrow = 3,byrow = T)
temp = c(25,2,9,28,8,
         5,1,
         5,2,
         13,1,1)/100

# probability of drug in each drug class
gammas <- matrix(c(temp[1:5]/sum(temp[1:5]), rep(0,7),  # SSRI
                   rep(0,5), temp[6:7]/sum(temp[6:7]), rep(0,5),          #SNRI
                   rep(0,7), temp[8:9]/sum(temp[8:9]), rep(0,3),    #OA
                   rep(0,9), temp[10:12]/sum(temp[10:12])),         #TCA
                 nrow = 4,byrow = T)

# drug effect
betas=c(-2.513784,-2.513784,-2.50111,-2.528043,-2.522842,
        -1.174818,-1.144003,
        -0.5319548,-0.5141345,
        0.4884806,0.4888055,0.5057889)

trt.drug.class <- as.factor(c(rep("Class 1",5),
                              rep("Class 2",2),
                              rep("Class 3",2),
                              rep("Class 4",3)))

# covariate effects on outcome
set.seed(100)
omega = abs(round(rnorm(22,0,1.5),2))*rep(c(1,0.2),11) 

set.seed(1)
# Different scenario settings 
# Scenario No. - # outcome / # treatment / # confounder

# Scenario E - 16/2/0
alphas1 <- alphas0
omega1 <- omega
alphas1[,c(1:6,8:13,15:22)] <- 0 
omega1[c(7,10,14,21)] <- 0
alphas1 <- cbind(c(-2.45,-2.33,-1.5),alphas1)

# Scenario F - 16/2/2
alphas2 <- alphas0
omega2 <- omega
alphas2[,c(2:8,10:22)] <- 0 
alphas2 <- cbind(c(-2.55,-2.68,-1.6),alphas2)
omega2[c(7,10,14,21)] <- 0

# Scenario G - 16/6/6
alphas3 <- alphas0
omega3 <- omega
alphas3[,c(2:8,10,12:16,18,20,21)] <- 0 
omega3[c(7,10,14,21)] <- 0
alphas3 <- cbind(c(-2.66,-2.70,-1.73),alphas3)

# Scenario H - 16/16/16
alphas4 <- alphas0
omega4 <- omega
alphas4[,c(7,10,14,21)] <- 0 
omega4[c(7,10,14,21)] <- 0
alphas4 <- cbind(c(-2.95,-2.15,-1.85),alphas4)

# Scenario A - 6/2/0
alphas5 <- alphas0
omega5 <- omega
alphas5[,c(1,3:6,8:13,15,17,19,21,22)] <- 0 
omega5[c(2:8,10,12:16,18,20,21)] <- 0
alphas5 <- cbind(c(-2.55,-1.9,-1.26),alphas5)

# Scenario B - 6/2/2
alphas6 <- alphas0
omega6 <- omega
alphas6[,c(2:8,10:22)] <- 0 
omega6[c(2:8,10,12:16,18,20,21)] <- 0
alphas6 <- cbind(c(-2.55,-2.68,-1.6),alphas6)

# Scenario C - 6/6/6
alphas7 <- alphas0
omega7 <- omega
alphas7[,c(2:8,10,12:16,18,20,21)] <- 0 
omega7[c(2:8,10,12:16,18,20,21)] <- 0
alphas7 <- cbind(c(-2.66,-2.70,-1.73),alphas7)

# Scenario D - 6/16/6
alphas8 <- alphas0
omega8 <- omega
alphas8[,c(7,10,14,21)] <- 0 
omega8[c(2:8,10,12:16,18,20,21)] <- 0
alphas8 <- cbind(c(-2.95,-2.15,-1.85),alphas8)


prevalence=c(4,8,15)
scenarios=c("a","b","c","d","e","f","g","h")
alphas=list(alphas1,alphas2,alphas3,alphas4,alphas5,alphas6,alphas7,alphas8)
omegas=list(omega1,omega2,omega3,omega4,omega5,omega6,omega7,omega8)
intercept_main=c(-2.82,-2.7,-2.8,-2.8,-2.7,-2.58,-2.65,-2.65, 
                 -1.9,-1.82,-1.88,-1.88,-1.78,-1.68,-1.73,-1.73,
                 -0.96,-0.89,-0.93,-0.93,-0.84,-0.76,-0.79,-0.79)



file_path = '~/ZWang-Lab/'
sim_matrix = read.csv(paste0(file_path, 'drug_similarity.csv'))
A = as.matrix(sim_matrix )
D=diag(colSums(A))
Lap_matrix <- D - A

# Define a grid of lambda values to test
lambda_grid <- c(1e-5,5e-5,1e-4,5e-4,0.001,0.005,0.01,0.05,0.1,0.5,1,5,10) # Adjust the range and length as needed

args = commandArgs(trailingOnly=TRUE)

for(l in 1:3){
  for(j in 1:8){
    for(k in c(0,50)){
      simu_results=NULL
      
      resimu=k
      seedi <- k+l*100+j*1000 + 10000*as.numeric(args[1])+resimu*100000
      cat("args",seedi, "\n")

      set.seed(seed=seedi)
      simu_datai=simu_data(alphas[[j]],omegas[[j]],intercept_main[(l-1)*8+j],gammas,4500)
      dim(simu_datai)
      Y=simu_datai[,1]
      
      Z <- cbind(simu_datai[,3:24],
                 simu_datai[,3]*simu_datai[,c(8:24)])
      X=simu_datai[,25:36]
      C=simu_datai[,37]
      Med=apply(X,1,function(x) which(x==1))
      Z  <- as.matrix(Z)
      X <- as.matrix(X)
      fit.cov = as.matrix(cbind(Z,X))
      all_Z = as.matrix(Z)
      colnames(all_Z) = paste("X",1:(dim(all_Z)[2]),sep="")
      colnames(fit.cov) = c(colnames(all_Z),paste("d",1:12,sep=""))
      
      dim(fit.cov) # 4500x68
      
      C=as.factor(C)
      Med=as.factor(Med)
      
      save.path=paste0("~/ZWang-Lab/result_",prevalence[l],scenarios[j])
      
      # LASSO
      outcome_Lasso_cv <- cv.glmnet(x=fit.cov,
                                    y=Y,
                                    family="binomial",
                                    intercept =F,
                                    type.measure = "auc",
                                    nfolds=5,
                                    penalty.factor=c(rep(1,dim(fit.cov)[2]-12),rep(0,12)))
      
      outcome_Lasso <- glmnet(x=fit.cov,
                              y=Y,
                              family="binomial",
                              intercept =F,
                              penalty.factor=c(rep(1,dim(fit.cov)[2]-12),rep(0,12)))
      
      first_sel_index = setdiff(which(outcome_Lasso$beta[,which(outcome_Lasso$lambda==outcome_Lasso_cv$lambda.1se)]!=0),
                                c((dim(fit.cov)[2]-11):(dim(fit.cov)[2])))
      
      # Group-Lasso
      trt_Lasso_cv <- cv.glmnet(x=all_Z, 
                                y=Med, 
                                family = "multinomial",
                                alpha=1,
                                type.multinomial="grouped",
                                nfolds=5)
      trt_Lasso <- glmnet(x=all_Z, 
                          y=Med, 
                          family = "multinomial",
                          alpha=1,
                          type.multinomial="grouped")
      
      second_sel_index = which(trt_Lasso$beta[[1]][,which(trt_Lasso$lambda==trt_Lasso_cv$lambda.1se)]!=0)
      
      single_sel_index = sort(union(first_sel_index,c((dim(fit.cov)[2]-11):(dim(fit.cov)[2]))))
      double_sel_index = sort(union(single_sel_index,second_sel_index))
      
      cv.seed = 123
      
      # Fit model with double selection 
      fit_PDS <- est_binary_laplacian_newton(X=fit.cov[,double_sel_index],
                                          Y=Y,
                                          Lap=0,
                                          lambdal=0,
                                          beta0=rep(0,dim(fit.cov[,double_sel_index])[2]),
                                          maxiter=10000,tol=10^(-6),rate=0.1)
      fit_PDS_Lap_cv <- est_binary_laplacian_cv(est_func = est_binary_laplacian_newton,
                                         X = fit.cov[,double_sel_index],
                                         Y = Y,
                                         Lap = Lap_matrix,
                                         lambda_list =  lambda_grid,
                                         k_folds=5,
                                         this.seed=cv.seed)
      
      fit_PDS_Lap <- est_binary_laplacian_newton(X=fit.cov[,double_sel_index],
                                          Y=Y,
                                          Lap=fit_PDS_Lap_cv$penalty_matrix,
                                          lambdal=fit_PDS_Lap_cv$lambda_grid[which.max(fit_PDS_Lap_cv$performance)],
                                          beta0=rep(0,dim(fit.cov[,double_sel_index])[2]),
                                          maxiter=10000,tol=10^(-6),rate=0.1)
      
      # Fit model - single selection 
      fit_PSS <- est_binary_laplacian_newton(X=fit.cov[,single_sel_index],
                                          Y=Y,
                                          Lap=0,
                                          lambdal=0,
                                          beta0=rep(0,dim(fit.cov[,single_sel_index])[2]),
                                          maxiter=10000,tol=10^(-6),rate=0.1)
      
      fit_PSS_Lap_cv <- est_binary_laplacian_cv(est_func = est_binary_laplacian_newton,
                                         X = fit.cov[,single_sel_index],
                                         Y = Y,
                                         Lap = Lap_matrix,
                                         lambda_list =  lambda_grid,
                                         k_folds=5,
                                         this.seed=cv.seed)
      
      fit_PSS_Lap <- est_binary_laplacian_newton(X=fit.cov[,single_sel_index],
                                          Y=Y,
                                          Lap=fit_PSS_Lap_cv$penalty_matrix,
                                          lambdal=fit_PSS_Lap_cv$lambda_grid[which.max(fit_PSS_Lap_cv$performance)],
                                          beta0=rep(0,dim(fit.cov[,single_sel_index])[2]),
                                          maxiter=10000,tol=10^(-6),rate=0.1)
      
      fit_SLM <- est_binary_laplacian_newton(X=fit.cov,
                                          Y=Y,
                                          Lap=0,
                                          lambdal=0,
                                          beta0=rep(0,dim(fit.cov)[2]),
                                          maxiter=10000,tol=10^(-6),rate=0.1)
      
      write.table(as.data.frame( t(c(as.numeric(args[1]) , resimu,fit_SLM[(length(fit_SLM)-11):length(fit_SLM)]))),
                  file=file.path(save.path,paste0("fitNoReg.txt")),
                  col.names = F,row.names = F,append = T)
      write.table(as.data.frame( t(c(as.numeric(args[1]) , resimu,fit_PSS[(length(fit_PSS)-11):length(fit_PSS)]))),
                  file=file.path(save.path,paste0("fitSingleSelNoReg.txt")),
                  col.names = F,row.names = F,append = T)
      write.table(as.data.frame(t(c(as.numeric(args[1]) , resimu, fit_PSS_Lap[(length(fit_PSS_Lap)-11):length(fit_PSS_Lap)]))),
                  file=file.path(save.path,paste0("fitSingleSelLap.txt")),
                  col.names = F,row.names = F,append = T)
      write.table(as.data.frame(t(c(as.numeric(args[1]) , resimu, fit_PDS[(length(fit_PDS)-11):length(fit_PDS)]))),
                  file=file.path(save.path,paste0("fitDoubleSelNoReg.txt")),
                  col.names = F,row.names = F,append = T)
      write.table(as.data.frame(t(c(as.numeric(args[1]) , resimu, fit_PDS_Lap[(length(fit_PDS_Lap)-11):length(fit_PDS_Lap)]))),
                  file=file.path(save.path,paste0("fitDoubleSelLap.txt")),
                  col.names = F,row.names = F,append = T)
         
    }
  } 
}

set.seed(12345)
interaction_alpha_low <- matrix(round(rnorm(3*2,0,1),3),
                                nrow = 3,byrow = T)
interaction_omega_low <- abs(round(rnorm(2,0,2),3))*c(0.2,1)

set.seed(1234)
interaction_alpha_medium <- matrix(round(rnorm(3*6,0,1),3),
                                   nrow = 3,byrow = T)
interaction_omega_medium <- abs(round(rnorm(6,0,1.5),3))*rep(c(1,0.2),3)
interaction_alpha_medium <- cbind(interaction_alpha_medium[,1],
                                  interaction_alpha_low[,1],
                                  interaction_alpha_medium[,2:5],
                                  interaction_alpha_low[,2],
                                  interaction_alpha_medium[,6])
interaction_omega_medium <- c(interaction_omega_medium[1],
                              interaction_omega_low[1],
                              interaction_omega_medium[2:5],
                              interaction_omega_low[2],
                              interaction_omega_medium[6])

set.seed(1)
# Scenario i interaction (i) + 2 interactions
high_int = F
alphas1 <- alphas0
omega1 <- omega
alphas1[,c(7,10,14,21)] <- 0 
alphas1 <- cbind(c(-2.95,-2.15,-1.9),alphas1,interaction_alpha_low)
omega1[c(7,10,14,21)] <- 0
omega1 <- c(omega1,interaction_omega_low)

# Scenario j interaction (ii) + 8 interactions
high_int = T
alphas2 <- alphas0
omega2 <- omega
alphas2[,c(7,10,14,21)] <- 0 
alphas2 <- cbind(c(-3.06,-2.2,-1.9),alphas2,interaction_alpha_medium)
omega2[c(7,10,14,21)] <- 0
omega2 <- c(omega2,interaction_omega_medium)


prevalence=c(4,8,15)
scenarios=c("i","j")
alphas=list(alphas1,alphas2)
omegas=list(omega1,omega2)
intercept_main=c(-2.87,-3.1, 
                 -1.9,-2.04,
                 -0.93,-1.01)
high_inti=list(FALSE,TRUE)

for(l in 1:3){
  for(j in 1:2){
    for(k in c(0,50)){
      simu_results=NULL
      
      resimu=k
      seedi <- k+l*100+j*1000 + 10000*as.numeric(args[1])+resimu*100000
      cat("args",seedi, "\n")

      set.seed(seed=seedi)
      simu_datai=simu_data_inter(alphas[[j]],omegas[[j]],intercept_main[(l-1)*2+j],high_inti[[j]],gammas,4500)
      dim(simu_datai)
      Y=simu_datai[,1]
      
      Z <- cbind(simu_datai[,3:24],
                 simu_datai[,3]*simu_datai[,c(8:24)])
      X=simu_datai[,25:36]
      C=simu_datai[,37]
      Med=apply(X,1,function(x) which(x==1))
      Z  <- as.matrix(Z)
      X <- as.matrix(X)
      fit.cov = as.matrix(cbind(Z,X))
      all_Z = as.matrix(Z)
      colnames(all_Z) = paste("X",1:(dim(all_Z)[2]),sep="")
      colnames(fit.cov) = c(colnames(all_Z),paste("d",1:12,sep=""))
      
      dim(fit.cov) # 4500x68
      
      C=as.factor(C)
      Med=as.factor(Med)
      
      save.path=paste0("~/ZWang-Lab/result_",prevalence[l],scenarios[j])
      
      # LASSO
      outcome_Lasso_cv <- cv.glmnet(x=fit.cov,
                                    y=Y,
                                    family="binomial",
                                    intercept =F,
                                    type.measure = "auc",
                                    nfolds=5,
                                    penalty.factor=c(rep(1,dim(fit.cov)[2]-12),rep(0,12)))
      
      outcome_Lasso <- glmnet(x=fit.cov,
                              y=Y,
                              family="binomial",
                              intercept =F,
                              penalty.factor=c(rep(1,dim(fit.cov)[2]-12),rep(0,12)))
      
      first_sel_index = setdiff(which(outcome_Lasso$beta[,which(outcome_Lasso$lambda==outcome_Lasso_cv$lambda.1se)]!=0),
                                c((dim(fit.cov)[2]-11):(dim(fit.cov)[2])))
      
      # Group-LASSO
      trt_Lasso_cv <- cv.glmnet(x=all_Z, 
                                y=Med, 
                                family = "multinomial",
                                alpha=1,
                                type.multinomial="grouped",
                                nfolds=5)
      trt_Lasso <- glmnet(x=all_Z, 
                          y=Med, 
                          family = "multinomial",
                          alpha=1,
                          type.multinomial="grouped")
      
      second_sel_index = which(trt_Lasso$beta[[1]][,which(trt_Lasso$lambda==trt_Lasso_cv$lambda.1se)]!=0)
      
      single_sel_index = sort(union(first_sel_index,c((dim(fit.cov)[2]-11):(dim(fit.cov)[2]))))
      double_sel_index = sort(union(single_sel_index,second_sel_index))

      cv.seed = 123
      # Fit model with double selection 
      fit_PDS <- est_binary_laplacian_newton(X=fit.cov[,double_sel_index],
                                          Y=Y,
                                          Lap=0,
                                          lambdal=0,
                                          beta0=rep(0,dim(fit.cov[,double_sel_index])[2]),
                                          maxiter=10000,tol=10^(-6),rate=0.1)
      
      fit_PDS_Lap_cv <- est_binary_laplacian_cv(est_func = est_binary_laplacian_newton,
                                         X = fit.cov[,double_sel_index],
                                         Y = Y,
                                         Lap = Lap_matrix,
                                         lambda_list =  lambda_grid,
                                         k_folds=5,
                                         this.seed=cv.seed)

      fit_PDS_Lap <- est_binary_laplacian_newton(X=fit.cov[,double_sel_index],
                                          Y=Y,
                                          Lap=fit_PDS_Lap_cv$penalty_matrix,
                                          lambdal=fit_PDS_Lap_cv$lambda_grid[which.max(fit_PDS_Lap_cv$performance)],
                                          beta0=rep(0,dim(fit.cov[,double_sel_index])[2]),
                                          maxiter=10000,tol=10^(-6),rate=0.1)

      
      # Fit model - single selection 
      fit_PSS <- est_binary_laplacian_newton(X=fit.cov[,single_sel_index],
                                          Y=Y,
                                          Lap=0,
                                          lambdal=0,
                                          beta0=rep(0,dim(fit.cov[,single_sel_index])[2]),
                                          maxiter=10000,tol=10^(-6),rate=0.1)

      fit_PSS_Lap_cv <- est_binary_laplacian_cv(est_func = est_binary_laplacian_newton,
                                         X = fit.cov[,single_sel_index],
                                         Y = Y,
                                         Lap = Lap_matrix,
                                         lambda_list =  lambda_grid,
                                         k_folds=5,
                                         this.seed=cv.seed)

      fit_PSS_Lap <- est_binary_laplacian_newton(X=fit.cov[,single_sel_index],
                                          Y=Y,
                                          Lap=fit_PSS_Lap_cv$penalty_matrix,
                                          lambdal=fit_PSS_Lap_cv$lambda_grid[which.max(fit_PSS_Lap_cv$performance)],
                                          beta0=rep(0,dim(fit.cov[,single_sel_index])[2]),
                                          maxiter=10000,tol=10^(-6),rate=0.1)

      fit_SLM <- est_binary_laplacian_newton(X=fit.cov,
                                          Y=Y,
                                          Lap=0,
                                          lambdal=0,
                                          beta0=rep(0,dim(fit.cov)[2]),
                                          maxiter=10000,tol=10^(-6),rate=0.1)
      
      write.table(as.data.frame( t(c(as.numeric(args[1]) , resimu,fit_SLM[(length(fit_SLM)-11):length(fit_SLM)]))),
                  file=file.path(save.path,paste0("fitNoReg.txt")),
                  col.names = F,row.names = F,append = T)
      write.table(as.data.frame( t(c(as.numeric(args[1]) , resimu,fit_PSS[(length(fit_PSS)-11):length(fit_PSS)]))),
                  file=file.path(save.path,paste0("fitSingleSelNoReg.txt")),
                  col.names = F,row.names = F,append = T)
      write.table(as.data.frame(t(c(as.numeric(args[1]) , resimu, fit_PSS_Lap[(length(fit_PSS_Lap)-11):length(fit_PSS_Lap)]))),
                  file=file.path(save.path,paste0("fitSingleSelLap.txt")),
                  col.names = F,row.names = F,append = T)
      write.table(as.data.frame(t(c(as.numeric(args[1]) , resimu, fit_PDS[(length(fit_PDS)-11):length(fit_PDS)]))),
                  file=file.path(save.path,paste0("fitDoubleSelNoReg.txt")),
                  col.names = F,row.names = F,append = T)
      write.table(as.data.frame(t(c(as.numeric(args[1]) , resimu, fit_PDS_Lap[(length(fit_PDS_Lap)-11):length(fit_PDS_Lap)]))),
                  file=file.path(save.path,paste0("fitDoubleSelLap.txt")),
                  col.names = F,row.names = F,append = T)
      
    }
  } 
}


save.path = ("~/ZWang-Lab")

pairwise_diff_unique <- function(a) {
  # Compute differences for all unique pairs
  differences <- combn(a, 2, FUN = function(x) x[1] - x[2])
  return(differences)
}

pairwise_diff_unique_name <- function(a) {
  # Compute differences for all unique pairs
  differences <- combn(a, 2, FUN = function(x) paste0(x[1], "v", x[2]))
  return(differences)
}

HR = pairwise_diff_unique_name(c(1:12))
true_CHR=pairwise_diff_unique(betas)

vs = c("1v1", "1v1", "1v1", "1v1", "1v2", "1v2", "1v3", "1v3", "1v4", "1v4", "1v4", 
       "1v1", "1v1", "1v1", "1v2", "1v2", "1v3", "1v3", "1v4", "1v4", "1v4",
       "1v1", "1v1", "1v2", "1v2", "1v3", "1v3", "1v4", "1v4", "1v4", 
       "1v1", "1v2", "1v2", "1v3", "1v3", "1v4", "1v4", "1v4", 
       "1v2", "1v2", "1v3", "1v3", "1v4", "1v4", "1v4", 
       "2v2", "2v3", "2v3", "2v4", "2v4", "2v4", 
       "2v3", "2v3", "2v4", "2v4", "2v4", 
       "3v3", "3v4", "3v4", "3v4", "3v4", "3v4", "3v4",
       "4v4", "4v4", "4v4")
	   
evaluation = data.frame(lambda=NA,
                        BIAS=NA,
                        BIAS_rel=NA,
                        MSE=NA,
                        Scenario=NA,
                        HR=NA,
                        vs=NA,
                        Prevalence=NA)
						
for(k_index in c(4,8,15)){
  for(s_index in c("a","b","c","d","e","f","g","h","i","j")){
    file.names = list.files(path=file.path(save.path,paste0("result_",k_index,s_index)), pattern="fit")

    result = list()
    for(s in 1:5){
      df = read.table(file.path(save.path,paste0("result_",k_index,s_index),file.names[s]),header = F)[,-c(1,2)]
      this.name = gsub("fit|.txt","",file.names[s])
      final_df=t(apply(df,1,pairwise_diff_unique))

      tmp = data.frame(lambda = this.name,
                       BIAS = colMeans(sapply(1:66, function(j) final_df[,j]-true_CHR[j])), 
                       BIAS_rel = colMeans(sapply(1:66, function(j) final_df[,j]-true_CHR[j]))/true_CHR, 
                       MSE = colMeans(sapply(1:66, function(j) (final_df[,j]-true_CHR[j])^2)),
                       Scenario=s_index,
                       HR=HR,
                       vs=vs,
                       Prevalence=paste0(k_index,"%"))
      evaluation = rbind(evaluation, tmp)
      print(paste0(k_index,s_index," - ",this.name,": ",nrow(df)))
    }
  }
}

evaluation = evaluation[-1,]
evaluation$lambda = factor(evaluation$lambda, 
                           levels = c("NoReg", 
                                      "SingleSelNoReg", "DoubleSelNoReg",
                                      "SingleSelLap","DoubleSelLap"),
                           labels = c("SLR", 
                                      "PSS", "PDS",
                                      "PSS-Lap","PDS-Lap"))

evaluation$Prevalence = factor(evaluation$Prevalence,
                         levels = c("4%","8%","15%"))

evaluation$Scenario = toupper(evaluation$Scenario)

evaluation$group_comp = ifelse(evaluation$vs%in%c("1v1","2v2","3v3","4v4"),
                               "Within group","Across group")

p_MSE = ggplot(evaluation, 
       aes(x=Scenario,y=MSE, fill = lambda)) +
  stat_boxplot(geom = "errorbar", width = 0.75) +
  geom_boxplot(outliers = F) +
  theme_classic()+
  theme(legend.position = "bottom",
        legend.box.spacing = unit(0, 'cm'),
        legend.text = element_text(size=10,color="black"),
        legend.title = element_blank(),
        axis.text = element_text(size=10,face="bold",color="black"),
        axis.title = element_text(size=14,color="black"),
        strip.text = element_text(color="black",family = "sans"))+
  ylab("MSE")+
  facet_grid(rows=vars(Prevalence), 
             scales="free",space="free_x",
             axes="all_x",
             labeller = labeller(Prevalence = label_both))+
  ggtitle(NULL)+
  scale_fill_manual(values = c("#FFF08C","#F28CA4","#A1B1EB", "#e6194b", "#4363d8"))

p_Bias = ggplot(evaluation, aes(x=Scenario,y=abs(BIAS),
                                fill = lambda)) +
  stat_boxplot(geom = "errorbar", width = 0.75) +
  geom_boxplot(outliers = F) +
  theme_classic()+
  theme(legend.position = "bottom",
        legend.box.spacing = unit(0, 'cm'),
        legend.text = element_text(size=10,color="black"),
        legend.title = element_blank(),
        axis.text = element_text(size=10,face="bold",color="black"),
        axis.title = element_text(size=14,color="black"),
        strip.text = element_text(color="black",family = "sans"))+
  ylab("Bias")+
  scale_y_log10()+
  facet_grid(rows=vars(Prevalence), 
             scales="free",space="free_x",
             axes="all_x",
             labeller = labeller(Prevalence = label_both))+
  ggtitle(NULL)+
  scale_fill_manual(values = c("#FFF08C","#F28CA4","#A1B1EB", "#e6194b", "#4363d8"))


evaluation$group_comp = ifelse(evaluation$group_comp=="Within group","Within drug class","Between drug class")
evaluation$group_comp = factor(evaluation$group_comp,
                               levels=c("Within drug class","Between drug class"))
							   
p_BIAS_group = ggplot(evaluation, aes(x=Scenario,y=abs(BIAS),
                                      fill = lambda)) +
  stat_boxplot(geom = "errorbar", width = 0.75) +
  geom_boxplot(outliers = F) +
  theme_classic()+
  theme(legend.position = "bottom",
        legend.box.spacing = unit(0, 'cm'),
        legend.text = element_text(size=10,color="black"),
        legend.title = element_blank(),
        axis.text = element_text(size=10,face="bold",color="black"),
        axis.title = element_text(size=14,color="black"),
        strip.text = element_text(color="black",family = "sans"))+
  ylab("Bias")+
  scale_y_log10()+
  facet_grid(rows=vars(Prevalence), 
             cols=vars(group_comp),
             scales="free",space="free_x",
             axes="all_x",
             labeller = labeller(Prevalence = label_both))+
  ggtitle(NULL)+
  scale_fill_manual(values = c("#FFF08C","#F28CA4","#A1B1EB", "#e6194b", "#4363d8"))
