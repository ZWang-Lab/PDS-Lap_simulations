# Generate data w/o interaction terms
simu_data=function(alphasi,omegai,intercept_outcome,gammas,sim.n = 4500){

  X1 = rep(1, sim.n) #intercept
  X2 = rnorm(sim.n,0,1) #age
  X3 = rbinom(sim.n,1,0.35) #sex
  X4 = t(rmultinom(sim.n,1,c(0.9,0.04,0.04,0.02))) #race
  
  # comorbidities 
  C1 = matrix(rbinom(sim.n*0.30*1,1,0.5),nrow=sim.n*0.30,ncol=1)
  C2 = matrix(rbinom(sim.n*0.30*1,1,0.2),nrow=sim.n*0.30,ncol=1)
  C3 = matrix(rbinom(sim.n*0.30*3,1,0.15),nrow=sim.n*0.30,ncol=3)
  C4 = matrix(rbinom(sim.n*0.30*3,1,0.1),nrow=sim.n*0.30,ncol=3)
  C5 = matrix(rbinom(sim.n*0.30*4,1,0.05),nrow=sim.n*0.30,ncol=4)
  C6 = matrix(rbinom(sim.n*0.30*1,1,0.01),nrow=sim.n*0.30,ncol=1)
  
  comorbidities = cbind(C1,C2,C3,C4,C5,C6)
  comorbidities = rbind(comorbidities, matrix(0,ncol=13,nrow=sim.n*0.70))
  comorbidities = comorbidities[sample(1:sim.n,sim.n),]
  
  C7 = matrix(rbinom(sim.n*2,1,0.2),nrow=sim.n,ncol=2)
  C8 = matrix(rbinom(sim.n*1,1,0.15),nrow=sim.n,ncol=1)
  C9 = matrix(rbinom(sim.n*1,1,0.02),nrow=sim.n,ncol=1)
  
  comorbidities = cbind(C7,C8,comorbidities,C9)
  covariates <- cbind(X1, X2,X3,X4[,-1],comorbidities)

  # simulate treatment (class)
  prob.base <- 1/(rowSums(exp(covariates%*%t(alphasi)))+1)
  prob.trt <- cbind(prob.base, (exp(as.matrix(covariates)%*%t(alphasi)))*(prob.base))
  
  trt.class.assign.dummy <- t(apply(prob.trt,1, function(x) rmultinom(1,1,x)))
  colMeans(trt.class.assign.dummy)
  
  trt.class.assign.cat <- max.col(trt.class.assign.dummy)
  table(trt.class.assign.cat)/sim.n
  
  # simulate treatment (drug)
  trt.drug.assign.dummy <- t(sapply(trt.class.assign.cat,function(x) rmultinom(1,1,gammas[x,])))
  trt.drug.assign.cat <- max.col(trt.drug.assign.dummy)

  # simulate outcome
  tmp = exp(trt.drug.assign.dummy%*%betas+covariates%*%c(intercept_outcome,omegai))
  y = sapply(tmp/(1+tmp),function(x) rbinom(1,1,x))
  
  while((sum(table(y,trt.drug.assign.cat)[2,]==0)+sum(sapply(3:ncol(covariates), 
                                                             function(i) table(y,covariates[,i])[2,2])<=1))>0){
    # simulate treatment (class)
    prob.base <- 1/(rowSums(exp(covariates%*%t(alphasi)))+1)
    prob.trt <- cbind(prob.base, (exp(as.matrix(covariates)%*%t(alphasi)))*(prob.base))
    trt.class.assign.dummy <- t(apply(prob.trt,1, function(x) rmultinom(1,1,x)))
    colMeans(trt.class.assign.dummy)
    
    trt.class.assign.cat <- max.col(trt.class.assign.dummy)
    table(trt.class.assign.cat)/sim.n
    
    # simulate treatment (drug)
    trt.drug.assign.dummy <- t(sapply(trt.class.assign.cat,function(x) rmultinom(1,1,gammas[x,])))
    trt.drug.assign.cat <- max.col(trt.drug.assign.dummy)

    # simulate outcome (survival)
    tmp = exp(trt.drug.assign.dummy%*%betas+covariates%*%c(intercept_outcome,omegai))
    y = sapply(tmp/(1+tmp),function(x) rbinom(1,1,x))
    
  }

  # Save simulation data 
  sim.data <- data.frame(cbind(y,covariates,
                               trt.drug.assign.dummy,trt.class.assign.cat))
  colnames(sim.data) <- c("y", paste("X",c(1:ncol(covariates)),sep = ""), 
                          paste("d",c(1:12),sep = ""),"Class")
  
  return(sim.data)
}

# Generate data w/ interaction terms
simu_data_inter=function(alphasi,omegai,intercept_outcome,high_inti,gammas,sim.n = 4500){
  y=rep(0,12*sim.n)
  # simulate covariates 
  X1 = rep(1, sim.n) #intercept
  X2 = rnorm(sim.n,0,1) #age
  X3 = rbinom(sim.n,1,0.35) #sex
  X4 = t(rmultinom(sim.n,1,c(0.9,0.04,0.04,0.02))) #race
  
  # comorbidities 
  C1 = matrix(rbinom(sim.n*0.30*1,1,0.5),nrow=sim.n*0.30,ncol=1)
  C2 = matrix(rbinom(sim.n*0.30*1,1,0.2),nrow=sim.n*0.30,ncol=1)
  C3 = matrix(rbinom(sim.n*0.30*3,1,0.15),nrow=sim.n*0.30,ncol=3)
  C4 = matrix(rbinom(sim.n*0.30*3,1,0.1),nrow=sim.n*0.30,ncol=3)
  C5 = matrix(rbinom(sim.n*0.30*4,1,0.05),nrow=sim.n*0.30,ncol=4)
  C6 = matrix(rbinom(sim.n*0.30*1,1,0.01),nrow=sim.n*0.30,ncol=1)
  
  comorbidities = cbind(C1,C2,C3,C4,C5,C6)
  comorbidities = rbind(comorbidities, matrix(0,ncol=13,nrow=sim.n*0.70))
  comorbidities = comorbidities[sample(1:sim.n,sim.n),]
  
  C7 = matrix(rbinom(sim.n*2,1,0.2),nrow=sim.n,ncol=2)
  C8 = matrix(rbinom(sim.n*1,1,0.15),nrow=sim.n,ncol=1)
  C9 = matrix(rbinom(sim.n*1,1,0.02),nrow=sim.n,ncol=1)
  
  comorbidities = cbind(C7,C8,comorbidities,C9)
  
  if(high_inti){
    interaction_terms <- X2*comorbidities[,c(1,4,7,8,10,13,14,17)]
  }else{
    interaction_terms <- cbind(X2*comorbidities[,4],X2*comorbidities[,14])
  }

  covariates <- cbind(X1, X2,X3,X4[,-1],comorbidities, interaction_terms)
  
  # simulate treatment (class)
  prob.base <- 1/(rowSums(exp(covariates%*%t(alphasi)))+1)
  prob.trt <- cbind(prob.base, (exp(as.matrix(covariates)%*%t(alphasi)))*(prob.base))
  
  trt.class.assign.dummy <- t(apply(prob.trt,1, function(x) rmultinom(1,1,x)))
  colMeans(trt.class.assign.dummy)
  
  trt.class.assign.cat <- max.col(trt.class.assign.dummy)
  table(trt.class.assign.cat)/sim.n
  
  # simulate treatment (drug)
  trt.drug.assign.dummy <- t(sapply(trt.class.assign.cat,function(x) rmultinom(1,1,gammas[x,])))
  trt.drug.assign.cat <- max.col(trt.drug.assign.dummy)
  
  # simulate outcome (survival)
  tmp = exp(trt.drug.assign.dummy%*%betas+covariates%*%c(intercept_outcome,omegai))
  y = sapply(tmp/(1+tmp),function(x) rbinom(1,1,x))
  
  while((sum(table(y,trt.drug.assign.cat)[2,]==0)+sum(sapply(3:ncol(cbind(X1, X2,X3,X4[,-1],comorbidities)), 
                                                             function(i) table(y,covariates[,i])[2,2])<=1))>0){
    # simulate treatment (class)
    prob.base <- 1/(rowSums(exp(covariates%*%t(alphasi)))+1)
    prob.trt <- cbind(prob.base, (exp(as.matrix(covariates)%*%t(alphasi)))*(prob.base))
    trt.class.assign.dummy <- t(apply(prob.trt,1, function(x) rmultinom(1,1,x)))
    colMeans(trt.class.assign.dummy)
    
    trt.class.assign.cat <- max.col(trt.class.assign.dummy)
    table(trt.class.assign.cat)/sim.n
    
    # simulate treatment (drug)
    trt.drug.assign.dummy <- t(sapply(trt.class.assign.cat,function(x) rmultinom(1,1,gammas[x,])))
    trt.drug.assign.cat <- max.col(trt.drug.assign.dummy)
    
    # simulate outcome (survival)
    tmp = exp(trt.drug.assign.dummy%*%betas+covariates%*%c(intercept_outcome,omegai))
    y = sapply(tmp/(1+tmp),function(x) rbinom(1,1,x))
  }

  # Save simulation data 
  sim.data <- data.frame(cbind(y,cbind(X1, X2,X3,X4[,-1],comorbidities),
                               trt.drug.assign.dummy,trt.class.assign.cat))
  colnames(sim.data) <- c("y", paste("X",c(1:ncol(cbind(X1, X2,X3,X4[,-1],comorbidities))),sep = ""), 
                          paste("d",c(1:12),sep = ""),"Class")
  
  return(sim.data)
  
}