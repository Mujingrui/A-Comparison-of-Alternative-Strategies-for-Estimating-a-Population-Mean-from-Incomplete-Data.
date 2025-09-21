library(MASS)
N <- 200
dim_z <- 4
R_z <- matrix(c(1,0,0,0,
                0,1,0,0,
                0,0,1,0,
                0,0,0,1), 
                nrow = 4, ncol = 4, byrow = TRUE)
Z <- mvrnorm(n = N, mu = rep(0,dim_z), Sigma = R_z)
beta <- c(210,27.4,13.7,13.7,13.7)
Y <- (cbind(rep(1,N),Z))%*%beta

expit <- function(x){
  return(1/(1+exp(-x)))
}

true_prs <- rep(0,N)
for(i in 1:N){
  true_prs[i] <- expit(Z[i,]%*%c(-1,0.5,-0.25,-0.1))
}
  
X <- matrix(data = NA, nrow = N, ncol = 4)
for(i in 1:N){
  X[i,1] <- exp(Z[i,1]/2)
  X[i,2] <- Z[i,2]/(1+exp(Z[i,1])) + 10
  X[i,3] <- ((Z[i,1]*Z[i,3])/25 + 0.6)^3
  X[i,4] <- (Z[i,2] + Z[i,4] + 20)^2
}

# respondents_index <- which(true_prs >= 0.5)
rs_index <- t <- rbinom(N,1,true_prs)
respondents_index <- which(rs_index == 1)
Y_respondents <- Y[respondents_index]
X_observed <- X[respondents_index,]

par(mfrow = c(2,2)) 
plot(X_observed[,1],Y_respondents, 
     main="",
     xlab = "x1",
     ylab = "y",
     col=rgb(0,100,0,50,maxColorValue=255), 
     pch=16)
plot(X_observed[,2],Y_respondents, 
     main="",
     xlab = "x2",
     ylab = "y",
     col=rgb(0,100,0,50,maxColorValue=255), 
     pch=16)
plot(X_observed[,3],Y_respondents, 
     main="",
     xlab = "x3",
     ylab = "y",
     col=rgb(0,100,0,50,maxColorValue=255), 
     pch=16)
plot(X_observed[,4],Y_respondents, 
     main="",
     xlab = "x4",
     ylab = "y",
     col=rgb(0,100,0,50,maxColorValue=255), 
     pch=16)

library(ggplot2)

x1 <- c(X_observed[,1],X[-respondents_index,1])
observed_x1 <- c(rep("1",97),rep("0",103))
dt_x1 <- data.frame(x1 = x1, observed = observed_x1)
ggplot(dt_x1, aes(x=observed_x1, y=x1)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
   # coord_cartesian(ylim=c(0.5,3.5))+
  xlab("t")+
  ylab("x1")

x2 <- c(X_observed[,2],X[-respondents_index,2])
observed_x2 <- c(rep("1",97),rep("0",103))
dt_x2 <- data.frame(x2 = x2, observed = observed_x2)
ggplot(dt_x2, aes(x=observed_x2, y=x2)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  coord_cartesian(ylim=c(8.5,11.5))+
  xlab("t")+
  ylab("x2")

x3 <- c(X_observed[,3],X[-respondents_index,3])
observed_x3 <- c(rep("1",97),rep("0",103))
dt_x3 <- data.frame(x3 = x3, observed = observed_x3)
ggplot(dt_x3, aes(x=observed_x3, y=x3)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  coord_cartesian(ylim=c(0.1,0.4))+
  xlab("t")+
  ylab("x3")

x4 <- c(X_observed[,4],X[-respondents_index,4])
observed_x4 <- c(rep("1",97),rep("0",103))
dt_x4 <- data.frame(x4 = x4, observed = observed_x4)
ggplot(dt_x4, aes(x=observed_x4, y=x4)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  coord_cartesian(ylim=c(250,600))+
  xlab("t")+
  ylab("x4")

# t <- rep(0, N)
# for(i in 1:N){
#   if (i %in% respondents_index){
#     t[i] <- 1
#   }
# }
# X1 = X[,1], X2 = X[,2], X3 = X[,3], X4 = X[,4]
data_XTY <- data.frame(X1 = X[,1], X2 = X[,2], X3 = X[,3], X4 = X[,4], Y = Y, t = t)

glm_tx <- glm(t~X1+X2+X3+X4, data = data_XTY, family = "binomial")

eta_fitted <- cbind(rep(1,N),X) %*% coefficients(glm_tx)

dt_eta <- data.frame(eta_fitted = eta_fitted, t = t)
ggplot(dt_eta, aes(x=as.factor(t), y=eta_fitted)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  # coord_cartesian(ylim=c(250,600))+
  xlab("t")+
  ylab("logit propensity score")

nreps<-1000 
ests<-matrix(0,nrow=nreps,ncol=4)
for(rep in 1:nreps){
  N <- 200
  dim_z <- 4
  R_z <- matrix(c(1,0,0,0,
                  0,1,0,0,
                  0,0,1,0,
                  0,0,0,1), 
                nrow = 4, ncol = 4, byrow = TRUE)
  Z <- mvrnorm(n = N, mu = rep(0,dim_z), Sigma = R_z)
  beta <- c(210,27.4,13.7,13.7,13.7)
  Y <- (cbind(rep(1,N),Z))%*%beta
  true_prs <- rep(0,N)
  for(i in 1:N){
    true_prs[i] <- expit(Z[i,]%*%c(-1,0.5,-0.25,-0.1))
  }
  X <- matrix(data = NA, nrow = N, ncol = 4)
  for(i in 1:N){
    X[i,1] <- exp(Z[i,1]/2)
    X[i,2] <- Z[i,2]/(1+exp(Z[i,1])) + 10
    X[i,3] <- ((Z[i,1]*Z[i,3])/25 + 0.6)^3
    X[i,4] <- (Z[i,2] + Z[i,4] + 20)^2
  }
  rs_index <- t <- rbinom(N,1,true_prs)
  respondents_index <- which(rs_index == 1)
  Y_respondents <- Y[respondents_index]
  X_observed <- X[respondents_index,]
  
  dt_XTYZ <- data.frame(Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], Z4 = Z[,4], X1 = X[,1], X2 = X[,2], X3 = X[,3], X4 = X[,4], Y = Y, t = t)
  eX_specified <- fitted(glm(t~Z1+Z2+Z3+Z4,data = dt_XTYZ, family=binomial))
  eX_misspecified <- fitted(glm(t~X1+X2+X3+X4,data = dt_XTYZ, family=binomial))
  
  w1<-t/eX_specified
  W1<-w1/sum(w1)
  w0<-(1-t)/(1-eX_specified)
  W0<-w0/sum(w0)
  r_1 <- length(respondents_index)/N
  ests[rep,1]<-sum(W1*Y)
  ests[rep,2]<-(1-r_1)*sum(W0*Y)+r_1*sum(W1*Y)
  
  w1_mis<-t/eX_misspecified
  W1_mis<-w1_mis/sum(w1_mis)
  w0_mis<-(1-t)/(1-eX_misspecified)
  W0_mis<-w0_mis/sum(w0_mis)
  
  ests[rep,3]<-sum(W1_mis*Y)
  ests[rep,4]<-(1-r_1)*sum(W0_mis*Y)+r_1*sum(W1_mis*Y)
  
}


nreps<-1000 
ests<-matrix(0,nrow=nreps,ncol=6)
for(rep in 1:nreps){
  N <- 200
  dim_z <- 4
  R_z <- matrix(c(1,0,0,0,
                  0,1,0,0,
                  0,0,1,0,
                  0,0,0,1), 
                nrow = 4, ncol = 4, byrow = TRUE)
  Z <- mvrnorm(n = N, mu = rep(0,dim_z), Sigma = R_z)
  beta <- c(210,27.4,13.7,13.7,13.7)
  Y <- (cbind(rep(1,N),Z))%*%beta
  true_prs <- rep(0,N)
  for(i in 1:N){
    true_prs[i] <- expit(Z[i,]%*%c(-1,0.5,-0.25,-0.1))
  }
  X <- matrix(data = NA, nrow = N, ncol = 4)
  for(i in 1:N){
    X[i,1] <- exp(Z[i,1]/2)
    X[i,2] <- Z[i,2]/(1+exp(Z[i,1])) + 10
    X[i,3] <- ((Z[i,1]*Z[i,3])/25 + 0.6)^3
    X[i,4] <- (Z[i,2] + Z[i,4] + 20)^2
  }
  rs_index <- t <- rbinom(N,1,true_prs)
  respondents_index <- which(rs_index == 1)
  Y_respondents <- Y[respondents_index]
  X_observed <- X[respondents_index,]
  
  dt_XTYZ <- data.frame(Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], Z4 = Z[,4], X1 = X[,1], X2 = X[,2], X3 = X[,3], X4 = X[,4], Y = Y, t = t)
  dt_XTYZ_res <- dt_XTYZ[respondents_index,]
  eX_specified <- fitted(glm(t~Z1+Z2+Z3+Z4,data = dt_XTYZ, family=binomial))
  eX_misspecified <- fitted(glm(t~X1+X2+X3+X4,data = dt_XTYZ, family=binomial))
  
  w1<-t/eX_specified
  W1<-w1/sum(w1)
  w0<-(1-t)/(1-eX_specified)
  W0<-w0/sum(w0)
  r_1 <- length(respondents_index)/N
  ests[rep,1]<-sum(W1*Y)
  ests[rep,2]<-(1-r_1)*sum(W0*Y)+r_1*sum(W1*Y)
  
  w1_mis<-t/eX_misspecified
  W1_mis<-w1_mis/sum(w1_mis)
  w0_mis<-(1-t)/(1-eX_misspecified)
  W0_mis<-w0_mis/sum(w0_mis)
  
  ests[rep,3]<-sum(W1_mis*Y)
  ests[rep,4]<-(1-r_1)*sum(W0_mis*Y)+r_1*sum(W1_mis*Y)
  
  # Y_specified <- fitted(lm(Y~Z1+Z2+Z3+Z4,data = dt_XTYZ_res))
  ests[rep,5] <- sum(coefficients(lm(Y~Z1+Z2+Z3+Z4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$Z1,dt_XTYZ$Z2,dt_XTYZ$Z3,dt_XTYZ$Z4)))/N
  
  # Y_misspecified <- fitted(lm(Y~X1+X2+X3+X4,data = dt_XTYZ_res))
  ests[rep,6] <- sum(coefficients(lm(Y~X1+X2+X3+X4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$X1,dt_XTYZ$X2,dt_XTYZ$X3,dt_XTYZ$X4)))/N
}

nreps<-1000 
ests_1000<-matrix(0,nrow=nreps,ncol=6)
for(rep in 1:nreps){
  N <- 1000
  dim_z <- 4
  R_z <- matrix(c(1,0,0,0,
                  0,1,0,0,
                  0,0,1,0,
                  0,0,0,1), 
                nrow = 4, ncol = 4, byrow = TRUE)
  Z <- mvrnorm(n = N, mu = rep(0,dim_z), Sigma = R_z)
  beta <- c(210,27.4,13.7,13.7,13.7)
  Y <- (cbind(rep(1,N),Z))%*%beta
  true_prs <- rep(0,N)
  for(i in 1:N){
    true_prs[i] <- expit(Z[i,]%*%c(-1,0.5,-0.25,-0.1))
  }
  X <- matrix(data = NA, nrow = N, ncol = 4)
  for(i in 1:N){
    X[i,1] <- exp(Z[i,1]/2)
    X[i,2] <- Z[i,2]/(1+exp(Z[i,1])) + 10
    X[i,3] <- ((Z[i,1]*Z[i,3])/25 + 0.6)^3
    X[i,4] <- (Z[i,2] + Z[i,4] + 20)^2
  }
  rs_index <- t <- rbinom(N,1,true_prs)
  respondents_index <- which(rs_index == 1)
  Y_respondents <- Y[respondents_index]
  X_observed <- X[respondents_index,]
  
  dt_XTYZ <- data.frame(Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], Z4 = Z[,4], X1 = X[,1], X2 = X[,2], X3 = X[,3], X4 = X[,4], Y = Y, t = t)
  dt_XTYZ_res <- dt_XTYZ[respondents_index,]
  eX_specified <- fitted(glm(t~Z1+Z2+Z3+Z4,data = dt_XTYZ, family=binomial))
  eX_misspecified <- fitted(glm(t~X1+X2+X3+X4,data = dt_XTYZ, family=binomial))
  
  w1<-t/eX_specified
  W1<-w1/sum(w1)
  w0<-(1-t)/(1-eX_specified)
  W0<-w0/sum(w0)
  r_1 <- length(respondents_index)/N
  ests_1000[rep,1]<-sum(W1*Y)
  ests_1000[rep,2]<-(1-r_1)*sum(W0*Y)+r_1*sum(W1*Y)
  
  w1_mis<-t/eX_misspecified
  W1_mis<-w1_mis/sum(w1_mis)
  w0_mis<-(1-t)/(1-eX_misspecified)
  W0_mis<-w0_mis/sum(w0_mis)
  
  ests_1000[rep,3]<-sum(W1_mis*Y)
  ests_1000[rep,4]<-(1-r_1)*sum(W0_mis*Y)+r_1*sum(W1_mis*Y)
  
  # Y_specified <- fitted(lm(Y~Z1+Z2+Z3+Z4,data = dt_XTYZ_res))
  ests_1000[rep,5] <- sum(coefficients(lm(Y~Z1+Z2+Z3+Z4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$Z1,dt_XTYZ$Z2,dt_XTYZ$Z3,dt_XTYZ$Z4)))/N
  
  # Y_misspecified <- fitted(lm(Y~X1+X2+X3+X4,data = dt_XTYZ_res))
  ests_1000[rep,6] <- sum(coefficients(lm(Y~X1+X2+X3+X4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$X1,dt_XTYZ$X2,dt_XTYZ$X3,dt_XTYZ$X4)))/N
}

nreps_AIPW<-1000 
ests_AIPW<-matrix(0,nrow=nreps,ncol=4)
for(rep in 1:nreps){
  N <- 200
  dim_z <- 4
  R_z <- matrix(c(1,0,0,0,
                  0,1,0,0,
                  0,0,1,0,
                  0,0,0,1), 
                nrow = 4, ncol = 4, byrow = TRUE)
  Z <- mvrnorm(n = N, mu = rep(0,dim_z), Sigma = R_z)
  beta <- c(210,27.4,13.7,13.7,13.7)
  Y <- (cbind(rep(1,N),Z))%*%beta
  true_prs <- rep(0,N)
  for(i in 1:N){
    true_prs[i] <- expit(Z[i,]%*%c(-1,0.5,-0.25,-0.1))
  }
  X <- matrix(data = NA, nrow = N, ncol = 4)
  for(i in 1:N){
    X[i,1] <- exp(Z[i,1]/2)
    X[i,2] <- Z[i,2]/(1+exp(Z[i,1])) + 10
    X[i,3] <- ((Z[i,1]*Z[i,3])/25 + 0.6)^3
    X[i,4] <- (Z[i,2] + Z[i,4] + 20)^2
  }
  rs_index <- t <- rbinom(N,1,true_prs)
  respondents_index <- which(rs_index == 1)
  Y_respondents <- Y[respondents_index]
  X_observed <- X[respondents_index,]
  
  dt_XTYZ <- data.frame(Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], Z4 = Z[,4], X1 = X[,1], X2 = X[,2], X3 = X[,3], X4 = X[,4], Y = Y, t = t)
  dt_XTYZ_res <- dt_XTYZ[respondents_index,]
  eX_specified <- fitted(glm(t~Z1+Z2+Z3+Z4,data = dt_XTYZ, family=binomial))
  eX_misspecified <- fitted(glm(t~X1+X2+X3+X4,data = dt_XTYZ, family=binomial))
  
  w1<-t/eX_specified
  W1<-w1/sum(w1)
  w0<-(1-t)/(1-eX_specified)
  W0<-w0/sum(w0)
  r_1 <- length(respondents_index)/N
  w1_mis<-t/eX_misspecified
  W1_mis<-w1_mis/sum(w1_mis)
  w0_mis<-(1-t)/(1-eX_misspecified)
  W0_mis<-w0_mis/sum(w0_mis)

  # ests[rep,1]<-sum(W1*Y)
  # ests[rep,2]<-(1-r_1)*sum(W0*Y)+r_1*sum(W1*Y)
  # 
  residulas_cor <- Y-t(coefficients(lm(Y~Z1+Z2+Z3+Z4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$Z1,dt_XTYZ$Z2,dt_XTYZ$Z3,dt_XTYZ$Z4)))
  ests_AIPW[rep,1] <- sum(coefficients(lm(Y~Z1+Z2+Z3+Z4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$Z1,dt_XTYZ$Z2,dt_XTYZ$Z3,dt_XTYZ$Z4)))/N + sum(W1*residulas_cor)
  residulas_incor <- Y-t(coefficients(lm(Y~X1+X2+X3+X4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$X1,dt_XTYZ$X2,dt_XTYZ$X3,dt_XTYZ$X4)))
  ests_AIPW[rep,2] <- sum(coefficients(lm(Y~X1+X2+X3+X4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$X1,dt_XTYZ$X2,dt_XTYZ$X3,dt_XTYZ$X4)))/N + sum(W1*residulas_incor)
  
  ests_AIPW[rep,3] <- sum(coefficients(lm(Y~Z1+Z2+Z3+Z4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$Z1,dt_XTYZ$Z2,dt_XTYZ$Z3,dt_XTYZ$Z4)))/N + sum(W1_mis*residulas_cor)
  ests_AIPW[rep,4] <- sum(coefficients(lm(Y~X1+X2+X3+X4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$X1,dt_XTYZ$X2,dt_XTYZ$X3,dt_XTYZ$X4)))/N + sum(W1_mis*residulas_incor)
  
  # w1_mis<-t/eX_misspecified
  # W1_mis<-w1_mis/sum(w1_mis)
  # w0_mis<-(1-t)/(1-eX_misspecified)
  # W0_mis<-w0_mis/sum(w0_mis)
  # 
  # ests[rep,3]<-sum(W1_mis*Y)
  # ests[rep,4]<-(1-r_1)*sum(W0_mis*Y)+r_1*sum(W1_mis*Y)
  # 
  # # Y_specified <- fitted(lm(Y~Z1+Z2+Z3+Z4,data = dt_XTYZ_res))
  # ests[rep,5] <- sum(coefficients(lm(Y~Z1+Z2+Z3+Z4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$Z1,dt_XTYZ$Z2,dt_XTYZ$Z3,dt_XTYZ$Z4)))/N
  # 
  # # Y_misspecified <- fitted(lm(Y~X1+X2+X3+X4,data = dt_XTYZ_res))
  # ests[rep,6] <- sum(coefficients(lm(Y~X1+X2+X3+X4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$X1,dt_XTYZ$X2,dt_XTYZ$X3,dt_XTYZ$X4)))/N
}

nreps_AIPW<-1000 
ests_AIPW_1000<-matrix(0,nrow=nreps,ncol=4)
for(rep in 1:nreps){
  N <- 1000
  dim_z <- 4
  R_z <- matrix(c(1,0,0,0,
                  0,1,0,0,
                  0,0,1,0,
                  0,0,0,1), 
                nrow = 4, ncol = 4, byrow = TRUE)
  Z <- mvrnorm(n = N, mu = rep(0,dim_z), Sigma = R_z)
  beta <- c(210,27.4,13.7,13.7,13.7)
  Y <- (cbind(rep(1,N),Z))%*%beta
  true_prs <- rep(0,N)
  for(i in 1:N){
    true_prs[i] <- expit(Z[i,]%*%c(-1,0.5,-0.25,-0.1))
  }
  X <- matrix(data = NA, nrow = N, ncol = 4)
  for(i in 1:N){
    X[i,1] <- exp(Z[i,1]/2)
    X[i,2] <- Z[i,2]/(1+exp(Z[i,1])) + 10
    X[i,3] <- ((Z[i,1]*Z[i,3])/25 + 0.6)^3
    X[i,4] <- (Z[i,2] + Z[i,4] + 20)^2
  }
  rs_index <- t <- rbinom(N,1,true_prs)
  respondents_index <- which(rs_index == 1)
  Y_respondents <- Y[respondents_index]
  X_observed <- X[respondents_index,]
  
  dt_XTYZ <- data.frame(Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], Z4 = Z[,4], X1 = X[,1], X2 = X[,2], X3 = X[,3], X4 = X[,4], Y = Y, t = t)
  dt_XTYZ_res <- dt_XTYZ[respondents_index,]
  eX_specified <- fitted(glm(t~Z1+Z2+Z3+Z4,data = dt_XTYZ, family=binomial))
  eX_misspecified <- fitted(glm(t~X1+X2+X3+X4,data = dt_XTYZ, family=binomial))
  
  w1<-t/eX_specified
  W1<-w1/sum(w1)
  w0<-(1-t)/(1-eX_specified)
  W0<-w0/sum(w0)
  r_1 <- length(respondents_index)/N
  w1_mis<-t/eX_misspecified
  W1_mis<-w1_mis/sum(w1_mis)
  w0_mis<-(1-t)/(1-eX_misspecified)
  W0_mis<-w0_mis/sum(w0_mis)
  
  # ests[rep,1]<-sum(W1*Y)
  # ests[rep,2]<-(1-r_1)*sum(W0*Y)+r_1*sum(W1*Y)
  # 
  residulas_cor <- Y-t(coefficients(lm(Y~Z1+Z2+Z3+Z4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$Z1,dt_XTYZ$Z2,dt_XTYZ$Z3,dt_XTYZ$Z4)))
  ests_AIPW_1000[rep,1] <- sum(coefficients(lm(Y~Z1+Z2+Z3+Z4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$Z1,dt_XTYZ$Z2,dt_XTYZ$Z3,dt_XTYZ$Z4)))/N + sum(W1*residulas_cor)
  residulas_incor <- Y-t(coefficients(lm(Y~X1+X2+X3+X4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$X1,dt_XTYZ$X2,dt_XTYZ$X3,dt_XTYZ$X4)))
  ests_AIPW_1000[rep,2] <- sum(coefficients(lm(Y~X1+X2+X3+X4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$X1,dt_XTYZ$X2,dt_XTYZ$X3,dt_XTYZ$X4)))/N + sum(W1*residulas_incor)
  
  ests_AIPW_1000[rep,3] <- sum(coefficients(lm(Y~Z1+Z2+Z3+Z4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$Z1,dt_XTYZ$Z2,dt_XTYZ$Z3,dt_XTYZ$Z4)))/N + sum(W1_mis*residulas_cor)
  ests_AIPW_1000[rep,4] <- sum(coefficients(lm(Y~X1+X2+X3+X4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$X1,dt_XTYZ$X2,dt_XTYZ$X3,dt_XTYZ$X4)))/N + sum(W1_mis*residulas_incor)
  
  # w1_mis<-t/eX_misspecified
  # W1_mis<-w1_mis/sum(w1_mis)
  # w0_mis<-(1-t)/(1-eX_misspecified)
  # W0_mis<-w0_mis/sum(w0_mis)
  # 
  # ests[rep,3]<-sum(W1_mis*Y)
  # ests[rep,4]<-(1-r_1)*sum(W0_mis*Y)+r_1*sum(W1_mis*Y)
  # 
  # # Y_specified <- fitted(lm(Y~Z1+Z2+Z3+Z4,data = dt_XTYZ_res))
  # ests[rep,5] <- sum(coefficients(lm(Y~Z1+Z2+Z3+Z4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$Z1,dt_XTYZ$Z2,dt_XTYZ$Z3,dt_XTYZ$Z4)))/N
  # 
  # # Y_misspecified <- fitted(lm(Y~X1+X2+X3+X4,data = dt_XTYZ_res))
  # ests[rep,6] <- sum(coefficients(lm(Y~X1+X2+X3+X4,data = dt_XTYZ_res))%*%t(cbind(rep(1,N),dt_XTYZ$X1,dt_XTYZ$X2,dt_XTYZ$X3,dt_XTYZ$X4)))/N
}


dt_XTYZ <- data.frame(Z1 = Z[,1], Z2 = Z[,2], Z3 = Z[,3], Z4 = Z[,4], X1 = X[,1], X2 = X[,2], X3 = X[,3], X4 = X[,4], Y = Y, t = t)

eX_specifed <- fitted(glm(t~Z1+Z2+Z3+Z4,data = dt_XTYZ, family=binomial))
eX_misspecified <- fitted(glm(t~X1+X2+X3+X4,data = dt_XTYZ, family=binomial))

w1<-t/eX_specifed
W1<-w1/sum(w1)
sum(W1[respondents_index]*Y_respondents)
