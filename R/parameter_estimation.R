library(raster)
library(fields)
library(GenSA)
library(rgenoud)


## Note: You can use the following 3 codes (M=1, 2 and 3) for your analysis by putting in values of z for SM, obs.grd for your coordinates
## and by defining your own X_wt (covariate matrix for covariance) and X_1 (covariate matrix for mean). For n number of SM observations 
## and p covariates for mean, the dimension of X_1 matrix will be n*p . Similarly for X_wt.
##Extension to other values of M is trivial.

## Uncomment the following if you want to do SMAPVEX12 analysis
#source("input_file_paper1.R")

z11 <- z ## Save the value of z (SM value) to z11 before scaling as it will help while doing predictions
z <- scale(z) ## Standardizing our SM data

obs.grd <- obs.grd/100000 ## We divide the x,y coordinate for easier estimation of the range parameter
dist.mat <- rdist(obs.grd)  ## distance matrix used for the analysis

###################################################################################################################
#####################################Parameter estimation###################################################################
###################################################################################################################

############################################ M = 1 #############################################################
par_start1 <- c(1,1,1)  ## We chose random values as starting parameters. Good values of starting parameters will
                        ## help in faster convergence
ell_exp1 <- function(p){
 print(p)
 Mat1<-  p[1]* Matern(dist.mat,range=p[2],nu=p[3])
 Sigma <- Mat1 #+diag(p[4], dim(Mat1)[1]) Uncomment if you want to include nugget
 rm(Mat1)
 Ch <- chol(Sigma)  ## We use the cholesky decomposition to calculate matrix inverse as it is quicker and more stable
 A1 <- chol2inv(Ch)  ## Matrix inverse C^(-1)
 A2 <- t(X_1) %*% A1 ## X'C^(-1)
 A3 <- solve(A2 %*% X_1) ## (X'C^(-1)X)^-1
 P <- A1 - (t(A2) %*% A3 %*% A2) ## C^(-1) - C^(-1)X %*% (X'C^(-1)X)^-1 %*% X'C^(-1)
 quad.form <- t(z) %*% P %*% z  ## quadratic part of the profile log-likelihood  z'Pz
 det.part <- 2*sum(log(diag(t(Ch))))  ## Determinant part of the profile log-likelihood log|C|
 L <-  det.part + quad.form ## twice the negative of profile log-likelihood
 beta <-  A3 %*% A2 %*% z ## Note that you dont need to find beta for every iteration but only for the final parameter values
                          ## We include it here for the sake of completeness
 print(L)
 return(L)
}

lower=c(0.00001,0.00001,0.00001) ## lower bounds for parameters
upper=c(1.5,10, 1.5) # upper bounds for the parameters
out1 <- GenSA(fn=ell_exp1,par=par_start1,lower=lower,upper=upper) ## For the isotropic case, both GenSA and optim
par_list1 <- c(out1$par,out1$value)                              ## using L-BFGS-B give the same results and optim
                                                                 ## is much faster
out1 

#out1 <- genoud(fn=ell_exp1,nvars=length(par_start1),boundary.enforcement=2,starting.values=par_start1,Domains=cbind(lower,upper))
# GenSA performs consistently better than genoud


############################################# M = 2###############################################################
par_start1 <- rep(1, 10) ## Good starting values will certainly help in faster convergence times, this was not explored
ell_exp <- function(p){
  print(p)
  w_denom <- 1+exp((X_wt)%*%c((p[1]),(p[2]),(p[3]),(p[4])))  ## denominator of the weighting function
  w1 <- 1/w_denom  ## first weighting function
  rm(w_denom)
  w2 <- (1-w1)  ## second weighting function
  
  
  ff_w1 <- function(a)
  {
    a*w1
  }
  
  ff_w2 <- function(a)
  {
    a*w2
  }
  
  Mat1 <-  matrix(NA, nrow=dim(X_wt)[1],ncol=dim(X_wt)[1])
  Mat1<-  p[5]*apply(w1,1,ff_w1)*Matern(dist.mat,range=p[6],nu=p[7])+ ## We use the Matern for individual C_j. 
          p[8]*apply(w2,1,ff_w2)*Matern(dist.mat,range=p[9],nu=p[10])
  Sigma <- Mat1 #+diag(p[11], dim(Mat1)[1]) Uncomment if you want to include nugget
  rm(Mat1)
  Ch <- chol(Sigma) ## The following steps are the same as M =1. look for documentation above 
  A1 <- chol2inv(Ch)  
  A2 <- t(X_1) %*% A1  
  A3 <- solve(A2 %*% X_1)  
  P <- A1 - (t(A2) %*% A3 %*% A2) 
  quad.form <- t(z) %*% P %*% z  
  det.part <- 2*sum(log(diag(t(Ch))))
  L <-  det.part + quad.form 
  beta <-  A3 %*% A2 %*% z 
  print(L)
  return(L)
}

lower=c(-20,-20,-20,-20,0.00001,0.00001,0.00001,0.00001 ,0.00001,0.00001) ## lower bounds for parameters
upper=c(20,20,20,20,10,1.5,1.5,10,1.5,1.5) # upper bounds for the parameters
out2 <- GenSA(fn=ell_exp,par=par_start1,lower=lower,upper=upper)  ## Note that GenSA keeps on iterating long after
par_list2 <- c(out2$par,out2$value)                               ## it finds the optimum value, so its best to keep
                                                                  ## on checking the output file or put constraints
out2                                                              ## on iteration times.

#out2 <- genoud(fn=ell_exp,nvars=length(par_start1),boundary.enforcement=2,starting.values=par_start1,Domains=cbind(lower,upper))
# GenSA performs consistently better than genoud


#################################################### M = 3 #######################################################

par_start1 <- rep(1, 17)  ## The following code is similar to M = 2, with an extra weighting function
                                        ## and an extra Matern function
ell_exp3 <- function(p){
  print(p)
  w_denom <- 1+exp((X_wt)%*%c((p[1]),(p[2]),(p[3]),(p[4]))) +
    exp((X_wt)%*%c((p[11]),(p[12]),(p[13]),(p[14])))
  w1_num <- 1
  w2_num <- exp((X_wt)%*%c((p[1]),(p[2]),(p[3]),(p[4])))
  w3_num <- exp((X_wt)%*%c((p[11]),(p[12]),(p[13]),(p[14])))
  
  w1 <- w1_num/w_denom; w2 <- w2_num/w_denom; w3 <- w3_num/w_denom
  
  rm(w_denom)
  
  
  ff_w1 <- function(a)
  {
    a*w1
  }
  
  ff_w2 <- function(a)
  {
    a*w2
  }
  
  ff_w3 <- function(a)
  {
    a*w3
  }
  
  Mat1<-  p[5]*apply(w1,1,ff_w1)*Matern(dist.mat,range=p[6],nu=p[7])+
    p[8]*apply(w2,1,ff_w2)*Matern(dist.mat,range=p[9],nu=p[10])+
    p[15]*apply(w3,1,ff_w3)*Matern(dist.mat,range=p[16],nu=p[17])
  
  Sigma <- Mat1 #+diag(p[18], dim(Mat1)[1])
  rm(Mat1)
  Ch <- chol(Sigma)  ## The following steps are the same as M =1. look for documentation above 
  A1 <- chol2inv(Ch)  
  A2 <- t(X_1) %*% A1
  A3 <- solve(A2 %*% X_1)
  P <- A1 - (t(A2) %*% A3 %*% A2)
  quad.form <- t(z) %*% P %*% z
  det.part <- 2*sum(log(diag(t(Ch))))
  L <-  det.part + quad.form 
  beta <-  A3 %*% A2 %*% z 
  print(L)
  return(L)
}

lower=c(-20,-20,-20,-20,0.00001,0.00001,0.00001,0.00001 ,0.00001,0.00001, -20,-20,-20,-20, 0.00001 ,0.00001,0.00001)
upper=c(20,20,20,20,10,1.5,1.5,10,1.5,1.5, 20, 20, 20, 20, 10, 1.5, 1.5)
out3 <- GenSA(fn=ell_exp3,par=par_start1,lower=lower,upper=upper)
par_list3 <- c(out3$par,out3$value)
out3 

