
# Refer to parameter-estimation file for the variables you need to input if for a generic study
# You need to run the parameter-estimate file if you want to find MLE parameters.

## Uncomment the following if you want to do SMAPVEX12 analysis
#source("input_file_paper1.R")



#######################################Prediction#############################################
# The code is for M=2, extending to other M values if trivial

sp1 <- 400 ## Assumed grid spacing for the prediction grid. The finer the grid  more will be the computation time

e_rs1 <- extent(min(f1d$X),max(f1d$X),min(f1d$Y),max(f1d$Y)) ## Extent for the fine-scale points. 

xgrd_s_p <- seq(from= e_rs1[1] , by=(sp1), to = e_rs1[2])
ygrd_s_p <- seq(from= e_rs1[4] , by=(-sp1), to = e_rs1[3])
grd_as_p <-  as.matrix(expand.grid(xgrd_s_p,ygrd_s_p)) ## defining a grid for the extent of the fine-scale points.
grd_as1_p <- extract(raster_clay_airborne,grd_as_p)   ## Since, we only want those points that are in the airborne
ind_grd_as_p <- which(is.na(grd_as1_p))               ## domain, we delete the points which fall outside the airborne
                                                      ## pixels
if(length(ind_grd_as_p)==0)            
{
  grd_as_del_p <- grd_as_p
}else
{
  
  grd_as_del_p <- grd_as_p[-ind_grd_as_p,]
  
}


pred.grd <- grd_as_del_p  ## points inside the airborne pixels on which we need to make predictions

LAI_pred <- extract(raster_LAI1,pred.grd)  ## Extracting covariates for the prediction grid. you need to have rasters/shapefiles for covariates for this to work
Clay_pred <- extract(raster_clay_airborne,pred.grd) ## you can change these covariates for your own study
Silt_pred <- extract(raster_silt_airborne,pred.grd)

##Standaradizing
LAI_pred <- (LAI_pred-mean(LAI))/sd(LAI)  ## We need to standardize the prediction covariates by the same mean and 
Clay_pred <- (Clay_pred-mean(Clay))/sd(Clay) ## sd that we took for the observations (refer input_parameter file)
Silt_pred <- (Silt_pred-mean(Silt))/sd(Silt)


Covariates_pred <- cbind(LAI_pred,Clay_pred,Silt_pred)  ## removing points for which there is no value of any covariate
ind_cov_pred <- which(is.na(Covariates_pred),TRUE)[,1]  ## For this study, those points were mostly zero
ind_cov1_pred <- unique(ind_cov_pred)
if(length(ind_cov1_pred)==0)
{
  Covariates_pred=Covariates_pred
}else
{
  
  Covariates_pred <- Covariates_pred[-ind_cov1_pred,]
  pred.grd <- pred.grd[-ind_cov1_pred,]
  
}


X_wt_pred <- cbind(rep(1,dim(Covariates_pred)[1]),Covariates_pred) # Defining the X_matrix for prediction for mean and
                                                                   #covariance as before
X_1_pred <- X_wt_pred
pred.grd <- pred.grd/100000

dist.mat_pred <- rdist(pred.grd,obs.grd) #distance matrix between observations and predictions
pred_mat <- rdist(pred.grd,pred.grd)  # distance matrix between predictions


theta <- rep(1, 10) # input the MLE parameters for the covariance, estimated from the parameter estimation file
                           # here, we are just inputting random values

Beta <- rep(1, 4) ## input value of beta corresponding to the p_covariance as calculated in the parameter-estimation file
                  ## here, we are just inputting random values
p <- as.numeric(c(theta, Beta))


library(gdata) 
keep(p, X_wt_pred, X_1_pred, X_wt, X_1, dist.mat, dist.mat_pred, pred_mat, z, z11, input, e_rs1, pred.grd, sure=TRUE)
detach("package:gdata", unload=TRUE) 

###################################################################################################################
#####################################Prediction code###################################################################
###################################################################################################################


w_denom <- 1+exp((X_wt)%*%c((p[1]),(p[2]),(p[3]),(p[4]))) ## Refer to parameter-estimation file for documentation
w1 <- 1/w_denom
rm(w_denom)
w2 <- (1-w1)


ff_w1 <- function(a)
{
  a*w1
}

ff_w2 <- function(a)
{
  a*w2
}


ff_w1_pred <- function(a)
{
  a*w1_pred
}

ff_w2_pred <- function(a)
{
  a*w2_pred
}

w_denom_pred <- 1+exp((X_wt_pred)%*%c((p[1]),(p[2]),(p[3]),(p[4])))
w1_pred <- 1/w_denom_pred
rm(w_denom_pred)
w2_pred <- (1-w1_pred)

Mat1 <-  matrix(NA, nrow=dim(X_wt)[1],ncol=dim(X_wt)[1]) #C(S,S)
Mat1<-  p[5]*apply(w1,1,ff_w1)*Matern(dist.mat,range=p[6],nu=p[7])+p[8]*apply(w2,1,ff_w2)*Matern(dist.mat,range=p[9],nu=p[10])
rm(dist.mat)
Ch <- chol(Mat1)
A1 <- chol2inv(Ch) 

Mean_obs <- (X_1) %*% c(p[11],p[12],p[13],p[14]) #Mean for the observations

Mat1_pred <-  matrix(NA, nrow=dim(X_wt_pred)[1],ncol=dim(X_wt_pred)[1]) #C(Sp,S)
Mat1_pred<-  p[5]*t(apply(w1_pred,1,ff_w1))*Matern(dist.mat_pred,range=p[6],nu=p[7])+p[8]*t(apply(w2_pred,1,ff_w2))*Matern(dist.mat_pred,range=p[9],nu=p[10])
rm(dist.mat_pred)

Mat_for_pred<-  matrix(NA, nrow=dim(X_wt_pred)[1],ncol=dim(X_wt_pred)[1])
Mat_for_pred<-  p[5]*t(apply(w1_pred,1,ff_w1_pred))*Matern(pred_mat,range=p[6],nu=p[7])+p[8]*t(apply(w2_pred,1,ff_w2_pred))*Matern(pred_mat,range=p[9],nu=p[10])

rm(pred_mat)

Mean_pred <- (X_1_pred) %*% c(p[11],p[12],p[13],p[14]) # Mean for the predictions

Mean_pred_total <- Mean_pred + Mat1_pred %*%  A1 %*% (z-Mean_obs) #Mean of posterior predictions at fine scale
Mean_pred_total <- Mean_pred_total *sd(z11)+mean(z11) ## Converting mean predictions to original scale (remember we
                                                      ## saved z to variable z11 before scaling)
Cov_pred <- Mat_for_pred - (Mat1_pred %*%  A1 %*% t(Mat1_pred)) # Covariance of posterior predictions at fine scale
Cov_pred <- Cov_pred*(sd(z11)^2) ## Converting Covariance to original scale

################################################################################################################
########################################################Airborne raster############################################
#################################################################################################################


## This is SMAPVEX12 specific analysis. You can input your own raster file for a generic study
sr <- "+proj=utm +zone=14N +ellps=GRS80 +datum=NAD83"
a1d<- input$Air    ##airborne data (refer input file)
colnames(a1d) <- c("Easting","Northing","SM")
grd_a <- (a1d[,1:2])
raster_airborne <- rasterFromXYZ(cbind(grd_a,a1d$SM),crs=sr) ## converting points to raster

r2 <- crop(raster_airborne,e_rs1) # Cropping the airborne raster to the study extent

## Mean upscaling
# Mean of posterior prediction at airborne scale
raster_pred_mean <- rasterize(pred.grd*100000,r2,Mean_pred_total,fun=mean,na.rm=FALSE) # Mean prediction is straightforward
                                                                                       # as it is simply the average
raster_pred_mean <- mask(raster_pred_mean,r2) # taking the same extent as domain
# Covariance upscaling
ss_a <- cellFromXY(r2,(pred.grd*100000)) #get cell numbers from X, Y coordinates of raster r2
ss_a[is.na(ss_a)] <- 0 ## Change NA to 0

Cov_upscale <- matrix(0, nrow=r2@ncols*r2@nrows,ncol=length(ss_a))

for(i in 1: (r2@ncols*r2@nrows))
{
  zz <- ss_a
  nn <- sum(zz==i)  # This counts the number of points in the i^(th) pixel
  zz[which(zz!=i)]=0  ## put every point = 0 which is not in the i^(th) pixel
  zz[which(zz!=0)]=1/nn # takes the number of points in a raster as explained in the manuscript
  Cov_upscale[i,] <- zz
}

pred_covar_air <- Cov_upscale %*% Cov_pred %*% t(Cov_upscale) ## By a simple matrix manipulation, double summation
                                                              ## can be replaced by a matrix multiplication
                                                              ## this speeds up computation in R as it is slow in for loops
pred_var_a3 <- diag(pred_covar_air)


coords <- xyFromCell(r2,(1:length(r2))) ##coordinates of the center of raster r2
r6 <- rasterize(coords,r2,pred_var_a3)  ##rasterizing pred var

# variance of posterior predictions at airborne scale
raster_pred_var <- mask(r6,r2) # taking the same extent as the domain    



