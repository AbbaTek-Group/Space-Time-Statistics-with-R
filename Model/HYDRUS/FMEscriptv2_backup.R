unloadNamespace("dream")
options(repos='http://cran.rstudio.com/')
library(dream)
library(hydrusR)
library(snow)
library(FME)
library(latticeExtra)


#HYDRUSpars1<-list(thr1=0.1096,ths1=0.56,Alfa1=0.024,n1=1.3,Ks1=1.2,l1=3)
#HYDRUSpars2<-list(GWLOL=2100, Aqh=-0.667, Bqh=-0.0083)
#HYDRUSpars3<-list(thr1=0.1096,ths1=0.5 ,Alfa1=0.0298,n1=1.2732,Ks1=0.0000001,l1=-3,
#                  thr2=0.1596,ths2=0.65,Alfa2=0.0598,n2=1.473 ,Ks2=0.0001   ,l2=-1,
#                  thr3=0.396 ,ths3=0.75,Alfa3=0.0698,n3=9     ,Ks3=0.01     ,l3=3)

# initialize
count2<-0
run<-1

#set.seed(123)
HYDRUSpar.range<-list(GWLOL=c(2100,3500),Aqh=c(-0.667,-0.001),Bqh=c(-0.01251,-0.0083),
                      thr1=c(0.1096,0.7),ths1=c(0.5,0.7),Alfa1=c(0.0298,20),n1=c(1.2732,9),Ks1=c(0.0000001,0.5),l1=c(-3,3),
                      thr2=c(0.1096,0.7),ths2=c(0.5,0.7),Alfa2=c(0.0298,20),n2=c(1.2732,9),Ks2=c(0.0000001,0.5),l2=c(-3,3),
                      thr3=c(0.1096,0.7),ths3=c(0.5,0.7),Alfa3=c(0.0298,20),n3=c(1.2732,9),Ks3=c(0.0000001,0.5),l3=c(-3,3))

Model.h<-function(HYDRUSpars){
  
  # read in observational data
  Qobs<-read.table("thetaObs2.txt",sep ="|", header = TRUE)[c(1:6000),1]

# script that updates Selector.in values with newly selected variables from mcmc
  system("awk 'NR<23' Selector.orig.in | awk '{print $0}' > Selector.in")
  HYDRUSparsform <- matrix(HYDRUSpars[1:3], nrow = 1, byrow=TRUE)
  write.table(HYDRUSparsform, "Selector.in", quote=FALSE,append =TRUE,row.names = FALSE,col.names = FALSE)
  system("awk 'NR>23 && NR<29' Selector.orig.in | awk '{print $0}' >> Selector.in")
  HYDRUSparsform <- matrix(HYDRUSpars[4:21], nrow = 3, byrow = TRUE)
  write.table(HYDRUSparsform,"Selector.in",quote=FALSE,append=TRUE,row.names=FALSE,col.names=FALSE)
  system("awk 'NR>31' Selector.orig.in | awk '{print $0}' >> Selector.in")

  # runs hydrus1d
  system("./H1D_CALC.EXE > H1DIO.out ")
  # cleans up hydrus1d output to match Qob.
  system("./h1dout2obs.sh")
 
  # is.nan statement and replace with blank to prevent crash
  nancheck.data <- read.table("OBS_NODE_mod.OUT", header=T, sep="", quote ="", 
                         dec = ".",nrows = 6000,skip=10,fill=T,stringsAsFactors = F)
  nancheck.theta<-nancheck.data[3]
  check<-as.numeric(unlist(nancheck.theta))
  if (is.nan(sum(check))){
    system("cp OBS_NODE_blank.OUT OBS_NODE_mod.OUT")
    print("is.nan==TRUE")
    count2<-count2+1
  }

  # this catches abrupt code exit and output shortfall
  system("./wccheck.sh")

  Qsim.data <- read.table("OBS_NODE_mod.OUT", header=T, sep="", quote ="", dec = ".",
                          nrows = 6000,skip=10,fill=T,stringsAsFactors = F)
  Qsim.all.theta<-Qsim.data[3]
  Qs<-as.numeric(unlist(Qsim.all.theta))

   
  # Calculate residual and preproc
  res<-Qobs-Qs
  mean.res<-mean(res)
  mean.abs.res<- mean(abs(res))
 
#  print(res)
#  print((res_k))
#  print((function(res_k)))

  logp<-sum(sapply(res,function(res_k)log((1/2*mean.abs.res))*exp(-abs(res_k)/(mean.abs.res))))

  logp<-abs(logp)
   
  # Generate output
  runid <- data.frame(run,matrix(HYDRUSpars,nrow=1,byrow=TRUE))
  results <- data.frame(run,mean.res,logp)
  write.table(runid[1,], file = "dreamparm.dat", append = TRUE, col.names=F, row.names=F)
  write.table(results[1,], file = "dreamout.dat", append = TRUE, col.names=F, row.names=F)

  run<-run+1
  print(paste0("NaN count :",count2))
  print(paste0("Success count :",run))

  return(logp) 

}
  
run<-1                      
count<-0                      
# Headers for output files
write("runID GWLOL Aqh Bqh thr1 ths1 Alfa1 n1 Ks1 l1 thr2 ths2 Alfa2 n2 Ks2 l2 thr3 ths3 Alfa3 n3 Ks3 l3",file = "dreamparm.dat", append = FALSE)
write("runID mean(res) logp",file = "dreamout.dat", append = FALSE)
 
#control<-list(nseq=14, ndraw=1000)
control<-list(nseq=4, ndraw=8)

dd<-dream(FUN=Model.h,pars=HYDRUSpar.range,func.type='logposterior.density',control=control)
  

print(paste0("NaN count :",count2))
print(paste0("Success count :",run))

print(summary(dd))
print(coef(dd))

pdf('ddplots.pdf')
plot(dd)
dev.off()
