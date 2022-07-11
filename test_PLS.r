install.packages('R.matlab')
library('R.matlab')
install.packages('pls')
library('pls')
rawdata <- readMat(file.path("F:/work_dir/AHBAenrich/test11.mat"));
typeof(rawdata$X)
pls1<-plsr(rawdata$Y~.,data=as.data.frame(rawdata$X),ncomp=15,validation="LOO",jackknife=TRUE,method="widekernelpls")
summary(pls1,what="all")
