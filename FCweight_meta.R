library(metafor)
b<-0
zdata<-as.data.frame(lapply(F,as.numeric))
study <- c("metaMDD-PKU","metaMDD-ZJU","metaMDD-ChMU","metaMDD-SXYH","metaMDD-CQMU","metaMDD-SWU","metaMDD-CapMU")   # study names
ni <- c(84,39,56,58,52,321,82)       # sample sizes
for (i in 1:264){
# input of data
ri<-c(zdata[i,1],zdata[i,2],zdata[i,3],zdata[i,4],zdata[i,5],zdata[i,6],zdata[i,7])
  # raw difference of accuracy
res<-rma(measure="ZCOR",ri,ni,method="DL",slab=study) # meta-for package for meta-analysis
b[i]<-as.data.frame(lapply(res[1],as.numeric))
}
write.csv(b,file="F.csv")