install.packages("metafor",repos='https://mran.microsoft.com/snapshot/2019-02-01/')
library(metafor)

# input of data
study <- c("metaMDD-PKU","metaMDD-ZJU","metaMDD-ChMU","metaMDD-SXYH","metaMDD-CQMU","metaMDD-SWU","metaMDD-CapMU")   # study names
ri <- c(-6.32,-9.7,-10.26,-5.74,-8.33,-1.08,-5.7)  # raw difference of accuracy
ni <- c(80,26,53,58,52,218,66)       # sample sizes

res <- rma(ri,ni,method="DL",slab=study) # meta-for package for meta-analysis
summary(res)
print(res)
forest(res,xlab = "accuracy difference in females") # forest plot
funnel(res,xlab = "accuracy difference in females") # plot the heterogeneity

ranktest(res) # Egger's test for heterogeneity
