######## Permutation test ##########

# Define parameters

## load file 
file=read.csv("comparison of score PHE.csv",h=T,sep=";",dec=".") 
attach(file)
## title of your analysis
analysis_title="test"
## define the factor tested
factor_tested=Variation.PEG
## define the parameter tested 
pheno=PHE.Genes
## define the levels of the two group 
gp_1="Overused"
gp_2="equallyused"
## define the number of replications
replications=1000
# test (don't change nothing!)

pdf(paste(analysis_title,"_permutation.pdf"))
mlm=cbind(factor_tested,pheno)
mlm=as.data.frame(mlm)
summary=NULL
summary2=NULL
pheno=subset(mlm$pheno,mlm$pheno!="NA")
var2=1:length(pheno)
TABLE=NULL
TABLE2=NULL
s=subset(mlm,as.factor(mlm$factor_tested)==gp_1)
s=subset(s,s$pheno!="NA")
s_mean=mean(as.numeric(s[,2]))
s_median=median(as.numeric(s[,2]))
delta=subset(mlm,as.factor(mlm$factor_tested)==gp_2)
delta=subset(delta,delta$pheno!="NA")
size=length(delta[,1])
delta_mean=mean(as.numeric(delta[,2]))
delta_median=median(as.numeric(delta[,2]))
delta_mean=delta_mean-s_mean
delta_median=delta_median-s_median
for(j in 1:replications)
  {
    sample=sample(1:length(pheno),size,replace = F)  
    sample2=setdiff(var2,sample)  
    sample2=pheno[sample2]
    sample=pheno[sample]
    sample_mean=mean(as.numeric(sample))
    sample2_mean=mean(as.numeric(sample2))
    sample_mean=sample_mean-sample2_mean
    TABLE=c(TABLE,sample_mean)
    sample_med=median(as.numeric(sample))
    sample2_med=median(as.numeric(sample2))
    sample_med=sample_med-sample2_med
    TABLE2=c(TABLE2,sample_med)
    }
med=median(TABLE)
t=sort(TABLE)
t1=t[0.025*replications]
t2=t[0.975*replications]
x=subset(t,t<=delta_mean)
signi1=length(x)/replications
x=subset(t,t>=delta_mean)
signi2=length(x)/replications
signi=min(c(signi1,signi2))
if (signi==0)
  {
    signi=1/replications
    signi=paste("<",signi)}
x=cbind(gp_1,gp_2,mean(as.numeric(s[,2])),mean(as.numeric(delta[,2])),delta_mean,t1,med,t2,signi)
summary=rbind(summary,x)
hist(TABLE,xlab="diff expected (mean)",ylab="frequency", main=paste(gp_2,"-",gp_1),col="grey",nclass=20,cex.main=1, cex.lab=1)
abline(v=med, col="red", lwd=3)
abline(v=t1, col="black", lwd=2)
abline(v=t2, col="black", lwd=2)
abline(v=delta_mean, col="black", lwd=1, lty=2)

med=median(TABLE2)
t=sort(TABLE2)
t1=t[0.025*replications]
t2=t[0.975*replications]
x=subset(t,t<=delta_median)
signi1=length(x)/replications
x=subset(t,t>=delta_mean)
signi2=length(x)/replications
signi=min(c(signi1,signi2))
if (signi==0)
{
  signi=1/replications
  signi=paste("<",signi)}
x=cbind(gp_1,gp_2,median(as.numeric(s[,2])),median(as.numeric(delta[,2])),delta_median,t1,med,t2,signi)
summary2=rbind(summary2,x)
hist(TABLE2,xlab="diff expected (median)",ylab="frequency", main=paste(gp_2,"-",gp_1),col="grey",nclass=20,cex.main=1, cex.lab=1)
abline(v=med, col="red", lwd=3)
abline(v=t1, col="black", lwd=2)
abline(v=t2, col="black", lwd=2)
abline(v=delta_median, col="black", lwd=1, lty=2)

summary=as.data.frame(summary)
colnames(summary)=c("first_group","second_group","Mean_first_group","Mean_second_group","Mean_delta_group","2.5_Mean_permutations","Median_Mean_permutations","97.5_Mean_permutations","Pvalue")
write.table(summary, paste(analysis_title,"_mean_permutation_test.csv"), sep=";", quote= FALSE)
summary2=as.data.frame(summary2)
colnames(summary2)=c("first_group","second_group","Median_first_group","Median_second_group","Median_delta_group","2.5_Median_permutations","Median_Median_permutations","97.5_Median_permutations","Pvalue")
write.table(summary2, paste(analysis_title,"_median_permutation_test.csv"), sep=";", quote= FALSE)

dev.off()
