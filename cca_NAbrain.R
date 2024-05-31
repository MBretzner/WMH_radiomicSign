install.packages("CCA", method = "curl")
install.packages("gdata", method = "curl")
install.packages("bestNormalize", method = "curl")
install.packages("ggplot2", method = "curl")
install.packages("GGally", method = "curl")
install.packages("CCP", method = "curl")
install.packages("rcc", method = "curl")
install.packages("vegan", method = "curl")


library(CCA)
library(gdata)
library(bestNormalize)
library(ggplot2)
library(GGally)
library(CCP)
library(rcc)
library(vegan)

EVTdata = read.csv("/PHShome/mi362/radiomics/MRI_GENIE/NAWM/Rstudio_analysis/radiomic_signature_scaled_forR.csv")
clinical_varnames <- c("Age","Sex","AF","DM","HTN","CAD","Smoking_ever")
clinical <- EVTdata[clinical_varnames]
radiomics <- EVTdata[,14:81]
radio_varnames <- c("WMH_Bvadj","Brain_volume","Ventricle_volume")
radio <- EVTdata[radio_varnames]

#visualize data
#ggpairs(clinical)
correl <- matcor(clinical, radiomics )
img.matcor(correl, type = 2)


#CCA
cc1 <- cc(clinical,radiomics)

# display the canonical correlations
cc1$cor
par(mfrow = c(1,2))
barplot(cc1$cor, main = "Canonical correlations for 'cancor()'", col = "gray")
barplot(cc1$cor, main = "Canonical correlations for 'cancor()'", col = "gray")

# raw canonical coefficients
cc1[3:4]

cc2 <- comput(clinical,radiomics,cc1)
cc2[3:6]

# tests of canonical dimensions
rho <- cc1$cor

#get dimensions
nx<-dim(clinical)[2]      
ny<-dim(radiomics)[2]
ncv<-min(nx,ny)
cvlab<-paste("CV",1:ncv)

#Get covariance matrices
use="complete.obs"
cxx<-cov(clinical,use=use)
cyy<-cov(radiomics,use=use)
cxy<-cov(clinical,radiomics,use=use)
cyx<-t(cxy)
#Find the projections
ey<-eigen(qr.solve(cyy,cyx)%*%qr.solve(cxx,cxy))
ex<-list(values=ey$values,vectors=qr.solve(cxx,cxy)%*%(ey$vec))
cc1$corr<-(ex$val[1:ncv])^0.5
names(cc1$corr)<-cvlab

# variance explained by structural correlations
cc1$xstructcorrsq<-cc1$xcoef^2 
cc1$ystructcorrsq<-cc1$ycoef^2 


# #Find the canonical communalities (total var explained)
cc1$xcancom<-apply(cc1$xstructcorrsq,1,sum) 
cc1$ycancom<-apply(cc1$ystructcorrsq,1,sum) 

#Find the canonical variate adequacies (Fraction of Total Variance Explained by Each CV, Within Sets):
cc1$xcanvad<-apply(cc1$xstructcorrsq,2,mean)      
cc1$ycanvad<-apply(cc1$ystructcorrsq,2,mean)

#Find the redundancy indices (Rd) for X|Y and Y|X
cc1$corrsq <- cc1$corr^2
cc1$xvrd<-cc1$xcanvad*cc1$corrsq  
cc1$yvrd<-cc1$ycanvad*cc1$corrsq
cc1$xrd<-sum(cc1$xvrd)
cc1$yrd<-sum(cc1$yvrd)

#plot var
plt.cc(cc1,d1=1,d2=3, var.label = FALSE, ind.names = NULL,type="v",Xnames=NULL,Ynames = NULL)

#indiv names add ind.names=EVTdata[,1] in

# Define number of observations, number of variables in first set, and number of variables in the second set.
n <- dim(clinical)[1]
p <- length(clinical)
q <- length(radiomics)

## Calculate p-values using the F-approximations of different test statistics:
cc1$wilks <- p.asym(rho, n, p, q, tstat = "Wilks")
p.asym(rho, n, p, q, tstat = "Pillai")

#stabdardized coefficients
std_coef1<-diag(sqrt(diag(cov(clinical))))
std_coef1%*%cc1$xcoef

std_coef2<-diag(sqrt(diag(cov(radiomics))))
std_coef2%*%cc1$ycoef

mod <- cca(clinical,radiomics)
plot(cc1$cor,type="b")
plt.cc(cc1)


#try vegan package
vegan_CCA <- CCorA(clinical,radiomics,stand.Y = TRUE,stand.X = TRUE,permutations = 999)
biplot(vegan_CCA,"biplots",xlabs=NA,int=0.5,)

#normalize variables
radiomics_BN <- lapply(radiomics,bestNormalize)
radiomics_xt <- lapply(radiomics_BN, '[[', 1)
radiomics_xt_df <- data.frame(radiomics_xt)
clinical_test <- clinical
Age_BN <- bestNormalize(clinical$Age)
clinical_test$Age <- Age_BN$x.t
clinical_BN <- lapply(clinical,bestNormalize)
clinical_xt <- lapply(clinical_BN, '[[', 1)
clinical_xt_df <- data.frame(clinical_xt)

vegan_CCA_xt <- CCorA(clinical_test,radiomics_xt_df,stand.Y = TRUE,stand.X = TRUE,permutations = 999)
vegan_CCA_xt_noperm <- CCorA(clinical_test,radiomics_xt_df,stand.Y = TRUE,stand.X = TRUE,permutations = 0)

cc1_xt <- cc(clinical_test,radiomics_xt_df)

rho <- cc1_xt$cor
n <- dim(clinical_xt_df)[1]
p <- length(clinical_xt_df)
q <- length(radiomics_xt_df)
cc1_xt$wilks <- p.asym(rho, n, p, q, tstat = "Wilks")
cc1_xt$wilks$p.value.FDRadj <- p.adjust(cc1_xt$wilks$p.value ,method='fdr')
cc1_xt$wilks$p.value.Bonferroni <- p.adjust(cc1_xt$wilks$p.value ,method='bonferroni')
biplot(vegan_CCA_xt,"biplots",xlabs=NA,int=0.5)

#Plot the data on each canonical variate
`plot.canoncorr` <-
  function(x,...){
    rbPal <- colorRampPalette(c('blue','red'))
    x$Color <- rbPal(10)[as.numeric(cut(x$Age,breaks = 10))]
    ncv = length(x)
    cc1_xt$Age <- clinical$Age
    for(i in 1:ncv){
      par(pty="s")
      plot(x$scores$xscores[,i],x$scores$yscores[,i],xlab="Clinical",ylab="Radiomics",main=paste("Canonical Variate Plot - Variate",i,sep=" "),
           cex=0.5,pch=20,col=x$Color,
           #col=rgb(red=30/255, green=144/255, blue=255/255, alpha=0.4),
           #sub = paste("r=",round(x$cor[i],digits=2),sep="")
      )
      legend("bottomright", bty="n", legend=paste("r=", format(cc1_xt$cor[i], digits=2)))
      abline(mean(x$scores$yscores[,i],na.rm=TRUE)-x$cor[i]*mean(x$scores$xscores[,i],na.rm=TRUE),x$cor[i])
      gradientLegend(valRange = c(18,94),pos = 0.5,color =alphaPalette(c('blue','red'), f.seq=seq(0,1, by=.1)),inside=TRUE)
      #text(mean(x$scores$xscores[,i],na.rm=TRUE),mean(x$scores$yscores[,i],na.rm=TRUE),label=paste("r=",round(x$cor[i],digits=2),sep=""),pos=3,srt=180/pi*atan(x$cor[i]))    
    }}

plot.canoncorr(cc1_xt)


plt.cc(cc1_xt,d1=1,d2=2, var.label = TRUE, ind.names = NULL,type="v",Xnames=NULL,Ynames = NULL)
cc2_xt <- cc(radiomics_xt_df,clinical_test)
plt.cc(cc2_xt,d1=1,d2=2, var.label = FALSE, ind.names = NULL,type="v",Xnames=NULL,Ynames = NULL)



library('formattable')
vegan_CCA_xt$Eigenvalues_perc = percent(vegan_CCA_xt$Eigenvalues)
barplot(vegan_CCA_xt$Eigenvalues, main="Scree Plot",xlab="Canonical Dimension",
        names.arg = cvlab,cex.names = 0.7)

vegan_CCA_xt_permmax <- CCorA(clinical_test,radiomics_xt_df,stand.Y = TRUE,stand.X = TRUE,permutations = 999999)
