reverseDirection <-
function(nSim=1000,n=100,MAF=c(0.5),
gamma0=0,gammaG=c(0.2),varX=1,
measurementError=F,delta0=0,deltaX=1,varME=1,
beta0=0,betaX=seq(from=0,to=1,length.out=4),
pleiotropy=F,betaG=c(0.2),varY=0.2,
unmeasuredConfounding=F,meanU=0,varU=1,gammaU=1,betaU=1,
sig.level=0.05,SEED=1,plot.pdf=T,plot.name="reverseDirection"){
 
################################################################################   
# load libraries
################################################################################
   # load library for Steiger test
   library(psych)

################################################################################    
# Set the seed
################################################################################
    set.seed(SEED)
    
################################################################################
# Simple Error checks
################################################################################
if(n<0 | n==0 | floor(n)!=ceiling(n) ){stop("n must be an integer greater than 0")}
if(nSim<0 | nSim==0 | floor(nSim)!=ceiling(nSim) ){stop("nSim must be an integer greater than 0")}
if(SEED<0 | SEED==0 | floor(SEED)!=ceiling(SEED) ){stop("SEED must be an integer greater than 0")}
if(length(sig.level)!=1| sig.level<0 | sig.level>1){stop("sig.level must be a single value greater than 0 and less than 1")}
if(length(unique(betaX))<2){stop("betaX must be a vector with at least two values")}

if(!varX>0){stop("varX must be greater than 0")}
if(length(varX)!=1){stop("length(varX) must equal 1")}
if(!varY>0){stop("varY must be greater than 0")}
if(length(varY)!=1){stop("length(varY) must equal 1")}
if(!varME>0){stop("varME must be greater than 0")} 
if(length(varME)!=1){stop("length(varME) must equal 1")}

if(length(gamma0)!=1){stop("length(gamma0), the intercept, must equal 1")}
if(length(delta0)!=1){stop("length(delta0), the intercept, must equal 1")}
if(length(beta0)!=1){stop("length(beta0), the intercept, must equal 1")}

if(length(MAF)!=length(gammaG)){stop("length(MAF) must equal length(gammaG) for the same number of SNPs")}
if(length(betaG)!=length(MAF)){stop("length(MAF) must equal length(betaG) for the same number of SNPs")}
if(length(deltaX)!=1){stop("length(deltaX) must equal 1")}


if(plot.pdf==T|plot.pdf=="T"|plot.pdf=="t"|plot.pdf=="True"|plot.pdf=="true"){plot.pdf<-T}  
if(plot.pdf==F|plot.pdf=="F"|plot.pdf=="f"|plot.pdf=="False"|plot.pdf=="false"){plot.pdf<-F} 
if(plot.pdf!=T & plot.pdf!=F){stop("plot.pdf must equal True or False. Please note that R is case sensitive.")}

if(pleiotropy==T|pleiotropy=="T"|pleiotropy=="t"|pleiotropy=="True"|pleiotropy=="true"){pleiotropy<-T}  
if(pleiotropy==F|pleiotropy=="F"|pleiotropy=="f"|pleiotropy=="False"|pleiotropy=="false"){pleiotropy<-F} 
if(pleiotropy!=T & pleiotropy!=F){stop("pleiotropy must equal True or False. Please note that R is case sensitive.")}

if(measurementError==T|measurementError=="T"|measurementError=="t"|measurementError=="True"|measurementError=="true"){measurementError<-T}  
if(measurementError==F|measurementError=="F"|measurementError=="f"|measurementError=="False"|measurementError=="false"){measurementError<-F} 
if(measurementError!=T & measurementError!=F){stop("measurementError must equal True or False. Please note that R is case sensitive.")}

if(unmeasuredConfounding==T|unmeasuredConfounding=="T"|unmeasuredConfounding=="t"|unmeasuredConfounding=="True"|unmeasuredConfounding=="true"){unmeasuredConfounding<-T}  
if(unmeasuredConfounding==F|unmeasuredConfounding=="F"|unmeasuredConfounding=="f"|unmeasuredConfounding=="False"|unmeasuredConfounding=="false"){unmeasuredConfounding<-F} 
if(unmeasuredConfounding!=T & unmeasuredConfounding!=F){stop("unmeasuredConfounding must equal True or False. Please note that R is case sensitive.")}
################################################################################
# Matrix to save Results
################################################################################
    #save results for type 1 error rate betaX=0 and power betaX>0
   matR <- matrix(0,ncol=9,nrow=length(betaX))
    colnames(matR) <- c("case1","case2","case3","Z+","Steiger","MR","corGX","corGY","corXY")
  
################################################################################
# cycle through the sideltalations
################################################################################  
    #cycle through the sideltalations
    for(ii in 1:nSim){
      printCut<-500
      if(floor(ii/printCut)==ceiling(ii/printCut)){print(paste(ii,"of",nSim,"simulations"))}
      
################################################################################
# cycle through values of betaX
################################################################################ 
      #cycle through values of betaX
      for(bX in 1:length(betaX)){
        
################################################################################
# generate data
################################################################################                            
        # Generate G (SNPs)
		nSNP<-length(MAF) # number of SNPs
        G<-matrix(0,nrow=n,ncol=nSNP)
		 for(mm in 1:nSNP){
		 G[,mm] <- rbinom(n,2,MAF[mm])
		 }
		
		# Generate unmeasured confounder U
         U<-rnorm(n,mean=meanU,sd=sqrt(varU))
    
        # Generate X (exposure/ intermediate phenotype)
        gammaG<-matrix(gammaG,nrow=length(gammaG),ncol=1)
        
        # Unmeasured confouding True or False
        if(unmeasuredConfounding==F){
        Xtrue <- rnorm(n,(gamma0 + G%*%gammaG),sqrt(varX)) } 
        if(unmeasuredConfounding==T){
        Xtrue <- rnorm(n,(gamma0 + G%*%gammaG+ gammaU*U),sqrt(varX)) } 
        
        # Measurment error True or False
		if(measurementError==T ){
			X<-rnorm(n,(delta0+deltaX*Xtrue),sqrt(varME))}
		if(measurementError==F ){
			X<-Xtrue}
	 
	 Xtrue<-(Xtrue- mean(Xtrue))/sd(Xtrue)
	 
		# Generate Y (outcome)
        if(pleiotropy==T & unmeasuredConfounding==F){
        	Y <- rnorm(n,(beta0 + betaX[bX]*Xtrue +G%*%betaG),sqrt(varY))}
        	if(pleiotropy==T & unmeasuredConfounding==T){
        	Y <- rnorm(n,(beta0 + betaX[bX]*Xtrue +G%*%betaG +betaU*U),sqrt(varY))}
        if(pleiotropy==F & unmeasuredConfounding==T){
         	Y <- rnorm(n,(beta0 + betaX[bX]*Xtrue +betaU*U),sqrt(varY))}
        if(pleiotropy==F & unmeasuredConfounding==F){
         	Y <- rnorm(n,(beta0 + betaX[bX]*Xtrue),sqrt(varY))}
	 
	Y<-(Y-mean(Y))/sd(Y)
        
        #Correlation with first SNP        
        matR[bX,"corGY"]<-matR[bX,"corGY"]+cor(G[,1],Y)
        matR[bX,"corGX"]<-matR[bX,"corGX"]+cor(G[,1],X)
        matR[bX,"corXY"]<-matR[bX,"corXY"]+cor(X,Y)

###############################################################
# 1 SNP based on paper
###############################################################
if(nSNP==1){
#xhat=BetaG_Hat*G where x=alphaG+BetaG*G+epsilonG where epsilonG~N(0,sigma^2)
xhat<-G*lm(X~G)$coef["G"]

#p-value for MR from y=betaMR*xHat+epsilonMR 
pMR<-summary(lm(Y~xhat))$coef["xhat",4]

# H0: BetaMR=0
if(pMR<sig.level){matR[bX,"MR"]<-matR[bX,"MR"]+1}

#william test of correlation from Steiger paper
Stest<-r.test(n, r12 = abs(cor(G,X)), r13 = abs(cor(G,Y)),r23=abs(cor(X,Y)),twotailed = TRUE)
Z<-Stest$t
pSteiger<-Stest$p

if(Z>0){matR[bX,"Z+"]<-matR[bX,"Z+"]+1}
if(pSteiger<sig.level){matR[bX,"Steiger"]<-matR[bX,"Steiger"]+1}
       
# case 1: X->Y if pSteiger<alpha and pMR<alpha and Z>0
# case 2: X<-Y if pSteiger<alpha and pMR<alpha and Z<0
# case 3: neither if pSteiger>alpha or pMR>alpha 
if(pSteiger<sig.level & pMR<sig.level){
if(Z>0){matR[bX,"case1"]<-matR[bX,"case1"]+1}
if(Z<0){matR[bX,"case2"]<-matR[bX,"case2"]+1}
}else{matR[bX,"case3"]<-matR[bX,"case3"]+1}

} #end if 1 SNP
################################################################################
# More than 1 SNP based on code from mr_steiger function in from TwoSampleMR package
# https://github.com/MRCIEU/TwoSampleMR
################################################################################
#more than 1 SNP
if(nSNP>1){
	
#xhat=BetaG_Hat*G where x=alphaG+BetaG*G+epsilonG where epsilonG~N(0,sigma^2)
xhat<-G%*%lm(X~G)$coef[2:(ncol(G)+1)]

#p-value for MR from y=betaMR*xHat+epsilonMR 
pMR<-summary(lm(Y~xhat))$coef["xhat",4]

# H0: BetaMR=0
if(pMR<sig.level){matR[bX,"MR"]<-matR[bX,"MR"]+1}
	
#Vector of absolute correlations for SNP-exposure
r_exp<-rep(0,nSNP)
for(re in 1:nSNP){r_exp[re]<-abs(cor(G[,re],X))}
        
#Vector of absolute correlations for SNP-outcome
r_out<-rep(0,nSNP)
for(ro in 1:nSNP){r_out[ro]<-abs(cor(G[,ro],Y))}
        
# Code from mr_steiger function
r_exp <- sqrt(sum(r_exp[!is.na(r_exp) | is.na(r_out)]^2))
r_out <- sqrt(sum(r_out[!is.na(r_exp) | is.na(r_out)]^2))

#william test of correlation from Steiger paper
Stest<-r.test(n, r12 = r_exp, r13 = r_out,r23=abs(cor(X,Y)),twotailed = TRUE)
Z<-Stest$t
pSteiger<-Stest$p

if(Z>0){matR[bX,"Z+"]<-matR[bX,"Z+"]+1}
if(pSteiger<sig.level){matR[bX,"Steiger"]<-matR[bX,"Steiger"]+1}
        
# case 1: X->Y if pSteiger<alpha and pMR<alpha and Z>0
# case 2: X<-Y if pSteiger<alpha and pMR<alpha and Z<0
# case 3: neither if pSteiger>alpha or pMR>alpha 
if(pSteiger<sig.level & pMR<sig.level){
if(Z>0){matR[bX,"case1"]<-matR[bX,"case1"]+1}
if(Z<0){matR[bX,"case2"]<-matR[bX,"case2"]+1}
}else{matR[bX,"case3"]<-matR[bX,"case3"]+1}

}# end if more than 1 SNP
###############################################################


################################################################################
# end loops
################################################################################    
      }#beta loop
    }#sim loop

################################################################################
# Plot and save results
################################################################################    
    mat_total <- matR/nSim
    
    if(plot.pdf){
    	
    	################### 
    	#case plot
    ################### 
      pdf(paste(plot.name,"Cases.pdf",sep=""))
      plot(-2,-2,xlim=c(min(betaX),max(betaX)),ylim=c(0,1.1),main="",xlab=expression(beta[X]),ylab="Proportion of Simulations")
      lines(betaX,mat_total[,"case1"],col=1,pch=1,lty=1,type="b")
      lines(betaX,mat_total[,"case2"],col=2,pch=2,lty=2,type="b")
      lines(betaX,mat_total[,"case3"],col=3,pch=3,lty=3,type="b")
      legend("topleft",c("case 1: X->Y","case 2: X<-Y","case 3: inconclusive"),col=c(1:3),pch=c(1:3),lty=c(1:3),horiz=T,x.intersp=0)
      dev.off()
      
      	################### 
      # p-value steiger plot
      	################### 
      
      pdf(paste(plot.name,"Components.pdf",sep=""))
      plot(-2,-2,xlim=c(min(betaX),max(betaX)),ylim=c(0,1.1),main="",xlab=expression(beta[X]),ylab="Proportion of Simulations")
      lines(betaX,mat_total[,"Z+"],col=1,pch=1,lty=1,type="b")
      lines(betaX,mat_total[,"Steiger"],col=2,pch=2,lty=2,type="b")
      lines(betaX,mat_total[,"MR"],col=3,pch=3,lty=3,type="b")
      legend("topleft",c("Test statistic >0","Steiger P < alpha","MR P < alpha"),
      col=c(1:3),pch=c(1:3),lty=c(1:3),horiz=T,x.intersp=0)
      dev.off()
      
      	################### 	
      #correlation plot
      	################### 
       pdf(paste(plot.name,"Correlation.pdf",sep=""))
      plot(-2,-2,xlim=c(min(betaX),max(betaX)),ylim=c(0,1.1),main="",xlab=expression(beta[X]),ylab="Average")
      lines(betaX,mat_total[,"corGX"],col=1,pch=1,lty=1,type="b")
      lines(betaX,mat_total[,"corGY"],col=2,pch=2,lty=2,type="b")
      lines(betaX,mat_total[,"corXY"],col=3,pch=3,lty=3,type="b")
      legend("topleft",c("Correlation G1 & X", "Correlation G1 & Y","Correlation X & Y"),col=c(1:3),pch=c(1:3),lty=c(1:3),horiz=T,x.intersp=0)
      dev.off()
      
    }
    
    # Print out the matrix but use list
     colnames(matR) <- c("case1:X->Y","case2:Y->X","case3:neither","Z+","Steiger","MR","corG1X","corG1Y","corXY")
  
    listM<-list(mat_total)
    names(listM)<-c("matrix")
   listM   
  }
