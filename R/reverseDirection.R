reverseDirection <-
function(nSim=1000,n=100,MAF=0.5,gamma0=0,gammaG=0.2,varX=1,
measurementError=F,delta0=0,deltaX=1,varME=1,
beta0=0,betaX=seq(from=0,to=1,length.out=4),pleiotropy=F,betaG=1,varY=0.2,
sig.level=0.05,SEED=1,plot.pdf=T,plot.name="reverseDirection.pdf"){
 
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
if(n<0 | n==0 | floor(n)!=ceiling(n) ){stop("n deltast be an integer greater than 0")}
if(nSim<0 | nSim==0 | floor(nSim)!=ceiling(nSim) ){stop("nSim deltast be an integer greater than 0")}
if(SEED<0 | SEED==0 | floor(SEED)!=ceiling(SEED) ){stop("SEED deltast be an integer greater than 0")}
if(length(MAF)!=1| min(MAF)<0 | max(MAF)>1){stop("MAF deltast be a single value greater than 0 and less than 1")}
if(length(sig.level)!=1| sig.level<0 | sig.level>1){stop("sig.level deltast be a single value greater than 0 and less than 1")}
if(length(unique(betaX))<2){stop("betaX deltast be a vector with at least two values")}
if(!varX>0){stop("varX deltast be greater than 0")}
if(!varY>0){stop("varY deltast be greater than 0")}
if(!varME>0){stop("varME deltast be greater than 0")} 

if(length(MAF)!=1){stop("length(MAF) deltast equal 1")}
if(length(gamma0)!=1){stop("length(gamma0) deltast equal 1")}
if(length(gammaG)!=1){stop("length(gammaG) deltast equal 1")}
if(length(varX)!=1){stop("length(varX) deltast equal 1")}

if(length(delta0)!=1){stop("length(delta0) deltast equal 1")}
if(length(deltaX)!=1){stop("length(deltaX) deltast equal 1")}
if(length(varME)!=1){stop("length(varME) deltast equal 1")}

if(length(beta0)!=1){stop("length() deltast equal 1")}
if(length(betaG)!=1){stop("length() deltast equal 1")}
if(length(varY)!=1){stop("length(varY) deltast equal 1")}
  
if(plot.pdf==T|plot.pdf=="T"|plot.pdf=="t"|plot.pdf=="True"|plot.pdf=="true"){plot.pdf<-T}  
if(plot.pdf==F|plot.pdf=="F"|plot.pdf=="f"|plot.pdf=="False"|plot.pdf=="false"){plot.pdf<-F} 
if(plot.pdf!=T & plot.pdf!=F){stop("plot.pdf deltast equal True or False. Please not that R is case sensitive.")}

if(pleiotropy==T|pleiotropy=="T"|pleiotropy=="t"|pleiotropy=="True"|pleiotropy=="true"){pleiotropy<-T}  
if(pleiotropy==F|pleiotropy=="F"|pleiotropy=="f"|pleiotropy=="False"|pleiotropy=="false"){pleiotropy<-F} 
if(pleiotropy!=T & pleiotropy!=F){stop("pleiotropy deltast equal True or False. Please not that R is case sensitive.")}

if(measurementError==T|measurementError=="T"|measurementError=="t"|measurementError=="True"|measurementError=="true"){measurementError<-T}  
if(measurementError==F|measurementError=="F"|measurementError=="f"|measurementError=="False"|measurementError=="false"){measurementError<-F} 
if(measurementError!=T & measurementError!=F){stop("measurementError deltast equal True or False. Please not that R is case sensitive.")}

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
      printCut<-100
      if(floor(ii/printCut)==ceiling(ii/printCut)){print(paste(ii,"of",nSim,"simulations"))}
      
################################################################################
# cycle through values of betaX
################################################################################ 
      #cycle through values of betaX
      for(bX in 1:length(betaX)){
        
################################################################################
# generate data
################################################################################                            
        # Generate G (1 SNP)
		 G <- rbinom(n,2,MAF)

        # Generate X (exposure/ intermediate phenotype)
        Xtrue <- rnorm(n,(gamma0 + gammaG*G),sqrt(varX))
		if(measurementError==T){X<-rnorm(n,(delta0+deltaX*Xtrue),sqrt(varME))}
		if(measurementError==F){X<-Xtrue}
		
		# Generate Y (outcome)
        if(pleiotropy==T){Y <- rnorm(n,(beta0 + betaX[bX]*Xtrue +betaG*G),sqrt(varY))}
         if(pleiotropy==F){Y <- rnorm(n,(beta0 + betaX[bX]*Xtrue),sqrt(varY))}
                
        matR[bX,"corGY"]<-matR[bX,"corGY"]+cor(G,Y)
        matR[bX,"corGX"]<-matR[bX,"corGX"]+cor(G,X)
        matR[bX,"corXY"]<-matR[bX,"corXY"]+cor(X,Y)

################################################################################
# MR p-value
################################################################################ 
#xhat=BetaG_Hat*G where x=alphaG+BetaG*G+epsilonG where epsilonG~N(0,sigma^2)
xhat<-G*lm(X~G)$coef["G"]

#p-value for MR from y=betaMR*xHat+epsilonMR 
pMR<-summary(lm(Y~xhat))$coef["xhat",4]

# H0: BetaMR=0
if(pMR<sig.level){matR[bX,"MR"]<-matR[bX,"MR"]+1}

################################################################################
#Steiger correlation
################################################################################
#william test of correlation from Steiger paper
Stest<-r.test(n, r12 = abs(cor(G,X)), r13 = abs(cor(G,Y)),r23=abs(cor(X,Y)),twotailed = TRUE)
Z<-Stest$t
pSteiger<-Stest$p

if(Z>0){matR[bX,"Z+"]<-matR[bX,"Z+"]+1}
if(pSteiger<sig.level){matR[bX,"Steiger"]<-matR[bX,"Steiger"]+1}
#if(abs(cor(G,X)) > abs(cor(G,Y))){matR[bX,"CorrectDirection"]<-matR[bX,"CorrectDirection"]+1}
        
# case 1: X->Y if pSteiger<alpha and pMR<alpha and Z>0
# case 2: X<-Y if pSteiger<alpha and pMR<alpha and Z<0
# case 3: neither if pSteiger>alpha or pMR>alpha 
if(pSteiger<sig.level & pMR<sig.level){
if(Z>0){matR[bX,"case1"]<-matR[bX,"case1"]+1}
if(Z<0){matR[bX,"case2"]<-matR[bX,"case2"]+1}
}else{matR[bX,"case3"]<-matR[bX,"case3"]+1}
         
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
      pdf(plot.name)
      plot(-2,-2,xlim=c(min(betaX),max(betaX)),ylim=c(0,1.1),main="",xlab=expression(beta[X]),ylab="Proportion of Simulations")
      lines(betaX,mat_total[,"case1"],col=1,pch=1,lty=1,type="b")
      lines(betaX,mat_total[,"case2"],col=2,pch=2,lty=2,type="b")
      lines(betaX,mat_total[,"case3"],col=3,pch=3,lty=3,type="b")
      legend("topleft",c("case 1: X->Y","case 2: X<-Y","case 3: inconclusive"),col=c(1:3),pch=c(1:3),lty=c(1:3),horiz=T,x.intersp=0)
      dev.off()
    }
    
    # Print out the matrix but use list
    listM<-list(mat_total)
    names(listM)<-c("matrix")
   listM   
  }
