reverseDirection <-
function(nSim=1000, n=100,
  nSNP=1,MAF=rep(0.5,1),gamma0=0,gammaX=rep(0.2,1),varM=1, 
  beta0=0,betaM=seq(from=0,to=1,length.out=4),varY=0.2,
  delta0=0,deltaX=rep(1,1),varU=1, gammaU=0, betaU=1,Uconfounder=T
  ,alpha=0.05,SEED=1,plot.pdf=T,plot.name="reverseDirection.pdf"){
    
    # load library for mr_steiger
    library(TwoSampleMR)
    
    # Set the seed
    set.seed(SEED)
    
    ################################################################################
    # Error checks
    ################################################################################
    if(n<0 | n==0 | floor(n)!=ceiling(n) ){stop("Error: n must be an integer greater than or equal to 1")}
    if(nSNP<0 | nSNP==0 | floor(nSNP)!=ceiling(nSNP) ){stop("Error: nSNP must be an integer greater than or equal to 1")}
    if(length(MAF)!=nSNP){stop("Error: nSNP must equal the length of MAF")}
    if(min(MAF)<0 | max(MAF)>1){stop("Error: MAF must be greater than 0 and less than 1")}
    if(length(unique(betaM))<2){stop("Error: betaM must be a vector with at least two values")}
    if(!varM>0){stop("Error: varM must be greater than 0")}
    if(!varY>0){stop("Error: varY must be greater than 0")}
    if(length(gammaX)!=nSNP){stop("Error: Length of gammaX does not equal nSNP")}
    if(length(deltaX)!=nSNP){stop("Error: Length of deltaX does not equal nSNP")}
  
    ################################################################################
    #matrix to save MR steiger
    ################################################################################
    #save results for type 1 error rate betaM=0 and power betaM>0
   matR <- matrix(0,ncol=10,nrow=length(betaM))
    colnames(matR) <- c("McausesY","McausesYadj","CorrectDirection","CorrectDirectionAdj","SteigerTest","SteigerTestAdj","SensitivityRatio","corX1M","corX1Y","corMY")
  
    ################################################################################
    # cycle through the simulations
    ################################################################################  
    #cycle through the simulations
    for(ii in 1:nSim){
      pcut<-50
      if(floor(ii/pcut)==ceiling(ii/pcut)){print(paste(ii,"of",nSim,"simulations"))}
      
      ################################################################################
      # cycle through values of betaM
      ################################################################################ 
      #cycle through values of betaM
      for(bM in 1:length(betaM)){
        
        ################################################################################
        # generate data
        ################################################################################                            
        # Generate X (SNPs)
        X <- matrix(NA,nrow=n,ncol=nSNP)
        # Create vector of SNPs
        for(jj in 1:nSNP){
          X[,jj] <- rbinom(n,2,MAF[jj])
        }
        
        # Generate U, M and Y
        if(Uconfounder==T){
        	U<-rnorm(n,delta0+ X%*%deltaX,sqrt(varU))
        	M <- rnorm(n,gamma0 + X%*%gammaX +gammaU*U,sqrt(varM))
        Y <- rnorm(n,beta0 + betaM[bM]*M +betaU*U,sqrt(varY))
        }
        if(Uconfounder==F){
        	U<-rnorm(n,delta0,sqrt(varU))
        	M <- rnorm(n,gamma0 + X%*%gammaX,sqrt(varM))
        Y <- rnorm(n,beta0 + betaM[bM]*M,sqrt(varY))
        	}
        
        matR[bM,"corX1Y"]<-matR[bM,"corX1Y"]+cor(X[,1],Y)
        matR[bM,"corX1M"]<-matR[bM,"corX1M"]+cor(X[,1],M)
        matR[bM,"corMY"]<-matR[bM,"corMY"]+cor(M,Y)
        
        ################################################################################
        # Get p-values from Y~X and M~X
        ################################################################################ 
        #Vector of p-values of SNP-exposure
        p_exp<-rep(0,nSNP)
        for(pe in 1:nSNP){
          p_exp[pe]<-summary(lm(M~X[,pe]))$coef[2,4]
        }
        
        #Vector of p-values of SNP-outcome
        p_out<-rep(0,nSNP)
        for(po in 1:nSNP){
          p_out[po]<-summary(lm(Y~X[,po]))$coef[2,4]
        }
        
        #Sample sizes for p_exp, p_out
        n_exp<-n
        n_out<-n
        
        #Vector of absolute correlations for SNP-exposure
        r_exp<-rep(0,nSNP)
        for(re in 1:nSNP){
          r_exp[re]<-abs(cor(X[,re],M))
        }
        
        #Vector of absolute correlations for SNP-outcome
        r_out<-rep(0,nSNP)
        for(ro in 1:nSNP){
          r_out[ro]<-abs(cor(X[,ro],Y))
        }
        
        #r_xxo Measurememt precision of exposure, default 1
        #r_yyo Measurement precision of outcome, default 1
        
        ################################################################################
        # MR Steiger correct way
        ################################################################################ 
        #A statistical test for whether the assumption that exposure causes outcome is valid
        mrs<-mr_steiger(p_exp, p_out, n_exp, n_out, r_exp, r_out, r_xxo = 1, r_yyo = 1)
        
        #correct_causal_direction: TRUE/FALSE 
        if(mrs$correct_causal_direction==TRUE){matR[bM,"CorrectDirection"]<-matR[bM,"CorrectDirection"]+1}

        #correct_causal_direction_adj: TRUE/FALSE, direction of causality for given measurement error parameters 
        if(mrs$correct_causal_direction_adj==TRUE){matR[bM,"CorrectDirectionAdj"]<-matR[bM,"CorrectDirectionAdj"]+1}

        #steiger_test: p-value for inference of direction 
        if(mrs$steiger_test<alpha){matR[bM,"SteigerTest"]<-matR[bM,"SteigerTest"]+1}

        #steiger_test_adj: p-value for inference of direction of causality for given measurement error parameters 
        if(mrs$steiger_test_adj<alpha){matR[bM,"SteigerTestAdj"]<-matR[bM,"SteigerTestAdj"]+1}

        #sensitivity_ratio: Ratio of vz1/vz0. Higher means inferred direction is less susceptible to measurement error - 
        matR[bM,"SensitivityRatio"]<-matR[bM,"SensitivityRatio"]+mrs$sensitivity_ratio
        
         # correct direction and Steiger test reject
        if(mrs$correct_causal_direction==TRUE & mrs$steiger_test<alpha){
        	matR[bM,"McausesY"]<-matR[bM,"McausesY"]+1}
        if(mrs$correct_causal_direction_adj==TRUE & mrs$steiger_test_adj<alpha){
        	matR[bM,"McausesYadj"]<-matR[bM,"McausesYadj"]+1}
         
        ################################################################################
        # end loops
        ################################################################################    
      }#beta loop
    }#sim loop
    
    mat_total <- matR/nSim
    
    if(plot.pdf){
      pdf(plot.name)
      plot(-2,-2,xlim=c(min(betaM),max(betaM)),ylim=c(0,1),main="",xlab=expression(beta),ylab="Correct Direction (% of simulations)")
      lines(betaM,mat_total[,"CorrectDirection"],col=1,pch=1,type="b")
      dev.off()
    }
    
    # Print out the matrix but use list
    listM<-list(mat_total)
    names(listM)<-c("matrix")
   listM
    
  
    
  }
