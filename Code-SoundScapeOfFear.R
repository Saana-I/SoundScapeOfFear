Sys.setenv(TZ='GMT')

## R tools

  library(geepack)
  library(mvtnorm)
  library(MuMIn)

  logit <- function(p) {return(log(p/(1-p)))}
  ilogit <- function(x) {return(exp(x)/(exp(x)+1))}

## Load data
  
  load("Data-SoundscapeOfFear.Rd")
  
## Edit modelled variables
  
  etab$Species <- factor(etab$Species, levels=c("SW", "PW", "HW", "BW"))
  etab$Species2 <- as.character(etab$Species)
  etab$Species2[etab$Species2=="BW"] <- "HW"
  etab$Species2 <- factor(etab$Species2, levels=c("SW", "PW", "HW"))
  
  etab$NS <- as.numeric(etab$Session=="No-sonar")
  etab$LFAS <- as.numeric(etab$Session=="LFAS")    # LFAS = pulsed active sonar (1-2 kHz in SW, PW, HW, 3-4 kHz in BW)
  etab$PB_BBN <- as.numeric(etab$Session=="PB-BBN")
  etab$PB_KWM <- as.numeric(etab$Session=="PB-mammal")
  
  etab$ind_num <- as.numeric(as.factor(etab$ind))
  
  etab$RI.pre <- (etab$Foraging1-etab$pre.Foraging1)/etab$abl.Foraging1
  etab$RI.abl <- (etab$Foraging1-etab$abl.Foraging1)/etab$abl.Foraging1
  etab$RI.kw <- (etab$akw.Foraging1-etab$abl.Foraging1)/etab$abl.Foraging1
  
  
## Fig S2 - individual responses
  
  etab$pch <- NA
  etab$pch[etab$Species=="PW"] <- 1
  etab$pch[etab$Species=="SW"] <- 2
  etab$pch[etab$Species=="HW"] <- 22
  etab$pch[etab$Species=="BW"] <- 23
  
  etab$pch[etab$Session!="Baseline" & etab$Species=="PW"] <- 16 # circles
  etab$pch[etab$Session!="Baseline" & etab$Species=="SW"] <- 17 # triangles
  etab$pch[etab$Session!="Baseline" & etab$Species=="HW"] <- 15 # squares
  etab$pch[etab$Session!="Baseline" & etab$Species=="BW"] <- 18 # diamonds
  
  tpval <- 0.5
  etab$col <- adjustcolor("black",alpha.f=tpval)
  etab$col[etab$Session=="LFAS"] <- adjustcolor("Blue",alpha.f=tpval)
  etab$col[etab$Session=="PB-mammal"] <- adjustcolor("Red",alpha.f=tpval)
  etab$col[etab$Session=="PB-fish"] <- adjustcolor("Orange",alpha.f=tpval)
  
  par(mfrow=c(1,1), mar=c(4,4,2,2))
  
  Session_num <- rep(0, length(etab$Session))
  Session_num[etab$Session=="Baseline"] <- 1.5
  Session_num[etab$Session=="No-sonar"] <- 1
  Session_num[etab$Session=="LFAS"] <- 0.5
  Session_num[etab$Session=="PB-mammal"] <- -0.5
  Session_num[etab$Session=="PB-fish"] <- -1

  x <- etab$RI.kw+Session_num/(80/(max(etab$RI.kw)-min(etab$RI.kw)))
  plot(-x*100, -etab$RI.pre*100, col=etab$col, pch=etab$pch,
       xlab="% Reduction during KW-mammal",
       ylab="% Reduction during all exposures", main="")
  abline(0,1, col=1, lty=2)
  grid(col="grey")
 
   
## Model fitting
  
  # Sample size for each exposure session
    table(etab$Session)
    # Baseline  No-sonar LFAS   PB-BBN    PB-mammal   PB-fish 
    # 202       24       26     13        18          4 
  
  # Remove data from during playbacks to fish-feeding killer whales
    etab <- etab[etab$Session!="PB-fish",]
  
  # Without species interactions
  fit00 <- geeglm(Foraging1 ~ Species + NS + LFAS + PB_BBN + PB_KWM, 
                  weights=overlap,id=ind_num, family=binomial(link="logit"),
                  corstr="independence",
                  data=etab, na.action = "na.fail")
  summary(fit00)
  anova(fit00)
  
  # with species interactions
  fit0 <- geeglm(Foraging1 ~ Species + NS + LFAS + PB_BBN + PB_KWM + 
                   Species2:NS + Species:LFAS + Species2:PB_BBN + Species:PB_KWM, 
                 weights=overlap,id=ind_num, family=binomial(link="logit"),
                 corstr="independence",
                 data=etab, na.action = "na.fail")
  summary(fit0)
  
  # Suppelementary Table 3
  # TableS3 <- rbind(data.frame(summary(fit00)$coefficients),data.frame(summary(fit0)$coefficients))
  
  
## Predict from both models
  
  speciesCode0 <- names(tapply(etab$Foraging1, etab$Species, unique))
  speciesCode2 <- speciesCode0
  speciesCode2[speciesCode2=="BW"] <- "HW"
  
  # Without species interactions
  
  # Reduction in odds during LFAS
    1-exp(summary(fit00)$coefficients$Estimate[6]) # 0.7952
  # Reduction in odds during KW-mammal
    1-exp(summary(fit00)$coefficients$Estimate[8]) # 0.8168
  # Percentage reduction during during sound exposures
    preddata1 <- data.frame(Species = speciesCode0)
    preddata1$NS <- 0; preddata1$LFAS <- 0; preddata1$PB_BBN <- 0; preddata1$PB_KWM <- 0; preddata1$overlap <- 1;
    preds_baseline <- predict(fit00, type="response", newdata=preddata1) 
    preddata1$PB_KWM <- 0; preddata1$LFAS <- 1;
    preds_LFAS <- predict(fit00, type="response", newdata=preddata1)
    preddata1$PB_KWM <- 1; preddata1$LFAS <- 0;
    preds_KWM <- predict(fit00, type="response", newdata=preddata1)
    (preds_LFAS-preds_baseline)/preds_baseline # -0.535 -0.766 -0.764 -0.678 
    (preds_KWM-preds_baseline)/preds_baseline # -0.570 -0.790 -0.788 -0.707 
    mean((preds_LFAS-preds_baseline)/preds_baseline) # -0.686
    mean((preds_KWM-preds_baseline)/preds_baseline)  # -0.714
  
  # With species interactions
  
  preddata1 <- data.frame(Species=speciesCode0[4:1], Species2=speciesCode2[4:1])
  preddata1$NS <- 0
  preddata1$LFAS <- 0
  preddata1$PB_BBN <- 0
  preddata1$PB_KWM <- 0
  preddata1$preds_baseline <- predict(fit0, type="response", newdata=preddata1)
  
  preddata1$LFAS <- 1
  preddata1$preds_LFAS <- predict(fit0, type="response", newdata=preddata1)
  preddata1$LFAS <- 0
  
  preddata1$NS <- 1
  preddata1$preds_NS <- predict(fit0, type="response", newdata=preddata1)
  preddata1$NS <- 0
  
  preddata1$PB_KWM <- 1
  preddata1$preds_KWM <- predict(fit0, type="response", newdata=preddata1)
  preddata1$PB_KWM <- 0
  
  preddata1$PB_BBN <- 1
  preddata1$preds_BBN <- predict(fit0, type="response", newdata=preddata1)
  preddata1$PB_BBN <- 0
  
  preddata1$RI_LFAS <- (preddata1$preds_LFAS-preddata1$preds_baseline)/preddata1$preds_baseline
  preddata1$RI_NS <- (preddata1$preds_NS-preddata1$preds_baseline)/preddata1$preds_baseline
  
  preddata1$RI_KWM <- (preddata1$preds_KWM-preddata1$preds_baseline)/preddata1$preds_baseline
  preddata1$RI_BBN <- (preddata1$preds_BBN-preddata1$preds_baseline)/preddata1$preds_baseline
  

  # bootstrap to obtain confidence intervals
  
  P1 <- matrix(NA, 1000, length(preddata1$LFAS))
  P2 <- P1
  P3 <- P1
  P4 <- P1
  
  source("Code/rmvnorm_fun.R")
  modelrandom <- function(fit0) {
    rfit <- fit0
    rcoefs <- rmvnorm2(1, coef(fit0), summary(fit0)$cov.scaled)
    rcoefs <- as.numeric(rcoefs)
    names(rcoefs) <- names(fit0$coefficients)
    rfit$coefficients <- rcoefs
    return(rfit)
  }
  
  for(j in 1:1000) {
    
    
    rfit0 <- modelrandom(fit0)
    
    # Baseline prediction
    preds1_baseline <- predict(rfit0, newdata=preddata1, type="response")
    
    # LFAS prediction
    preddata1$LFAS <- 1
    preds1_LFAS <- predict(rfit0, newdata=preddata1, type="response")
    preddata1$LFAS <- 0
    
    # NS prediction
    preddata1$NS <- 1
    preds1_NS <- predict(rfit0, newdata=preddata1, type="response")
    preddata1$NS <- 0
    
    # PB_KWM prediction
    preddata1$PB_KWM <- 1
    preds1_PB_KWM <- predict(rfit0, newdata=preddata1, type="response")
    preddata1$PB_KWM <- 0
    
    # PB_BBN prediction
    preddata1$PB_BBN <- 1
    preds1_PB_BBN <- predict(rfit0, newdata=preddata1, type="response")
    preddata1$PB_BBN <- 0
    
    # Store RI
    P1[j,] <- (preds1_LFAS-preds1_baseline)/preds1_baseline
    P2[j,] <- (preds1_PB_KWM-preds1_baseline)/preds1_baseline
    P3[j,] <- (preds1_NS-preds1_baseline)/preds1_baseline
    P4[j,] <- (preds1_PB_BBN-preds1_baseline)/preds1_baseline
    
  }
  
  preddata1$RI_LFAS.u <- apply(P1,2,quantile,0.025)
  preddata1$RI_LFAS.l <- apply(P1,2,quantile,0.975)
  
  preddata1$RI_KWM.u <- apply(P2,2,quantile,0.025)
  preddata1$RI_KWM.l <- apply(P2,2,quantile,0.975)
  
  preddata1$RI_NS.u <- apply(P3,2,quantile,0.025)
  preddata1$RI_NS.l <- apply(P3,2,quantile,0.975)
  
  preddata1$RI_BBN.u <- apply(P4,2,quantile,0.025)
  preddata1$RI_BBN.l <- apply(P4,2,quantile,0.975)
  
 
   
## Fig 3
  
  par(mfrow=c(1,1), mar=c(5,5,0.5,0.5))
  
  mypch <- c(18,15,16,17)
  tpval <- 0.5
  
  mycol <- c(adjustcolor("Red",alpha.f=tpval),   # BW
             adjustcolor("Orange",alpha.f=tpval), # HW
             adjustcolor("Cyan",alpha.f=tpval), # PW
             adjustcolor("Blue",alpha.f=tpval)) # SW
  mycol2 <- c("Red", "Orange", "Cyan", "Blue")

  etab$col2 <- NA
  etab$col2[etab$Species=="BW"] <-  adjustcolor("Red",alpha.f=tpval)
  etab$col2[etab$Species=="HW"] <- adjustcolor("Orange",alpha.f=tpval)
  etab$col2[etab$Species=="PW"] <- adjustcolor("Cyan",alpha.f=tpval)
  etab$col2[etab$Species=="SW"] <- adjustcolor("Blue",alpha.f=tpval)

  plot(1,1, col=NA, yaxt="n", xaxt="n",
       xlab="% Reduction during KW-mammal",
       ylab="% Reduction during 1-4 kHz sonar",
       ylim=c(-50,120),xlim=c(-50,120), main="")
  axis(1, at=c(-300,-250,-200,-150,-100,-50,0,50,100))
  axis(2, at=c(-300,-250,-200,-150,-100,-50,0,50,100))
  
  # Rug plot
  
  radd <- 4
  tBool <- etab$Session=="LFAS"
  tx <- rep(120, sum(tBool))
  ty <- -(etab$RI.abl[tBool])*100
  tx[ty > 99.9] <- tx[ty > 99.9] + runif(length(tx[ty > 99.9]),-1*radd,radd)
  ty[ty > 99.9] <- ty[ty > 99.9] + runif(length(ty[ty > 99.9]),-1*radd,radd)
  points(tx, ty, pch=etab$pch[tBool], 
         col=etab$col2[tBool], cex=1)
  
  tBool <- etab$Session=="PB-mammal"
  tx <- rep(120, sum(tBool))
  ty <- -(etab$RI.abl[tBool])*100
  tx[ty > 99.9] <- tx[ty > 99.9] + runif(length(tx[ty > 99.9]),-1*radd,radd)
  ty[ty > 99.9] <- ty[ty > 99.9] + runif(length(ty[ty > 99.9]),-1*radd,radd)
  #tvec <- -120+runif(sum(tBool),0,radd)
  points(ty, tx, 
         pch=etab$pch[tBool], 
         col=etab$col2[tBool], cex=1)
  
  abline(0,1, col="grey", lty=1)
  
  segments(x0=-preddata1$RI_KWM*100, 
           x1=-preddata1$RI_KWM*100, 
           y0=-preddata1$RI_LFAS.l*100, 
           y1=-preddata1$RI_LFAS.u*100, 
           lwd=2, col=mycol)
  
  segments(x0=-preddata1$RI_KWM.l*100, 
           x1=-preddata1$RI_KWM.u*100, 
           y0=-preddata1$RI_LFAS*100, 
           y1=-preddata1$RI_LFAS*100, 
           lwd=2, col=mycol)
  
  points(-preddata1$RI_KWM*100, -preddata1$RI_LFAS*100, 
         pch=mypch, col=mycol2, cex=1.5)
  
  
## Correlation between species response to sonar vs killer whale playbacks
  
  cor.test(preddata1$RI_LFAS, preddata1$RI_KWM)
  
  
## Predicting sonar responses from average species response to killer whale playbacks
  
  etab_temp <- etab[etab$Session!= "No-sonar" & etab$Session!= "PB-mammal" & etab$Session!= "PB-BBN",]
  
  # Equal response to sonar across the four species / study populations
  fit1 <- geeglm(Foraging1 ~ Species + LFAS, weights=overlap,
                 id=ind_num, family=binomial(link="logit"),
                 corstr="independence",
                 data=etab_temp, na.action = "na.fail")
  
  # Response to sonar mediated by response to killer whale playbacks
  fit2 <- geeglm(Foraging1 ~ Species + RI.kw:LFAS, weights=overlap,
                 id=ind_num, family=binomial(link="logit"),
                 corstr="independence",
                 data=etab_temp, na.action = "na.fail")
  
  # Response to sonar explained by a species-specific response
  fit3 <- geeglm(Foraging1 ~ Species + Species:LFAS, weights=overlap,
                 id=ind_num, family=binomial(link="logit"),
                 corstr="independence",
                 data=etab_temp, na.action = "na.fail")

  
  QIC(fit1,fit2, fit3) # lowest QIC model: fit2
  QICu(fit1,fit2, fit3) # lowest QICu model: fit2
  
  
  # Supplementary Table 4
  # rbind(data.frame(anova(fit1)), data.frame(anova(fit2)), data.frame(anova(fit3)))