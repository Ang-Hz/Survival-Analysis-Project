
###SURVIVAL DATA ANALYSIS 

# Packages
install.packages("arsenal") 
install.packages("survival") # our main package for survival analysis
install.packages("SurvRegCensCov")

# Libraries
library(arsenal) 
library(survival)
library(SurvRegCensCov)


# Load data
esrd <- read.csv("esrd.csv", sep=";")
head(esrd)

# Factor variables
esrd$sex <- factor(esrd$sex,levels = c("F","M"))

indVar <- c("copd","dm","hypert","heart","liver","neoplasia","vascular","dth")
for (j in indVar)
{
  esrd[esrd[,j]==0,j] <- "No"
  esrd[esrd[,j]=="1",j] <- "Yes"
  esrd[,j] <- factor(esrd[,j],levels = c("No","Yes"))
}
esrd$treatmnt[esrd$treatmnt==0] <- "conservative"
esrd$treatmnt[esrd$treatmnt==1] <- "hemodialysis"
esrd$treatmnt <- factor(esrd$treatmnt,levels = c("conservative","hemodialysis"))

# Labels using library(arsenal)
attr(esrd$id,"label") <- "Unique ID number" 
attr(esrd$age,"label") <- "Age in years"
attr(esrd$sex,"label") <- "Sex"
attr(esrd$copd,"label") <- "Chronic obstructive pulmonary disease"
attr(esrd$dm,"label") <- "Diabetes mellitus"
attr(esrd$hypert,"label") <- "Hypertension"
attr(esrd$heart,"label") <- "Heart disease"
attr(esrd$liver,"label") <- "Liver disease"
attr(esrd$neoplasia,"label") <- "Neoplasia"
attr(esrd$vascular,"label") <- "Vascular disease"
attr(esrd$treatmnt,"label") <- "Treatment modality"
attr(esrd$dth,"label") <- "Death (from any cause)"
attr(esrd$dthtime,"label") <- "Time to death or to the end of follow up"

#Replace zero follow-up times by 0.001 years for the analysis
summary(esrd$dthtime) # there are some zero times 
esrd$dthtime[esrd$dthtime==0] <- 0.001 


### Descriptive characteristics by group (Table 1)  ###

# using the arsenal library
tbl1 <- tableby(treatmnt ~ dth + dthtime + age + sex + copd + hypert + heart + 
                liver + neoplasia + vascular, data = esrd,
                control = tableby.control(total = T, numeric.stats = c("meansd"), 
                                          digits = 1,cat.stats = "countpct",
                                          numeric.test = "anova", cat.test = "chisq"))
summary(tbl1, text=TRUE)



### Relative risk ratio and hemodialysis effectiveness ###


xx <- prop.table(table(esrd$treatmnt,esrd$dth),1)
print(xx)
RR <- xx[2,2]/xx[1,2] # death risk in hemodialysis / death risk in conservative therapy
RR # relative risk
(1 - RR)*100 # Hemodialysis effectiveness 


### Incidence rate ratio and hemodialysis effectiveness ###

# Poisson regression
fitPoiss <- glm(I(dth=="Yes") ~ I(treatmnt=="hemodialysis") + offset(log(dthtime)),
                data = esrd,family = poisson)
summary(fitPoiss)

round(exp(coef(fitPoiss)[1] + coef(fitPoiss)[2]),2) # Incidence rate for hemodialysis
round(exp(coef(fitPoiss)[1]),2) # Incidence rate for conservative treatment
round(exp(coef(fitPoiss)[2]),2) # IRR point estimate
round(exp(confint(fitPoiss)[2,]),2) # conf. limits for IRR
round(100*(1-exp(coef(fitPoiss)[2])),0) # Hemodialysis effectiveness point estimate
rev(round(100*(1-exp(confint(fitPoiss)[2,])),0))  # Hemodialysis effectiveness CI



### Kaplan-Meier Survival Curves  ###

# library(survival)

fitKM <- survfit( Surv(dthtime,dth=="Yes") ~ treatmnt,data = esrd,conf.type = "log-log")

plot(fitKM,xlab = "Years of follow-up",ylab = "Survival probability",col = c("black","red"),
     conf.int = T)
    legend("topright",col = c("black","red"),legend = levels(esrd$treatmnt),lty = 1,bty = "n")
  logRank <- survdiff(Surv(dthtime,dth=="Yes")~treatmnt,data = esrd)
    print(logRank)
    plogRank <- 1 - pchisq(logRank$chisq, length(logRank$n) - 1)
  peto <- survdiff(Surv(dthtime,dth=="Yes")~treatmnt,data = esrd,rho = 1) 
    print(peto)
    ppeto <- 1 - pchisq(peto$chisq, length(peto$n) - 1)
  text(x = 13,y = 0.8,label = paste0("Logrank test, p = ",round(plogRank,3))) 
  text(x = 13,y = 0.75,label = paste0("Peto-Peto test, p = ",round(ppeto,3)))


### KM survival probabilities & conditional survival probabilities ###

fitKM <- survfit( Surv(dthtime,dth=="Yes") ~ treatmnt,data = esrd,conf.type = "log-log")
  summary(fitKM,times = c(0.5,1,3,5))

  fitKM_2 <- survfit( Surv(dthtime,dth=="Yes") ~ treatmnt,data = esrd,conf.type = "log-log",
                    subset = esrd$dthtime > 0.5)
summary(fitKM_2,times = 0.5 + c(0.5,1,1.5,2))


### Effect of age ###

fitAgeL <- coxph(Surv(dthtime,dth=="Yes")~age,data = esrd)

esrd$ageCat <- cut(esrd$age,breaks = c(-Inf,
                                       quantile(esrd$age,probs = c(0.25,0.50,0.75)),
                                       Inf))
fitAgeCat <- coxph(Surv(dthtime,dth=="Yes")~ageCat,data = esrd)
fitAgeQ <- coxph(Surv(dthtime,dth=="Yes")~age+I(age^2),data = esrd)
# Likelihoods
fitAgeL$loglik
fitAgeCat$loglik
fitAgeQ$loglik
# AIC
AIC(fitAgeL) # can also use extractAIC(fitAgeL)
AIC(fitAgeCat)
AIC(fitAgeQ)
# BIC
BIC(fitAgeL)
BIC(fitAgeCat)
BIC(fitAgeQ)

# Investigation of functional form using martingale residuals
esrd$mgNULL <- residuals(coxph(Surv(dthtime,dth=="Yes")~1,data = esrd),type = "martingale")
plot(mgNULL~age,data = esrd,ylab = "Martingale residual",
     xlab = "Age (years)")
lines(smooth.spline(esrd$age, esrd$mgNULL, df = 3), col = "red", lwd = 2)
# lowess is similar
esrd$mgNULL <- residuals(coxph( Surv(dthtime,dth=="Yes") ~ 1,data = esrd),
                         type = "martingale")
plot(mgNULL ~ age, data = esrd,ylab = "Martingale residual")
lines(lowess(esrd$age,esrd$mgNULL),col = "red")

# The Grambsch-Therneau  PH Test 
zphage <- cox.zph(fitAgeL,transform = "identity",) # use t as the function of time (Stata default)
print(zphage)
plot(zphage[1], lwd=2, col="red")

zphage <- cox.zph(fitAgeL) #KM estimator S^(t) as the function of time (default)
print(zphage)
plot(zphage[1], lwd=2, col="red") # small P value, but plot looks sufficiently linear


# Effect per 1 year increase
summary(fitAgeL)

# Effect per 10 years increase
esrd$age10 <- esrd$age/10 
fitAgeL10 <- coxph(Surv(dthtime,dth=="Yes")~age10,data = esrd)
summary(fitAgeL10)
#alternatively
fitAge <- coxph( Surv(dthtime,dth=="Yes") ~ I(age/10),data = esrd)
summary(fitAge)


### Univariate Cox model for treatment mode ###

fitCOX <- coxph( Surv(dthtime,dth=="Yes") ~ treatmnt,data = esrd)
summary(fitCOX)

exp(coef(fitCOX)) # hazard ratio
round(exp(confint(fitCOX)),2)

round(100*(1-exp(coef(fitCOX))),0) # Effectiveness %
rev(round(100*(1-exp(confint(fitCOX))),0))



### Proportional hazards plots and tests for treatment mode ### 

fitKM <- survfit( Surv(dthtime,dth=="Yes") ~ treatmnt,data = esrd,conf.type = "log-log")
fitCOX <- coxph( Surv(dthtime,dth=="Yes") ~ treatmnt,data = esrd)

# Log-log survival plot
plot(fitKM, fun = "cloglog", col = c("blue","red"), 
     xlab = "Follow-up time (years in log scale)", ylab = "Cumulative hazard of death (log scale)" )
legend("topleft",col = c("blue","red"),legend = levels(esrd$treatmnt),lty = 1,bty = "n")


# KM vs Cox
plot(fitKM, xlab = "Years of follow-up",ylab = "Survival probability",col = c("blue","red"))
xx <- survfit(fitCOX,newdata = data.frame(treatmnt = levels(esrd$treatmnt)))
lines(xx,col = c("blue","yellow"),lty = 2)
legend("topright",col = c("blue","red","black","yellow"),
       legend = paste0(c("Observed (KM)","Observed (KM)","Predicted (Cox PH)","Predicted (Cox PH)"),": ",
                       rep(levels(esrd$treatmnt),2)), lty = c(1,1,2,2),bty = "n")

# Weighted Schoenfeld residuals
zphage <- cox.zph(fitCOX,transform = "identity")
zphage
plot(zphage[1], lwd=2, col="red")
plot(zphage[1],xlim = c(0,1), lwd=2, col="red") # restrict the data to the first year 
#using the defualt
# Grambsch-Therneau  PH Test 
zphageb <- cox.zph(fitCOX)
zphageb
plot(zphageb[1], lwd=2, col="red")
plot(zphageb[1],xlim = c(0,1), lwd=2, col="red") # restrict the data to the first year 



### Treatment mode effect using a time-dependent variable ### 

esrd$trt <- esrd$treatmnt=="hemodialysis" # esrd$trt = 1*(esrd$treatmnt=="hemodialysis")

# linear-time interaction
fitcox_tvc <- coxph(Surv(dthtime,dth=="Yes") ~ trt + tt(trt),data = esrd,
                    tt = function(x,t,...){x*t})
summary(fitcox_tvc) # the model is HR(t) = exp(b1*trt + b2*trt*t)

fitcox_tvc$coefficients["trtTRUE"] # b1 coefficient
fitcox_tvc$coefficients["tt(trt)"] # b2 coefficient

# plot the hazard ratio estimates (assuming linear time interaction)
tm <- seq(min(esrd$dthtime), 3, by = 0.001) # time values
hrtm <- exp(fitcox_tvc$coefficients["trtTRUE"] + fitcox_tvc$coefficients["tt(trt)"] * tm)
plot(tm, hrtm, type = "l", col = "blue", lwd = 2,
     xlab = "Time since beginning of treatment (years)", 
     ylab = "Estimated HR (haemodialysis/conservative)",
     xlim =c(0, 2),      ylim=c(0,18),
     main = "1")
abline(h = 1, col = "black", lty = 2)

# convert the plot to show haemodialysis effectiveness %
hetm <- 100 * (1 - exp(fitcox_tvc$coefficients["trtTRUE"] + fitcox_tvc$coefficients["tt(trt)"] * tm))
plot(tm, hetm, type = "l", col = "blue", lwd = 2,
     xlab = "Time since beginning of treatment (years)", 
     ylab = "Estimated HR (haemodialysis/conservative)",
     xlim =c(0, 1.5),      ylim=c(-300,100),
     main = "2")
abline(h = 0, col = "black", lty = 2)

# log-time interaction

fitcox_tvc <- coxph(Surv(dthtime,dth=="Yes") ~ trt + tt(trt),data = esrd,
                    tt = function(x,t,...){x*log(t)})
summary(fitcox_tvc) # the model is HR(t) = exp(b1*trt + b2*trt*t)

fitcox_tvc$coefficients["trtTRUE"] # b1 coefficient
fitcox_tvc$coefficients["tt(trt)"] # b2 coefficient

# plot the hazard ratio estimates (assuming log-time interaction)
tm <- seq(min(esrd$dthtime), 3, by = 0.001) # time values
hrtm <- exp(fitcox_tvc$coefficients["trtTRUE"]) * tm^(fitcox_tvc$coefficients["tt(trt)"])
plot(tm, hrtm, type = "l", col = "blue", lwd = 2,
     xlab = "Time since beginning of treatment (years)", 
     ylab = "Estimated HR (haemodialysis/conservative)",
     xlim =c(0, 2),      ,
     main = "3")
abline(h = 1, col = "black", lty = 2)
# convert the plot to show haemodialysis effectiveness %
hetm <- 100 * (1 - exp(fitcox_tvc$coefficients["trtTRUE"]) * tm^(fitcox_tvc$coefficients["tt(trt)"]))
plot(tm, hetm, type = "l", col = "blue", lwd = 2,
     xlab = "Time since beginning of treatment (years)", 
     ylab = "Estimated HR (haemodialysis/conservative)",
     xlim =c(0, 2),      
     main = "4")
abline(h = 0, col = "black", lty = 2)


### Multivariable Cox model ###

esrd$trt = 1*(esrd$treatmnt=="hemodialysis")

# Main-effects MV model
fitMultMain <- coxph( Surv(dthtime,dth=="Yes") ~ trt + I(age/10) + sex
                  + copd + dm + hypert +  heart + liver + neoplasia + vascular,data = esrd)
summary(fitMultMain)

# Interactions, loop
varlist <- c("age","sex","copd","dm","hypert","heart","liver","neoplasia","vascular")
for (xvar in varlist) {
  formula_string <-paste0("Surv(dthtime,dth=='Yes') ~ trt : ", xvar, " + trt + age + sex + copd + dm + hypert +  heart + liver +neoplasia + vascular")

  fitMultInter1 <- coxph(as.formula(formula_string),data = esrd)
  print(paste("Results for predictor:", xvar))
  # print(summary(fitMultage)) # Wald test
  print(anova(fitMultMain,fitMultInter1) ) # LR test
}

# Main-effects + 4 interaction terms
fitMultInt <- coxph( Surv(dthtime,dth=="Yes") ~ trt + I(age/10) + sex
                      + copd + dm + hypert +  heart + liver + neoplasia + vascular
                      + trt:dm + trt:heart + trt:hypert + trt:liver
                      ,data = esrd)
summary(fitMultInt) # trt:dm, trt:liver not significant

# Main-effects + 2 interaction terms
fitMultInt <- coxph( Surv(dthtime,dth=="Yes") ~ trt + I(age/10) + sex
                        + copd + dm + hypert +  heart + liver + neoplasia + vascular
                        +trt:heart + trt:hypert 
                        ,data = esrd)
summary(fitMultInt)

# Cox-snell
esrd$death = 1*(esrd$dth=="Yes")
esrd$cs <- esrd$death - residuals(fitMultInt,type = "martingale")
fitcs <- survfit(Surv(cs,death) ~ 1,data = esrd)
plot(fitcs,fun = "cloglog",xlab = "Cox-Snell (log scale)",
     ylab = "Log cumulative hazard of Cox-Snell",conf.int = F,col = "blue")
xx <- range(esrd$cs)
lines(xx,log(xx),type = "l",col = "red")
legend("topleft",col = c("blue","red"),legend = c("Observed","Ideal"),lty = 1)

# Grambsch-Therneau PH Test
print(cox.zph(fitMultInt))

# test PH with tvc
fitMultInt_tvc <- coxph( Surv(dthtime,dth=="Yes") ~ tt(trt) + I(age/10) + sex
                         + copd + dm + hypert +  heart + liver + neoplasia + vascular
                         +trt:heart + trt:hypert ,
                         tt = function(x,t,...){x*t},
                      data = esrd)
summary(fitMultInt_tvc) 
fitMultInt_tvc$coefficients["tt(trt)"] # b coefficient for the time-by-treatment interaction 

# Address non-PH: fit piecewise Cox model, before and after 1 year, adjusting for covariates 
# time-stratified Cox model with time-by-treatment interaction
# assume time-fixed interactions with covariates
splidat1 <- survSplit(Surv(dthtime,dth=="Yes") ~ ., data = esrd, cut = c(1), episode = "interval")
fitMultInt_pcwseCox <- coxph(Surv(dthtime,event==1) ~ trt:strata(interval) 
                             + I(age/10) + sex + copd + dm + hypert +  heart + liver + neoplasia + vascular
                             +trt:heart + trt:hypert , 
                             data = splidat1) # 
summary(fitMultInt_pcwseCox)


### Weibull model ###

# assess fit of Weibull distribution
# construct a log-log plot using the KM estimates 
fitKM1 <- survfit( Surv(dthtime,dth=="Yes") ~ 1,data = esrd)
plot(fitKM1, fun = "cloglog", mark.time = F, conf.int = FALSE, main = "",
     xlab = "log analysis time (log years)", ylab = "Estimated log(-logS(t))")

# estimate λ and γ by linearly regressing
loglogS <- log(-log(fitKM1$surv))
logtm <- log(fitKM1$time)
linreg <- lm(loglogS ~ logtm,data = esrd)
summary(linreg)
coef(linreg)[2] # estimate of γ (shape)
exp(coef(linreg)[1]) # estimate of λ (scale)

# estimate λ and γ by MLE using survreg
fitWei1 = survreg(Surv(dthtime,dth=="Yes")  ~ 1 , dist = "weibull",data = esrd)
summary(fitWei1)
# need to extract from an AFT parameterisation
fitWei1$scale # AFT scale (σ) #  1/σ = γ
log(fitWei1$scale) # output shows its log : Log(scale)=0.1353 
1 / fitWei1$scale  #shape γ= 1/(AFT scale) = 1/σ
coef(fitWei1) # shown as intercept
exp(-coef(fitWei1)/fitWei1$scale) # logλ = -intercept/σ


# MV Weibull model
esrd$trt = 1*(esrd$treatmnt=="hemodialysis")

fitMultIntWei <- survreg(Surv(dthtime,dth=="Yes") ~ trt + I(age/10) + sex
                     + copd + dm + hypert +  heart + liver + neoplasia + vascular
                     + trt:heart + trt:hypert ,
                     dist = "weibull",data = esrd)
summary(fitMultIntWei)
# Coefficients in PH form
exp(-coef(fitMultIntWei)/fitMultIntWei$scale) # b(PH) = - b(AFT)/σ
exp(-confint(fitMultIntWei)/fitMultIntWei$scale)


# Define the two comparison groups (at1 and at2) for a covariate pattern
at1 <- data.frame(trt=0, age=74, sex="M", copd="Yes",dm="Yes",hypert="Yes",heart="Yes", 
                  liver="Yes", neoplasia="Yes", vascular="Yes")
at2 <- data.frame(trt=1, age=74, sex="M", copd="Yes",dm="Yes",hypert="Yes",heart="Yes", 
                  liver="Yes", neoplasia="Yes", vascular="Yes")
pct <- 1:98 / 100 # Generate percentiles to plot (1% to 99%)
ptime1 <- predict(fitMultIntWei, newdata = at1, type = "quantile", p = pct)
ptime2 <- predict(fitMultIntWei, newdata = at2, type = "quantile", p = pct)

# Predicted Weibull Survival Curves: Plot survival (1-pct) against predicted times
plot(ptime1, 1 - pct, type = "l", col = "blue", lwd = 2, 
     xlab = "Time", ylab = "Survival Probability", 
     main = "Weibull regression", xlim = c(0, 17))
lines(ptime2, 1 - pct, col = "red", lwd = 2)
legend("topright", legend = c("Treatment 0", "Treatment 1"), 
       col = c("blue", "red"), lty = 1, bty = "n")

# Adress non-PH : piecewise Weibull model
# time-stratified Weibull model with time-by-treatment interaction
splidat1 <- survSplit(Surv(dthtime,dth=="Yes") ~ ., data = esrd, cut = c(1), episode = "interval")
piecewise_wei <- survreg(Surv(dthtime,event==1) ~ trt:strata(interval) + age + sex
                         + copd + dm + hypert +  heart + liver + neoplasia + vascular
                         + trt:heart + trt:hypert ,
                         dist = "weibull",data = splidat1)
summary(piecewise_wei)
print(ConvertWeibull(piecewise_wei)$HR)





