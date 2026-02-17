
## Intoduction ----

#import all libraries needed to the analysis
library(survival)


#import data 
data=read.csv('/Users/petros/Library/Mobile Documents/com~apple~CloudDocs/MSc Biostatistics/Survival Data Analysis/Kritsotakis Notes/HomeWork/hmwrk2025/esrd.csv',sep=";") 



#If dthtime is 0, we assume that the patient didn't die right away, but 0,5 day later

data$dthtime[data$dthtime==0]=0.5/365.25

## 1 ----
##Descriptive Statistics

#For Continues Variables

summary(data[data$treatmnt==1,c(2,13)])
summary(data[data$treatmnt==0,c(2,13)])

#Number of patient in each treatment group
nrow(data[data$treatmnt==1,])
nrow(data[data$treatmnt==0,])

#Number of patient died by the end of the study in each treatment group
nrow(data[data$dth==1 & data$treatmnt==1,])
nrow(data[data$dth==1 & data$treatmnt==0,])

#descriptives for age
sd(data$age[data$treatmnt==1])
sd(data$age[data$treatmnt==0])
sd(data$age)
mean(data$age[data$treatmnt==0])
mean(data$age[data$treatmnt==1])
mean(data$age)

#distribution of age
hist(data$age[data$treatmnt==1])
hist(data$age[data$treatmnt==0])

#distribution of dhtime
hist(data$dthtime[data$treatmnt==1])
hist(data$dthtime[data$treatmnt==0])

#assymetric distribution => median not mean
median(data$dthtime[data$treatmnt==0]) 
median(data$dthtime[data$treatmnt==1])
median(data$dthtime[data$treatmnt==0])*365.25
median(data$dthtime)
IQR(data$dthtime)
IQR(data$dthtime[data$treatmnt==0])
IQR(data$dthtime[data$treatmnt==1])

#For Categorical Variables

#Frequencies of each covariate in each treatment group
lapply(data[data$treatmnt==1,c(3:10)], table)
lapply(data[data$treatmnt==0,c(3:10)], table)

#Percentagies of each covariate in each treatment group 
lapply(data[data$treatmnt==1, c(3:10)], function(x) round(prop.table(table(x)) * 100,1))
lapply(data[data$treatmnt==0, c(3:10)], function(x) round(prop.table(table(x)) * 100,1))

#Frequencies and Percentagies of each covariate in overall study population
lapply(data[,c(3:10)], table)
lapply(data[, c(3:10)], function(x) round(prop.table(table(x)) * 100,1))


#Times of follow-up per treatment

#Median and IQR
tapply(data$dthtime, data$treatmnt, median)
tapply(data$dthtime, data$treatmnt, IQR)

#Min and Max 
tapply(data$dthtime, data$treatmnt, min)
tapply(data$dthtime, data$treatmnt, max)

## 2 ----

#Overall Death Risks per treatment group

#Calculating the overall population of each treatment group
N0=sum(data$treatmnt==0)
N1=sum(data$treatmnt==1)
N=N0+N1

#Calculating the events in each treatment group
events_0=sum(data$treatmnt==0 & data$dth==1)
events_1=sum(data$treatmnt==1 & data$dth==1)
total_deaths=sum(data$dth==1)

#Calculating the overall death risks on each treatment group
risk_0 = events_0 / N0
risk_1 = events_1 / N1

#Calculating CI of the point estimates

ci_0_lower = risk_0 - qnorm(0.975) * sqrt(risk_0 * (1 - risk_0) / N0)
ci_0_upper = risk_0 + qnorm(0.975) * sqrt(risk_0 * (1 - risk_0) / N0)

ci_1_lower = risk_1 - qnorm(0.975) * sqrt(risk_1 * (1 - risk_1) / N1)
ci_1_upper = risk_1 + qnorm(0.975) * sqrt(risk_1 * (1 - risk_1) / N1)

# Significance examination

chis=table(data$treatmnt,data$dth)
test=chisq.test(chis,correct = F)



#Estimation of Effectiveness of haemodialysis

#Relative Risk calculation

RR=risk_1/risk_0

se_logRR = sqrt((1-risk_1)/(risk_1*N1) + (1-risk_0)/(risk_0*N0))
logRR_lower = log(RR) - qnorm(0.975) * se_logRR
logRR_upper = log(RR) + qnorm(0.975) * se_logRR

RR_lower = exp(logRR_lower)
RR_upper = exp(logRR_upper)

# Effectiveness
eff=1-RR
eff_up=1-RR_upper
eff_low=1-RR_lower

## 3 ----

# Estimations of Mortality Rate

t0=sum(data$dthtime[data$treatmnt==0])
t1=sum(data$dthtime[data$treatmnt==1])

Rate0=events_0/t0
Rate1=events_1/t1

#Estimation of IRR, to check the significance

IRR=Rate1/Rate0
se_logIRR=sqrt((1/events_0)+(1/events_1))
logIRR_lower = log(IRR) - qnorm(0.975) * se_logIRR
logIRR_upper = log(IRR) + qnorm(0.975) * se_logIRR
IRR_lower = exp(logIRR_lower)
IRR_upper = exp(logIRR_upper) 

#Calculating Effectiveness

eff1=1-IRR
eff1_up=1-IRR_lower
eff1_low=1-IRR_upper


## 4 ----


#Make treatment as a factor 
data$treatmnt=factor(data$treatmnt,levels=c(0,1),labels=c("conservative","haemodialysis"))

#Run KM for each treatment group
KM.model=survfit(Surv(dthtime,dth)~treatmnt,data=data,conf.type="log-log")
quantile( KM.model, probs = c(0.25,0.5,0.65,0.75,0.9,0.95),conf.int = F )

#1-KM.model$surv[length(KM.model$surv)]
s=summary(KM.model)

#Plot of Surv Difference between the 2 groups
plot(KM.model,mark.time = F, main = "Comparison of Treatments for
Renal Disease",ylab = "Survival Probabilities",
     xlab = "Time from Remission to Relapse (Years)",lty = 1:2,col =
       c("black","red"))
lines(KM.model, conf.int = TRUE, lty = 3, lwd = 0.5, col = c("grey40","pink"))

legend("topright",lty = c(1,2,3,3),col = c("black","red","grey40","pink"),
       legend = c("Conservative (N=174)","Haemodialysis (N=8622)","Conservative CI","Haemodialysis CI"),bty = "n")


#Log-Rank Test - Wilcoxon Test
survdiff(Surv(dthtime,dth)~treatmnt,data=data,rho = 1)
survdiff(Surv(dthtime,dth)~treatmnt,data=data)
model1=survdiff(Surv(dthtime,dth)~treatmnt,data=data)


## 5 ----


## a ##

#Estimating Survival Probabilities for 0.5,1,3,5 years
summary(KM.model,times = c(0.5,1,3,5))


## b ##

#make a new data set for people that survived at least 6 months
data1=subset(data, dthtime > 0.5)

#Fit the Kaplan-Meier on this new dataset
km.cond=survfit(Surv(dthtime, dth) ~ treatmnt, data = data1, conf.type="log-log")

#Get the estimates for 1, 1.5, 2, 2.5 years after the first 6 months
summary(km.cond, times = c(1, 1.5, 2,2.5))


## 6 ----

#Applying an univariate cox PH model with age
model1=coxph(Surv(dthtime,dth)~age,data=data)
summary(model1)


## 7 ----

#Applying an univariate cox PH model with treatment
model2=coxph(Surv(dthtime,dth)~treatmnt,data=data)
summary(model2)
Efc=1-exp(model2$coefficients)


## 8 ----

#1st graph way to check proportional hazard assumption
plot(KM.model,mark.time = F,fun="cloglog", main = "Evaluation of PH Assumption",ylab = "log(-log(Survival Probabilities))",
     xlab = "Time from Remission to Relapse (Years) on log scale",lty = 1:2,col =
       c("blue","red")
     )
legend("topleft",lty = 1:2,col = c("blue","red"),
       legend = c("Conservative (N=174)","Haemodialysis (N=8622)"),bty = "n")


#2nd graph way to check proportional hazard assumption 

SchMar=cox.zph(model2,transform = "identity")
plot(SchMar[1],main = "scaled Schoenfeld residuals for Treatment status(Identity)",col="lightblue",cex=1)
abline(h = 0,col = "red", lty = 1)
legend("center",bty="n",lty = c(NA,1,1),col=c("black","lightblue","red"),pt.cex=1,pch=c(1,NA,NA),legend=c("Schoenfeld residual","Smoothed Curve","y=0"))

SchMar1=cox.zph(model2)
plot(SchMar1[1],main = "scaled Schoenfeld residuals for Treatment status(Default)",col="lightblue",cex=1)
abline(h = 0,col = "red", lty = 1)
abline(v=1)
legend("center",bty="n",lty = c(NA,1,1),col=c("black","lightblue","red"),pt.cex=1,pch=c(1,NA,NA),legend=c("Schoenfeld residual","Smoothed Curve","y=0"))

#Grambsch & Therneau test
SchMar
SchMar1

## 9 ----

#Re-red the data to avoid factors
data=read.csv("esrd.csv",sep=";") 
data$dthtime[data$dthtime==0]=0.5/365.25

#Run the cox model with time-dependent variable
tdcox= coxph(Surv(dthtime,dth)~treatmnt+tt(treatmnt),data=data,tt=function(x,t,...){x*t})
summary(tdcox)

#Efficacy plot

#Extract coef of the variables
b1=coef(tdcox)["treatmnt"]
b2=coef(tdcox)["tt(treatmnt)"]

#Setting the time
t=seq(0, 2, length.out = 200)

#Calculating HR(t) and VA(t)
HR_t=exp(b1 + b2 * t)
Eff2=1 - HR_t

#Calculating CI of Effectiveness
V=vcov(tdcox)[c("treatmnt", "tt(treatmnt)"),
                 c("treatmnt", "tt(treatmnt)")]
var_logHR=V[1,1] +  t^2 * V[2,2] +  2 * t * V[1,2]
se_logHR=sqrt(var_logHR)
HR_low=exp(b1 + b2 * t - 1.96 * se_logHR)
HR_up=exp(b1 + b2 * t + 1.96 * se_logHR)

Eff_low=1 - HR_up
Eff_up=1 - HR_low

#Plot 
plot(t, Eff2,
     type = "l",
     lwd = 2,
     col = "blue",
     ylim = c(-1, 1),
     xlab = "Time (years)",
     ylab = "Effectiveness: 1 − HR(t)",
     main = "Time-varying effectiveness of haemodialysis")
polygon(
  c(t, rev(t)),
  c(Eff_low, rev(Eff_up)),
  col = rgb(0, 0, 1, alpha = 0.2),
  border = NA
)
abline(h = 0, lty = 2)
legend("topright",
       legend = c("Effectiveness", "95% CI (delta method)"),
       lty = c(1, NA),
       lwd = c(2, NA),
       pch = c(NA, 15),
       pt.cex = c(NA, 2),
       col = c("blue", rgb(0, 0, 1, 0.2)),
       bty = "n")

## 10 ----



#Run non PH piecewise cox model

#Split Data on the cut points I chose
data2=survSplit(Surv(dthtime,dth)~treatmnt,data=data,cut=c(0.5,1))


x=rep(0,length(data2$dth))
data2$x1=x
data2$x2=x
data2$x1[data2$tstart>=0.5 & data2$dthtime<=1 & data2$treatmnt==1]=1
data2$x2[data2$tstart>=1 & data2$treatmnt==1]=1
pwcox= coxph(Surv(tstart,dthtime,dth)~treatmnt+x1+x2,data=data2)

summary(pwcox)

#Create CI

# coefficients
b=coef(pwcox)

# variance-covariance matrix
V=vcov(pwcox)

#Interval 2: treatmnt + x1
logHR_2=b["treatmnt"] + b["x1"]
var_2=V["treatmnt","treatmnt"]+V["x1","x1"]+2*V["treatmnt","x1"]

se_2=sqrt(var_2)

HR2=exp(logHR_2)
CI2=exp(c(logHR_2 - qnorm(0.975)*se_2,logHR_2 + qnorm(0.975)*se_2))

# Interval 3: treatmnt + x2
logHR_3=b["treatmnt"] + b["x2"]
var_3=V["treatmnt","treatmnt"] +V["x2","x2"] +2 * V["treatmnt","x2"]

se_3=sqrt(var_3)

HR3=exp(logHR_3)
CI3=exp(c(logHR_3 - qnorm(0.975)*se_3,logHR_3 + qnorm(0.975)*se_3))

#Point estimates and CI
HR2; CI2
HR3; CI3



## 11 ----

data_f=data

#Turn all variables into factors
data_f$treatmnt=factor(data_f$treatmnt,levels=c(0,1),labels=c("conservative","haemodialysis"))
data_f$sex=factor(data_f$sex,levels=c("F","M"),labels=c("Female","Male"))

vartofactor=c("copd","dm","hypert","heart","liver","neoplasia","vascular")

data_f[vartofactor] <- lapply(data_f[vartofactor], function(x) {
  factor(x, levels = c(0, 1), labels = c("No", "Yes"))
})

#Run the base multivariate coxph model
base_model=coxph(Surv(dthtime,dth)~treatmnt+sex+age+copd+dm+hypert+heart+liver+neoplasia+vascular,data=data_f)
summary(base_model)

#Covariate Patterns of "typical" patient on Conventional Treatment Group
nd0=data.frame(
  treatmnt   = factor("conservative", levels = levels(data_f$treatmnt)),
  sex        = factor("Male", levels = levels(data_f$sex)),
  age        = mean(data_f$age, na.rm = TRUE),
  copd       = factor("No", levels = levels(data_f$copd)),
  dm         = factor("No", levels = levels(data_f$dm)),
  hypert     = factor("No", levels = levels(data_f$hypert)),
  heart      = factor("No", levels = levels(data_f$heart)),
  liver      = factor("No", levels = levels(data_f$liver)),
  neoplasia  = factor("No", levels = levels(data_f$neoplasia)),
  vascular   = factor("No", levels = levels(data_f$vascular))
)
#Covariate Patterns of "typical" patient on Haemodialysis Treatment Group

nd1 <- nd0
nd1$treatmnt <- factor("haemodialysis", levels = levels(data_f$treatmnt))

#Fit Cox model
sf_cox <- survfit(base_model, newdata = rbind(nd0, nd1))


plot(KM.model,mark.time = F,main = "Predicted survival for multivariate COX model VS KM",
     xlab = "Length of Stay (years)",ylab = "Survival probability",
     col = c("lightblue","red"))
# Cox predicted survival
lines(sf_cox$time, sf_cox$surv[,1],
      lwd = 2, lty = 3, col = "blue")

lines(sf_cox$time, sf_cox$surv[,2],
      lwd = 2, lty = 3, col = "black")
legend("topright",bty = "n",lty = c(1,1,3,3),col = c("lightblue","red","blue","black"),
       legend = c("KM: Conservative","KM: Haemodialysis","COX: Conservative","COX: Haemodialysis"),ncol = 2)


#Run all possible models containing the base model and an interaction of treatment and a different baseline covariate
covariates=c("sex","age","copd","dm","hypert","heart","liver","neoplasia","vascular")
interaction_results=data.frame(Variable = character(), P_Value = numeric(), Significant = character(),stringsAsFactors = FALSE)

for (var in covariates) {
  #Update the model: Add 'treatmnt * var'
  #The update function is smart: .~. means "same outcome ~ same predictors"
  #We just add the interaction term.
  new_model=update(base_model, as.formula(paste(". ~ . + treatmnt *", var)))
  
  #Extract the P-value for the interaction term
  #We look at the summary coefficients. 
  #The interaction term usually appears at the bottom.
  #We use grep to find the row name containing both "treatmnt" and the variable name.
  coef_matrix=summary(new_model)$coefficients
  interaction_row=grep(paste0(":", var, "|", var, ":"), rownames(coef_matrix))
  
  #Safety check: if the interaction term exists
  if (length(interaction_row) > 0) {
    p_val=coef_matrix[interaction_row, "Pr(>|z|)"]
    
  #Store result
    interaction_results[nrow(interaction_results) + 1, ] <- list(
      Variable = var, 
      P_Value = round(p_val, 4),
      Significant = ifelse(p_val < 0.05, "YES", "No")
    )
  }
}

#View the final table
print(interaction_results)

#Cox-Snell Residuals

#Cox-Snell residual=delta-Martingale residual
data_f$csres = data_f$dth - residuals(base_model,type = "martingale")

#Fit K-M on Cox-Snell residual
fitcs = survfit(Surv(csres,dth)~1,data = data_f)

#Plot

#Order data throught CS-residual, so the line can be straight
data_f1=data_f
data_f1=data_f1[order(data_f1$csres),]

plot(fitcs,fun = "cloglog",conf.int = F,mark.time = F, 
     xlab = "Cox-Snell residuals(log scaled)",ylab = "Log(-log S(cs))",xaxt = "n",col="blue")
abline(h = 0,col = "black", lty = 2) #add y=0
abline(v = 1,col = "black", lty = 2) #add x=1
#Add the ticks
x_ticks=c(0.001, 0.01, 0.1, 1, 10)
axis(1, at = x_ticks, labels = c("0.001", "0.01", "0.1", "1", "10"))
#Add the line
lines(data_f1$csres,log(data_f1$csres),col="red",lty=3)
#Add legend
legend("topleft",bty = "n",lty=c(1,3),col = c("blue","red"),legend = c("Observed","Expected"))




#Re-examine PH assumption

#Schoenfeld Residual
SchMar2=cox.zph(base_model,transform = "identity")
plot(SchMar2[1],main = "Scaled Schoenfeld Residuals for Treatment(Identity)",col="lightblue",cex=0.25)
abline(h = 0,col = "red", lty = 1) #add y=0
legend("bottomright",bty="n",lty = c(NA,1,1),col=c("black","lightblue","red"),pt.cex=0.25,pch=c(1,NA,NA),legend=c("Schoenfeld residual","Smoothed Curve with CI","y=0"))


SchMar3=cox.zph(base_model)
plot(SchMar3[1],main = "Scaled Schoenfeld Residuals for Treatment(Default)",col="lightblue",cex=0.25)
abline(h = 0,col = "red", lty = 1) #add y=0
legend("center",bty="n",lty = c(NA,1,1),col=c("black","lightblue","red"),pt.cex=0.25,pch=c(1,NA,NA),legend=c("Schoenfeld residual","Smoothed Curve with CI","y=0"))


#Grambsch & Therneau test
SchMar2
SchMar3

#Address the violation of PH

#Re-read the data, to avoid the variables being functions


data_pw <- survSplit(Surv(dthtime,dth)~treatmnt+sex+age+copd+dm+hypert+heart+liver+neoplasia+vascular,data=data,cut = c(0.5, 1),  episode = "interval")
x=rep(0,length(data_pw$dth))
data_pw$x1=x
data_pw$x2=x
data_pw$x1[data_pw$interval==2& data_pw$treatmnt==1]=1
data_pw$x2[data_pw$interval==3 & data_pw$treatmnt==1]=1
pwcox= coxph(Surv(tstart,dthtime,dth)~treatmnt+x1+x2+sex+age+copd+dm+hypert+heart+liver+neoplasia+vascular,data=data_pw)

summary(pwcox)

#Create CI

# coefficients
b=coef(pwcox)

# variance-covariance matrix
V=vcov(pwcox)

#Interval 2: treatmnt + x1
logHR_2=b["treatmnt"] + b["x1"]
var_2=V["treatmnt","treatmnt"]+V["x1","x1"]+2*V["treatmnt","x1"]

se_2=sqrt(var_2)

HR2=exp(logHR_2)
CI2=exp(c(logHR_2 - qnorm(0.975)*se_2,logHR_2 + qnorm(0.975)*se_2))

# Interval 3: treatmnt + x2
logHR_3=b["treatmnt"] + b["x2"]
var_3=V["treatmnt","treatmnt"] +V["x2","x2"] +2 * V["treatmnt","x2"]

se_3=sqrt(var_3)

HR3=exp(logHR_3)
CI3=exp(c(logHR_3 - qnorm(0.975)*se_3,logHR_3 + qnorm(0.975)*se_3))

#Point estimates and CI
HR2; CI2
HR3; CI3


## 12a ----

KM.model1=survfit(Surv(dthtime,dth)~1,data=data)

km_df=data.frame(
  time   = KM.model1$time,
  surv   = KM.model1$surv
)

# drop invalid points (avoid log(0), log(-log(1)), log(-log(0)))
km_df=subset(km_df, time > 0 & surv > 0 & surv < 1)

km_df$x=log(km_df$time)
km_df$y=log(-log(km_df$surv))
lm_weib=lm(y~x,data = km_df)

gamma=lm_weib$coefficients[2] #Shape
l=exp(lm_weib$coefficients[1]) #Scale

#Fit Weibull model
weib_model=survreg(Surv(dthtime, dth)~1,data = data,dist ="weibull")

summary(weib_model)
#Because survreg run ATF Format we need to transform log(time ratios) to HR
shape = 1 / weib_model$scale #gamma
lambda = exp(-weib_model$coefficients/weib_model$scale) #scale weib PH

#Plot
plot(y ~ x, data = km_df,
     pch=16, cex=0.6, xlab="log(t)", ylab="log(-log(S(t)))",
     main="Weibull graphical fit check",col="darkblue")
abline(lm_weib, lwd=2,col="red",lty=1)
abline(log(lambda),shape,lwd=2,col="green",lty=2)
legend("topleft",bty="n",col=c("blue","red","green"),legend = c("Observed","Graphically estimated parameters","MLE parameters"),lty=c(NA,1,2),pch = c(16,NA,NA),pt.cex = c(0.6,NA,NA))


## 12b ----



weib_model1=survreg(Surv(dthtime,dth)~treatmnt+sex+age+copd+dm+hypert+heart+liver+neoplasia+vascular,data = data_f,dist ="weibull")
summary(weib_model1)

#BIC(pwcox)-BIC(weib_model1)

shape_m=1 / weib_model1$scale #gamma
lambda_m=exp(-weib_model1$coefficients/weib_model1$scale) 

#scale weib PH

plot(KM.model,mark.time = F,main = "Predicted survival for Weibull model VS KM VS Predicted Cox",
     xlab = "Length of Stay (years)",ylab = "Survival probability",
     col = c("lightblue","orange"))
# Add curves fitted by the Weibull model
pct = seq(0,1,by = 0.001)
lines(predict(weib_model1,newdata = nd0,type = "quantile",
              p = pct),1-pct,lty = 2,col = "blue")
lines(predict(weib_model1,newdata = nd1,type = "quantile",
              p = pct),1-pct,lty = 2,col = "red")
lines(sf_cox$time, sf_cox$surv[,1],
      lwd = 2, lty = 4, col = "darkblue")
lines(sf_cox$time, sf_cox$surv[,2],
      lwd = 2, lty = 4, col = "red4")
legend("topright",bty = "n",lty = c(1,1,2,2,4,4),col = c("lightblue","orange","blue","red","darkblue","red4"),
       legend = c("KM: Conservative","KM: Haemodialysis","WEI: Conservative","WEI: Haemodialysis","COX: Conservative","COX: Haemodialysis"),ncol = 3,cex=0.7)

#Check if the coef diff of the 2 models is significant differenct



# Cox
beta_cox=coef(base_model)
se_cox=sqrt(diag(vcov(base_model)))

# Weibull PH parameterization
beta_wei=-coef(weib_model1)[-1]/weib_model1$scale
se_wei=sqrt(diag(vcov(weib_model1)))[c(-1,-12)]/weib_model1$scale

#Calculate pvalue
z=(beta_cox-beta_wei)/sqrt(se_cox^2+se_wei^2)
p=2*(1-pnorm(abs(z)))

results=data.frame(Z = z,p_value = p)
#Results
round(results,2)





