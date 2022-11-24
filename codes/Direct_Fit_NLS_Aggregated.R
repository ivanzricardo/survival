#############################################################################################
#Fitting Survival Models From Aggregate Data
#Authors: Solange Borges (solangeborges@gmail.com), Ivan Zimmermann (ivanzricardo@gmail.com)
#Date: 04/13/2021
#############################################################################################

#Instal and call packages  
install.packages("flexsurv")
library(flexsurv)

#Load survival data
data <- read.csv2('https://raw.githubusercontent.com/ivanzricardo/survival/main/data/OS_anastrozole.csv')

#Save time and survival vectors
time <- data$time
survival <- data$survival

#Fit desired survival functions with NLS
reg_exp <- nls(survival ~ exp(-l_e*time),data, start=list(l_e=0.1)) #Exponential 
reg_wei <- nls(survival ~ exp(-l_w*(time^g_w)),data, start=list(l_w=1,g_w=0.01)) #Weibull
reg_logl <- nls(survival ~ 1/(1 + l_l*time^g_l),data, start=list(l_l=0.1,g_l=0.01)) #Log-Logistic 
reg_gomp <- nls(survival ~ exp(-((g_g)/l_g)*(exp(l_g*time) - 1)),data, start=list(l_g=0.1,g_g=0.01), algorithm="port", lower=c(-Inf,-Inf)) #Gompertz
reg_lnorm <- nls(survival ~ 1- pnorm(((log(time)-mu_ln)/sd_ln),mu_ln,sd_ln),data, start=list(mu_ln=1,sd_ln=1)) #Lognormal
reg_gam_gen <- nls(survival ~ 1 - pgengamma(time, mu_g, sigma, Q),data, start=list(sigma=1,mu_g=1, Q=1), algorithm="port", lower=c(-Inf,-Inf,-Inf)) #Generalized Gama 

#Store summary fitting
s_reg_exp <- summary(reg_exp)
s_reg_wei<- summary(reg_wei)
s_reg_logl <- summary(reg_logl)
s_reg_gomp <- summary(reg_gomp)
s_reg_lnorm <- summary(reg_lnorm)
s_reg_gam_gen <- summary(reg_gam_gen)

#Store curve parameters
p_reg_exp <- data.frame(s_reg_exp[["parameters"]])
  p_reg_exp$Model <- c("Exponential")
p_reg_wei <- data.frame(s_reg_wei[["parameters"]])
  p_reg_wei$Model <- c("Weibull")
p_reg_logl <- data.frame(s_reg_logl[["parameters"]])
  p_reg_logl$Model <- c("Log-logistic")
p_reg_gomp <- data.frame(s_reg_gomp[["parameters"]])
  p_reg_gomp$Model <- c("Gompertz")
p_reg_lnorm <- data.frame(s_reg_lnorm[["parameters"]])
  p_reg_lnorm$Model <- c("Lognormal")
p_reg_gam_gen <- data.frame(s_reg_gam_gen[["parameters"]])
  p_reg_gam_gen$Model <- c("Generalized Gamma")

Parameters <- rbind(p_reg_exp,p_reg_wei,p_reg_logl,p_reg_gomp,p_reg_lnorm, p_reg_gam_gen)

#Add fitted data to observed survival 
horizon <- 10 #time beyond clinical trial data
data$time_pred <- data$time
data[(nrow(data)+1):(nrow(data)+horizon),"time_pred"] <- seq(from = round(as.numeric(data[nrow(data),1]+1), digits = 0), 
                                                             to = round(as.numeric(data[nrow(data),1]+horizon), digits = 0), by = 1)
time_pred <-data$time_pred
data$ajuste_exp <- exp(-p_reg_exp[1,1]*time_pred)
data$ajuste_wei <- exp(-p_reg_wei[1,1]*(time_pred^p_reg_wei[2,1]))
data$ajuste_logl <- 1/(1 + p_reg_logl[1,1]*time_pred^p_reg_logl[2,1])
data$ajuste_gomp <- exp(-((p_reg_gomp[2,1])/p_reg_gomp[1,1])*(exp(p_reg_gomp[1,1]*time_pred) - 1))
data$ajuste_lnorm <- 1- pnorm(((log(time_pred)-p_reg_lnorm[1,1])/p_reg_lnorm[2,1]),p_reg_lnorm[1,1],p_reg_lnorm[2,1])
data$ajuste_gam_gen <- 1 - pgengamma(time_pred, p_reg_gam_gen[2,1], p_reg_gam_gen[1,1], p_reg_gam_gen[3,1])

#Plot Kaplan-Meier data and fitted curves
plot(data$time_pred,data$survival, main="Visual Inspection", xlab = "Time", ylab = "Survival", type = "l", col="black",lty=1,lwd=2)
line_exp <- lines(data$time_pred,data$ajuste_exp,col="red",lty=1,lwd=1)
line_wei <- lines(data$time_pred,data$ajuste_wei,col="green3",lty=1,lwd=1)
line_logl <- lines(data$time_pred, data$ajuste_logl,col="blue",lty=1,lwd=1)
line_gomp <- lines(data$time_pred,data$ajuste_gomp,col="purple",lty=1,lwd=1)
line_lnorm <- lines(data$time_pred,data$ajuste_lnorm,col="yellow",lty=1,lwd=1)
line_gam_gen <- lines(data$time_pred,data$ajuste_gam_gen,col="gray",lty=1,lwd=1)
legend("topright",
       c("Kaplan-Meier","Exponential", "Weibull", "Log-logistic", "Gompertz", "Lognormal", "Generalized Gamma"),
       fill=c("black","red", "green3", "blue", "purple","yellow","gray"),
       border = F,
       cex=0.8, pt.cex = 1,
       bty = "n"
)

#Goodness-of-fit
#The lower the value for AIC, the better the fit of the model. 
#The absolute value of the AIC value is not important. 
#It can be positive or negative.
AIC(reg_exp, reg_wei, reg_logl, reg_gomp, reg_lnorm, reg_gam_gen)
BIC(reg_exp, reg_wei, reg_logl, reg_gomp, reg_lnorm, reg_gam_gen)
AIC_BIC <- cbind(AIC(reg_exp, reg_wei, reg_logl, reg_gomp, reg_lnorm, reg_gam_gen),BIC(reg_exp, reg_wei, reg_logl, reg_gomp, reg_lnorm, reg_gam_gen))

#Export results
write.csv2(Parameters,"Parameters.csv")
write.csv2(AIC_BIC,"AIC_BIC.csv")

######################################## THIS IS THE END ##############################################

