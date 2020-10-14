#install packages
install.packages("minpack.lm")
install.packages("boot")
library(minpack.lm)
library(boot)

# Exsampledata of mortalitybase, absolute frequency of deaths in age x
dx <- c(355,27,16,13,12,9,10,8,8,8,8,9,8,10,12,15,26,29,38,41,45,
        44,45,46,48,48,47,52,51,57,57,68,69,70,79,87,87,92,103,107,
        123,130,146,160,175,199,212,241,269,303,338,370,422,473,527,
        578,638,699,777,844,908,983,1073,1143,1232,1307,1369,1479,
        1513,1638,1742,1854,1939,2044,2170,2315,2473,2608,2803,3022,
        3210,3432,3600,3778,3887,3966,3951,3853,3706,3478,3220,2978,
        2631,2220,1883,1430,1116,830,588,399,260,381)
        
# Data of mortalitybase, absolute frequency of survivers in age x
lx <- rev(cumsum(rev(dx)))

#qx = obserevd probability of death between the ages of x and age x+1
qx <- dx/lx
logqx  <- log(dx/lx)

#w = limiting age of the mortality table/ database
w <- length(dx)-1
x<- 0:90
data <- data.frame(Y=rep(0:w, dx))

#Models:
in this example 4 models are used: 

#Gompertzfunction S(x)
GompertzFuncS <- function(x, m, s) { 
  return(exp(exp(-m/s)-exp((x-m)/s)))
}
#Inverse Gompertzfunction S(X)
InvGompertzFuncS <- function(x, m, s) { 
  return((1-exp(-exp(-(x-m)/s)))/(1-exp(-exp(m/s))))
}
#Weibullfunction S(X)
WeibullFuncS <- function(x, m, s) {
  return(exp(-(x/m)^(m/s)))
}
#Inverse Weibullfunction S(X)
InvWeibullFuncS <- function(x, m, s) {
  return(1-exp(-(x/m)^(-m/s)))
}

#Aggregation of models:
#Function to combine all models in one S(X)
S <- function(x, psi1, psi2, psi3, psi4, m1, m2, m3, m4, s1, s2, s3, s4) {
  return(psi1 * GompertzFuncS(x, m1, s1) +
           psi2 * InvGompertzFuncS(x, m2, s2) +
           psi3 * WeibullFuncS(x, m3, s3) +
           psi4 * InvWeibullFuncS(x, m4, s4))
}
myS <- function(x, psi1, psi2, psi3, psi4, m1, m2, m3, m4, s1, s2, s3, s4) {
  return(psi1 * GompertzFuncS(x, m1, s1) +
           psi2 * InvGompertzFuncS(x, m2, s2) +
           psi3 * WeibullFuncS(x, m3, s3) +
           psi4 * InvWeibullFuncS(x, m4, s4))
}
#Function to define probability of death  ages of x and age x+1 X as Q(X) = 1- S(X+1)/S(x)
Q <- function(x, psi1, psi2, psi3, m1, m2, m3, m4, s1, s2, s3, s4) {
  psi4 <- 1 - psi1 - psi2 - psi3
  return(1 - (S(x+1, psi1, psi2, psi3, psi4, m1, m2, m3, m4, s1, s2, s3, s4)/
                S(x, psi1, psi2, psi3, psi4, m1, m2, m3, m4, s1, s2, s3, s4)))
}
#Function to return the logarithmic probability of death  ages of x and age x+1 as log(Q(X))
logQ <- function(x, psi1, psi2, psi3, m1, m2, m3, m4, s1, s2, s3, s4) {
  psi4 <- 1 - psi1 - psi2 - psi3
  return(log(1 - (myS(x+1, psi1, psi2, psi3, psi4, m1, m2, m3, m4, s1, s2, s3, s4)/
                    myS(x, psi1, psi2, psi3, psi4, m1, m2, m3, m4, s1, s2, s3, s4))))
}
#fitting of parameters on the statistical values qx
qx <- qx[0:91]
PaperFit <- nlsLM(qx  ~ Q(x, psi1, psi2, psi3, m1, m2, m3, m4, s1, s2, s3, s4),
                  start = list(psi1=0.97, psi2=0.01, psi3=0.01, 
                               m1=85, m2=42, m3=.26, m4=23, 
                               s1=11, s2=14, s3=1, s4=8),
                  weights = 1/abs(qx) , control = nls.control(maxiter=1000))
#fitting of parameters on the logarithmic statistical values qx using logQ
logqx <- logqx[0:91]
logPaperFit <- nlsLM(logqx ~ logQ(x, psi1, psi2, psi3, m1, m2, m3, m4, s1, s2, s3, s4),
                     start = list(psi1=0.97, psi2=0.01, psi3=0.01, 
                                  m1=85, m2=42, m3=.26, m4=23, 
                                  s1=11, s2=14, s3=1, s4=8), weights = 1/abs(logqx), control = nls.control(maxiter=1000))
qxPaperFit <- predict(PaperFit, data=data.frame(x=0:90))
logqxPaperFit <- predict(logPaperFit, data=data.frame(x=0:90))
##Visualization:
plot(x, logSterbewahrscheinlichkeit)
lines(x, logqxPaperFit, col="red")
#fitting of parameters on the logarithmic statistical values non parametric Bootstrap
bootfit <- function(data, indices) {
  d <- data[indices, ]
  # Bestimme für d die logqx Werte
  dx <- as.numeric(table(d)) 
  #lx = observed number of people alive, relative to an original cohort, at age x
  lx <- rev(cumsum(rev(dx)))
  #qx = obserevd probability of death between the ages of x and age x+1
  qx <- dx/lx
  #w = imiting age of the mortality table/ database
  w <- length(dx)-1
  #x= age 
  x <- 0:90
  qx <- qx[0:91]
  PaperFit <- nlsLM(Sterbewahrscheinlichkeit ~ Q(x, psi1, psi2, psi3, m1, m2, m3, m4, s1, s2, s3, s4),
                       start = list(psi1=0.97, psi2=0.01, psi3=0.01, 
                                    m1=85, m2=42, m3=.26, m4=23, 
                                    s1=11, s2=14, s3=1, s4=8),
                       weights = 1/abs(qx), control = nls.control(maxiter=1000))
  qxPaperFit <- predict(PaperFit, data=data.frame(x=0:90))
  return(PaperFit$m$getPars())
}
logbootfit <- function(data, indices) {
  d <- data[indices, ]
  # Bestimme für d die logqx Werte
  dx <- as.numeric(table(d))
  #lx = observed number of people alive, relative to an original cohort, at age x
  lx <- rev(cumsum(rev(dx)))
  #qx = obserevd probability of death between the ages of x and age x+1
  qx <- dx/lx
  logqx  <- log(dx/lx)
  #w = imiting age of the mortality table/ database
  w <- length(dx)-1
  #x= age 
  x <- 0:90
  logqx<- logqx[0:91]
  logPaperFit <- nlsLM(logqx ~ logQ(x, psi1, psi2, psi3, m1, m2, m3, m4, s1, s2, s3, s4),
                       start = list(psi1=0.97, psi2=0.01, psi3=0.01, 
                                    m1=85, m2=42, m3=.26, m4=23, 
                                    s1=11, s2=14, s3=1, s4=8),
                       weights = 1/abs(logqx), control = nls.control(maxiter=1000))
  logqxPaperFit <- predict(logPaperFit, data=data.frame(x=0:90))
  return(logPaperFit$m$getPars())
}
#call non parametric bootstrap
set.seed(2)
np_boot <- boot(data=data, statistic=bootfit, R=1000)
set.seed(1)
lognp_boot <- boot(data=data, statistic=logbootfit, R=1)
set.seed(2)

#plot results
lines(x, logqxPaperFit, col="red")
lines(x, logQ(x,9.862212e-01, 2.189063e-03, 7.514800e-03, 8.311480e+01, 3.440280e+01, 2.878187e+01, 2.196682e+01, 9.342605e+00, 6.470701e+00, 2.004074e+02, 5.208513e+00) 
, col="blue")
lines(x, logQ(x, 0.917213634, 0.067941919, 0.004328894, 86.424663344,  57.935706039, 0.379915365, 34.177258127, 7.591399243, 4.847847541, 0.928726117, 4.373326151  ), col="red")
