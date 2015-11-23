library(foreign)
library(hdrcde)
setwd("Homework_3_779")
data <- read.dta(file="arsenicrice2.dta")
Y <- read.dta(file="arsenicrice2.dta")

#We set weakly informative prior values
n <-  nrow(Y)
eta0 <- 1; t20 <- 3
mu0 <- mean(Y$arsenic) 
g20  <- var(Y$arsenic)
a = 1; b = 3; alpha = 1

# mu0 <- 100 #scenario 1
# a <- 100*s20; b <- 100*b #scenario 2
t20 <- 100*t20; eta0 <- 100 #scenario 3
#Setup starting values 

m <- length(unique(Y$food_num)) #number of groups
n <- s2 <- ybar <- rep(NA,m) 
for (i in 1:m) 
{
  n[i] <- sum(Y$food_num==i)
  s2[i] <- var(Y$arsenic[which(Y$food_num==i)])
  ybar[i] <- mean(Y$arsenic[which(Y$food_num==i)])
}
theta <- ybar; s2_0 <- mean(sv)
mu <- mean(theta); tau2 <- var(theta)
nu0 <- 1

#Setup MCMC
set.seed(0808)
S <- 10000
THETA <- S2 <- matrix(nrow=S, ncol=m)
RES <- matrix(nrow=S, ncol=4)
ALL <- matrix(nrow=S, ncol=4+2*m)
NUMAX <- 50

### Functions sampling from posteriors
newTheta <- function(n, ybar, s2, tau2, mu)
{
  v = 1/(n/s2 +1/tau2)
  e = v * (ybar*n/s2 + mu/tau2)
  new <- rnorm(length(n), e, sqrt(v))
  return(new)
}
newSigma2 <- function(m, n, nu0, s20, theta, y)
{
  nun = n + nu0
  ss <- nu0 * s20 + sum((y - rep(theta, times=n))^2) 
  sigma2 <- 1/rgamma(m, nun/2, ss/2)
  return(sigma2)
}
newMu <- function(m, theta, tau2, g20)
{
  v = 1/(m/tau2 + 1/g20)
  e = v *(m*mean(theta)/tau2 + mu0/g20)
  mu <- rnorm(1, e, sqrt(v))
  return(mu)
}
newTau2 <- function(m, eta0, t20, theta, mu)
{
  etam = eta0 + m
  ss = eta0*t20 + sum( (theta-mu) ^2 )
  tau2 <- 1/rgamma(1, etam/2, ss/2)
  return(tau2)
}
newS20 <- function(m, nu0, s2, a, b)
{
  L = a + 0.5*m*nu0
  R = b + 0.5*sum(1/s2^2)
  s20 <- rgamma(1, L, R)
}

### Run MCMC
for(i in 1:S)
{
  #Get new values for parameters
  theta <- newTheta(n, ybar, s2, tau2, mu)
  s2 <- newSigma2(m, n, nu0, s20, theta, Y$arsenic)
  mu <- newMu(m, theta, tau2, g20)
  tau2 <- newTau2(m, eta0, t20, theta, mu)
  s20 <- newS20(m, nu0, s2, a, b)
  
  x <- 1:NUMAX
    lpnu0<-m*(.5*x*log(s20*x/2)-lgamma(x/2))+
      (x/2-1)*sum(log(1/s2))+
      -x*(alpha+.5*s20*sum(1/s2))
  nu0<-sample(x,1,prob=exp(lpnu0-max(lpnu0)))
  
  
  #Store in chain
  THETA[i,] <- theta
  S2[i,] <- s2
  RES[i,] <- c(mu,tau2, nu0, s20)
  ALL[i,] <- c(theta,s2,mu,tau2, nu0, s20)
}

# graphdata <- data.frame("Iteration"=c(1:S), "Mu"=RES[,1], "Sigma2"=RES[,2], "Tau"=RES[,3], 
#                         "Theta_1" = THETA[,1], "Theta_2" = THETA[,2], "Theta_3" = THETA[,3],
#                         "Theta_4" = THETA[,4], "Theta_5" = THETA[,])

graphdata <- data.frame("Iteration"=c(1:S), "Mu"=RES[,1], "Tau2"=RES[,2], "Nu_0"=RES[,3], "Sigma2_0"=RES[,4], 
                        "Theta_1" = THETA[,1], "Theta_2" = THETA[,2], "Theta_3" = THETA[,3],
                        "Theta_4" = THETA[,4], "Theta_5" = THETA[,5], 
                        "Sigma_1" = S2[,1], "Sigma_2" = S2[,2], "Sigma_3" = S2[,3],
                        "Sigma_4" = S2[,4], "Sigma_5" = S2[,5])
# 
# ggplot(graphdata,aes(x=Iteration,y=Mu)) +
#   theme_minimal(base_family = "") + geom_line(colour="wheat3")
# ggplot(graphdata,aes(x=Iteration,y=Tau2)) +
#   theme_minimal(base_family = "") + geom_line(colour="wheat3")
# ggplot(graphdata,aes(x=Iteration,y=Nu_0)) +
#   theme_minimal(base_family = "") + geom_line(colour="wheat3")
# ggplot(graphdata,aes(x=Iteration,y=Sigma2_0)) +
#   theme_minimal(base_family = "") + geom_line(colour="wheat3")
# 
# ggplot(graphdata,aes(x=Iteration,y=Theta_1)) +
#   theme_minimal(base_family = "") + geom_line(colour="wheat3")
# ggplot(graphdata,aes(x=Iteration,y=Theta_2)) +
#   theme_minimal(base_family = "") + geom_line(colour="wheat3")
# ggplot(graphdata,aes(x=Iteration,y=Theta_3)) +
#   theme_minimal(base_family = "") + geom_line(colour="wheat3")
# ggplot(graphdata,aes(x=Iteration,y=Theta_4)) +
#   theme_minimal(base_family = "") + geom_line(colour="wheat3")
# ggplot(graphdata,aes(x=Iteration,y=Theta_5)) +
#   theme_minimal(base_family = "") + geom_line(colour="wheat3")
# 
# ggplot(graphdata,aes(x=Iteration,y=Sigma_1)) +
#   theme_minimal(base_family = "") + geom_line(colour="wheat3")
# ggplot(graphdata,aes(x=Iteration,y=Sigma_2)) +
#   theme_minimal(base_family = "") + geom_line(colour="wheat3")
# ggplot(graphdata,aes(x=Iteration,y=Sigma_3)) +
#   theme_minimal(base_family = "") + geom_line(colour="wheat3")
# ggplot(graphdata,aes(x=Iteration,y=Sigma_4)) +
#   theme_minimal(base_family = "") + geom_line(colour="wheat3")
# ggplot(graphdata,aes(x=Iteration,y=Sigma_5)) +
#   theme_minimal(base_family = "") + geom_line(colour="wheat3")

for (i in (1:14))
{
  print(round(unname(quantile(ALL[,i], probs=c(0.025, 0.5, 0.975))),3))
}
qmat=apply(THETA,2,quantile,probs=c(0.025, 0.5,0.975))
mu_ci = quantile(RES[,1], probs=(c(0.025, 0.5, 0.975)))
res <- data.frame("Rice"=c("Non-Basmati", "Basmati", "Beverages", "Cakes", "Cereal"),
                  "l95"=qmat[1,], "median"=qmat[2,], "u95"=qmat[3,], "mean"=ybar)
g <- ggplot(res, aes(x = Rice, group=Rice, colour=Rice)) + 
  labs(x="Rice Products", y="Arsenic concentration, mcg/serving") + 
  theme(legend.position="none", panel.background =element_rect(colour = "black")) +
  scale_y_continuous(breaks=seq(0, 7.5, 1)) +
  geom_hline(aes(yintercept=c(mu_ci[2])), size=0.7) +
  geom_hline(aes(yintercept=c(mu_ci[1])), linetype="dashed") +
  geom_hline(aes(yintercept=c(mu_ci[3])), linetype="dashed") +
  geom_errorbar(aes(ymin=l95, ymax=u95), width=.3, size=0.8) +
  geom_point(aes(y=median), fill="white", shape=21, size=5)  +
  geom_point(aes(y=mean), fill="red", shape=21, size=3)
g


