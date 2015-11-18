library(foreign)
library(hdrcde)
data <- read.dta(file="arsenicrice2.dta")
Y <- read.dta(file="arsenicrice2.dta")

#We set weakly informative prior values
n <-  nrow(Y)
nu0 <- 1; eta0 <- 1; t20 <- 3;
mu0 <- mean(Y$arsenic); 
g20 <- s20 <- var(Y$arsenic)


# mu0 <- 100 #scenario 1
# s20 <- 100*s20; nu0 <- 100 #scenario 2
# t20 <- 100*t20; eta0 <- 100 #scenario 3
#Setup starting values 

m <- length(unique(Y$food_num)) #number of groups
n <- sv <- ybar <- rep(NA,m) #create empty vectors for group descriptions
for (i in 1:m) 
{
  n[i] <- sum(Y$food_num==i)
  sv[i] <- var(Y$arsenic[which(Y$food_num==i)])
  ybar[i] <- mean(Y$arsenic[which(Y$food_num==i)])
}
theta <- ybar; s2 <- mean(sv)
mu <- mean(theta); tau2 <- var(theta)

#Setup MCMC
set.seed(0808)
S <- 10000
THETA <- matrix(nrow=S, ncol=m)
RES <- matrix(nrow=S, ncol=3)
ALL <- matrix(nrow=S, ncol=3+m)

### Functions sampling from posteriors
newTheta <- function(n, ybar, s2, tau2, mu)
{
  v = 1/(n/s2 +1/tau2)
  e = v * (ybar*n/s2 + mu/tau2)
  new <- rnorm(1, e, sqrt(v))
  return(new)
}
newSigma2 <- function(m, n, nu0, s20, theta, Y)
{
  nun = nu0 + sum(n)
  ss <- nu0 * s20
  for(i in 1:m) ss = ss+sum((Y$arsenic[which(Y$food_num==i)] - theta[j])^2)
  sigma2 <- 1/rgamma(1, nun/2, ss/2)
  return(sigma2)
}
newMu <- function(m, theta, tau2, g20)
{
  v = 1/(m/tau2 + 1/g20)
  e = v *(m*mean(theta)/tau2 + mu0/g20)
  mu <- rnorm(1, e, v)
  return(mu)
}
newTau2 <- function(m, eta0, t20, theta, mu)
{
  etam = eta0 + m
  ss <- eta0*t20 + sum( (theta-mu) ^2 )
  tau2 <- 1/rgamma(1, etam/2, ss/2)
  return(tau2)
}

### Run MCMC
for(i in 1:S)
{
  #Get new values for parameters
  for(j in 1:m) theta[j] <- newTheta(n[j], ybar[j], s2, tau2, mu)
  s2 <- newSigma2(m, n, nu0, s20, theta, Y)
  mu <- newMu(m, theta, tau2, g20)
  tau2 <- newTau2(m, eta0, t20, theta, mu)
  
  #Store in chain
  THETA[i,] <- theta
  RES[i,] <- c(mu,s2,tau2)
  ALL[i,] <- c(theta,mu,s2,tau2)
}

graphdata <- data.frame("Iteration"=c(1:S), "Mu"=OTH[,1], "Sigma2"=OTH[,2], "Tau"=OTH[,3], 
                        "Theta_1" = THETA[,1], "Theta_2" = THETA[,2], "Theta_3" = THETA[,3],
                        "Theta_4" = THETA[,4], "Theta_5" = THETA[,])
# ggplot(graphdata,aes(x=Iteration,y=Sigma2)) +
#   theme_minimal(base_family = "") + geom_line(colour="wheat3")
# ggplot(graphdata,aes(x=Iteration,y=Tau2)) +
#   theme_minimal(base_family = "") + geom_line(colour="wheat3")
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

for (i in (1:8))
{
  # k <- c(median(ALL[,i]),  hdr(ALL[,i])$hdr[2,1:2])
  k <- quantile(ALL[,i], probs=c(0.025, 0.5, 0.975))
  print(round(k,3))
}