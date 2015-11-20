library(foreign)
library(hdrcde)
data <- read.dta(file="arsenicrice2.dta")
Y <- read.dta(file="arsenicrice2.dta")

#We set weakly informative prior values
n <-  nrow(Y)
nu0 <- 1; omega0 <- 1; 
s2_eta0 <- 10; s2_y0 <- 10
mu0 <- mean(Y$arsenic); 
g20 <- var(Y$arsenic)

# mu0 <- 100 #scenario 1
# s20 <- 100*s20; nu0 <- 100 #scenario 2
# t20 <- 100*t20; eta0 <- 100 #scenario 3

#Setup starting values 
m <- length(unique(Y$food_num)) #number of groups
n <- sv_y <- ybar <- rep(NA,m) #create empty vectors for group descriptions
for (i in 1:m) 
{
  n[i] <- sum(Y$food_num==i)
  sv_y[i] <- var(Y$arsenic[which(Y$food_num==i)])
  ybar[i] <-  mean(Y$arsenic[which(Y$food_num==i)])
}
eta <- ybar - mean(Y$arsenic); s2_y <- mean(sv_y)
mu <- mean(Y$arsenic); s2_eta <- var(eta)
xi <- 1

#Setup MCMC
set.seed(1)
S <- 1000
ETA <- matrix(nrow=S, ncol=m)
RES <- matrix(nrow=S, ncol=4)
ALL <- matrix(nrow=S, ncol=4+m)

sigsq.y = s2_y
sigsq.n = s2_eta
y.bar = ybar
gam0 = g20
y = Y$arsenic
delta0 = omega0
phi.sq0 = s2_eta0
psi.sq0 = s2_y0
N = sum(n)
ALL2 = ALL

#FUNCTIONS
newS2eta <- function(m, s2_eta0, eta, omega0)
{
  ss <- s2_eta0*omega0 + sum( eta^2 )
  s2_eta <- 1/rgamma(1, (omega0 + m)/2, ss/2)
  # s2_eta <- 1/rgamma(1, shape=((omega0 + m)/2), rate=(ss/2))
  return(s2_eta)
}
newS2y <- function(nu0, n, s2_y0, y, m, eta, xi, mu)
{
  nun = nu0 + sum(n)
  ss = sum((y - mu - xi*rep(eta, times=n))^2) + nu0*s2_y0
  s2_y <- 1/rgamma(1, nun/2, ss/2)
  # s2_y <- 1/rgamma(1, shape=(nun/2), rate=(ss/2))
  return(s2_y)
}
newMu <- function(y, xi, eta, s2_y, g20, mu0, n)
{
  v = ( sum(n)/s2_y + 1/g20 )
  ss = sum(y-rep(eta,times=n)*xi)/s2_y + mu0/g20
  e = v*(ss/s2_y + mu0/g20)
  mu = rnorm(1, ss/v, sqrt(1/v))
  return(mu)
}
newEta <- function(xi, mu, ybar, n, s2_y, s2_eta)
{
  v = 1/(n*xi^2/s2_y + 1/s2_eta)
  # v = 1/(n/xi^2/s2_y + 1/s2_eta) 
  e = v*xi*(n*ybar - mu*n)/s2_y
  # e = v*(n*ybar - mu*n)/xi/s2_y
  eta = rnorm(m, e, sqrt(v))
  return(eta)
}
newXi <- function(eta, m, y, n, mu, s2_y)
{
  v = (sum(n*eta^2)/s2_y + 1)
  # v = (sum(n/eta^2)/s2_y + 1)
  e = sum((y-mu)*rep(eta,times=n))/s2_y
  # e = sum((y-mu)/rep(eta,times=n))/s2_y
  xi = rnorm(1, e/v, sqrt(1/v))
  return(xi)
}
#RUN MCMC
for(i in 1:S)
{
  #Get new values for parameters
  eta <- newEta(xi, mu, ybar, n, s2_y, s2_eta)
  mu <- newMu(Y$arsenic, xi, eta, s2_y, g20, mu0, n)
  xi <- newXi(eta, m, Y$arsenic, n, mu, s2_y)
  s2_eta <- newS2eta(m, s2_eta0, eta, omega0)
  s2_y <- newS2y(nu0, n, s2_y0, Y$arsenic, m, eta, xi, mu)
  #Store in chain
  ETA[i,] <- eta
  RES[i,] <- c(mu,s2_y, s2_eta, xi)
  ALL[i,] <- c(eta, mu, xi, s2_eta, s2_y)
}
#CHECK
m <- length(unique(Y$food_num)) #number of groups
n <- sv_y <- ybar <- rep(NA,m) #create empty vectors for group descriptions
for (i in 1:m) 
{
  n[i] <- sum(Y$food_num==i)
  sv_y[i] <- var(Y$arsenic[which(Y$food_num==i)])
  ybar[i] <-  mean(Y$arsenic[which(Y$food_num==i)])
}
eta <- ybar - mean(Y$arsenic); s2_y <- mean(sv_y)
mu <- mean(Y$arsenic); s2_eta <- var(eta)
xi <- 0

#Setup MCMC
set.seed(0808)
S <- 1000
ETA <- matrix(nrow=S, ncol=m)
RES <- matrix(nrow=S, ncol=4)

sigsq.y = s2_y
sigsq.n = s2_eta
y.bar = ybar
gam0 = g20
y = Y$arsenic
delta0 = omega0
phi.sq0 = s2_eta0
psi.sq0 = s2_y0
N = sum(n)
ALL2 = matrix(nrow=S, ncol=4+m)
for(i in 1:S)
{
  A <- n*xi^2/sigsq.y+1/sigsq.n
  B <- xi*(n*y.bar-n*mu)/sigsq.y
  eta <- rnorm(m,B/A,sqrt(1/A))
  A <- N/sigsq.y + 1/gam0
  B <- sum((y-xi*rep(eta,times=n)))/sigsq.y + mu0/gam0
  mu <- rnorm(1,B/A,sqrt(1/A))
  A <- sum(n*eta^2)/sigsq.y + 1
  B <- sum((y-mu)*rep(eta,times=n))/sigsq.y
  xi <- rnorm(1,B/A,sqrt(1/A))
  A <- (m+delta0)/2
  B <- .5*(sum(eta^2) + delta0*phi.sq0)
  sigsq.n <- 1/rgamma(1,shape=A,rate=B)
  A <- (N+nu0)/2
  sum.sq <- sum((y-mu-xi*rep(eta,times=n))^2)
  B <- .5*(sum.sq + nu0*psi.sq0)
  sigsq.y <- 1/rgamma(1,shape=A,rate=B)
  ALL2[i,] <- c(eta, mu, xi, sigsq.n, sigsq.y)
}
graphdata <- data.frame(
  "Iteration"=c(1:S), "Mu"=RES[,1], "s2_y"=RES[,2], "s2_eta"=RES[,3], "xi"=RES[,4],
  "Eta_1" = ETA[,1], "Eta_2" = ETA[,2], "Eta_3" = ETA[,3], 
  "Eta_4" = ETA[,4], "Eta_5" = ETA[,5])

# ggplot(graphdata,aes(x=Iteration,y=Eta_5)) +
#   theme_minimal(base_family = "") + geom_line(colour="wheat3")
