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
xi <- 0

#Setup MCMC
set.seed(0808)
S <- 10000
ETA <- matrix(nrow=S, ncol=m)
RES <- matrix(nrow=S, ncol=4)
ALL <- matrix(nrow=S, ncol=4+m)

#FUNCTIONS
newS2eta <- function(m, s2_eta, s2_eta0, eta, omega0)
{
  nomega = omega0 + m
  ss <- s2_eta0*omega0 + sum( eta^2 )
  s2_eta <- 1/rgamma(1, nomega/2, ss/2)
  return(s2_eta)
}
newS2y <- function(nu0, n, s2_y0, Y, m, eta, xi, mu)
{
  nun = nu0 + sum(n)
  ss = nu0 * s2_y0
  for(i in 1:m) ss = ss+sum((Y$arsenic[which(Y$food_num==i)] - mu - xi*eta[j])^2)
  s2_y <- 1/rgamma(1, nun/2, ss/2)
  return(s2_y)
}
newMu <- function(Y, xi, eta, s2_y, g20, mu0, n)
{
  v = 1/(sum(n)/s2_y + 1/g20)
  sumd = 0
  for(i in 1:m) sumd = sumd + sum(Y$arsenic[which(Y$food_num==i)] - xi*eta[i])
  e = v* sumd/s2_y + mu0/g20
  mu = rnorm(1, e, v)
  return(mu)
}
newEta <- function(xi, mu, Y, n, s2_y, s2_eta)
{
  v = 1/(n*xi^2/s2_y + 1/s2_eta)
  e = v*xi/s2_y*sum(Y$arsenic[which(Y$food_num==i)] - mu)
  eta = rnorm(1, e, v)
  return(eta)
}
newXi <- function(eta, m, Y, n, mu, s2_y)
{
  v = 1/(sum(n*eta^2)/s2_y + 1)
  sumd = 0
  for(i in 1:m) sumd = sumd + sum((Y$arsenic[which(Y$food_num==i)] - mu)*eta[i])
  e = v * sumd/s2_y
  xi = rnorm(1, e, v)
  return(xi)
}
#RUN MCMC
for(i in 1:S)
{
  #Get new values for parameters
  for(j in 1:m) eta[j] <- newEta(xi, mu, Y, n[j], s2_y, s2_eta)
  s2_y <- newS2y(nu0, n, s2_y0, Y, m, eta, xi, mu)
  s2_eta <- newS2eta(m, s2_eta, s2_eta0, eta, omega0)
  mu <- newMu(Y, xi, eta, s2_y, g20, mu0, n)
  xi <- newXi(eta, m, Y, n, mu, s2_y)
  
  #Store in chain
  ETA[i,] <- eta
  RES[i,] <- c(mu,s2_y, s2_eta, xi)
  ALL[i,] <- c(eta, mu,s2_y, s2_eta, xi)
}
graphdata <- data.frame(
  "Iteration"=c(1:S), "Mu"=RES[,1], "s2_y"=RES[,2], "s2_eta"=RES[,3], "xi"=RES[,4],
  "Eta_1" = ETA[,1], "Eta_2" = ETA[,2], "Eta_3" = ETA[,3], 
  "Eta_4" = ETA[,4], "Eta_5" = ETA[,5])

ggplot(graphdata,aes(x=Iteration,y=Eta_5)) +
  theme_minimal(base_family = "") + geom_line(colour="wheat3")
