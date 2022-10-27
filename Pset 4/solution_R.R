# packages ----------------------------------------------------------------
library(here)
library(data.table)
library(clipr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(dplyr)
library(tictoc)
library(expm) 
library(purrr)
library(dlm)
library(lubridate)
# Problem 1 ---------------------------------------------------------------

set.seed(6314)
phi <- 0.6
t <- 5
N <- 100
B <- 5000
# functions to generate the data
generate_y0 <- function(alphai,phi) {
  return(rnorm(n = 1,mean = alphai/(1-phi),sd =1/sqrt(1-phi^2)))
}     
data_generation <- function(phi,N,t) {
  alpha_i <- rnorm(N,0,1)
  epsilon <- matrix(rnorm(N*t,0,1),nrow = N,ncol = t)
  y0 <- sapply(alpha_i,generate_y0,phi=phi)
  data <- matrix(0,nrow = N,t+1)
  data[,1] <- y0
  for (i in 2:ncol(data)) {
    data[,i] <- phi*data[,i-1]+alpha_i+epsilon[,i-1]
  }
  return(data)
  }

# ml solver for phi
ml_solver <- function(data,A) {
  y <- c(t(data[,-1]))
  y_ <- c(t(data[,-ncol(data)]))
  phi_ml <- solve(t(y_)%*%A%*%y_)%*%t(y_)%*%A%*%y
  return(phi_ml)
}
# create A matrix outside loop
A_t <- diag(t)-1/t*matrix(1,t,t)
A <- kronecker(diag(N),A_t)
phis_ml <- rep(NA, B)
# monte carlo for ML estimation
for (i in 1:B) {
  data <- data_generation(phi,N,t) 
  phis_ml[i] <- ml_solver(data,A)
}
# bias for ML estimation
bias_ml <- mean(phis_ml)-phi
rmse_ml <- sqrt(sum((phis_ml-phi)^2)/B)

print(paste0('bias ML= ',bias_ml))
print(paste0('RMSE ML= ',rmse_ml))

# indirect inference part
theta_tilde_h <- function(phi_choice,h,N,t,A) {
  phi_hat <- c()
  for (i in 1:h) {
    data_h <- data_generation(phi_choice,N,t)
    phi_hat[i] <- ml_solver(data_h,A)
  }
  return(mean(phi_hat))
}

obj_function <- function(phi_choice,theta_data,h,N,t,A) {
  mean_phi_hat <- theta_tilde_h(phi_choice,h,N,t,A)
  return((theta_data-mean_phi_hat)^2)
}

h<- 10
phis_ii <- rep(NA, B)
for (i in 1:B) {
  data <- data_generation(phi,N,t)
  theta_data <- ml_solver(data,A)
  phis_ii[i] <- optimize(obj_function,h=h,theta_data=theta_data,N=N,t=t,A=A,interval = c(-1,1))$minimum
}

#bias II
bias_ii <- mean(phis_ii)-phi
rmse_ii <- sqrt(sum((phis_ii-phi)^2)/B)
print(paste0('Bias II = ',bias_ii))
print(paste0('RMSE II = ',rmse_ii))


save.image(file =here("problem1.RData"))


# histograms Monte Carlo simulation
plt1 <- ggplot(data.frame('Estimates'=phis_ii),aes(x=Estimates,y=..density..))+
  geom_histogram()+
  geom_vline(aes(xintercept=phi),
             color="blue", linetype="dashed", size=1)+
  ggtitle('Histogram II')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
  

plt2 <- ggplot(data.frame('Estimates'=phis_ml),aes(x=Estimates,y=..density..))+
  geom_histogram()+
  geom_vline(aes(xintercept=phi),
             color="blue", linetype="dashed", size=1)+
  ggtitle('Histogram ML')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))


plt3 <- ggarrange(plt1,plt2,ncol=2 )
plt3
ggsave(plot=plt3,here('histogram.png'),width = 14,height = 8,dpi=200)

  

# Problem 2 ---------------------------------------------------------------

# identity matrix function
I <- function(n) {
  return(diag(1,n))
}

# return matrices for the state space representation given the vector of parameters of the model
state_space_representation <- function(par) {
  # y = Z * beta_t 
  # beta_t = f* beta_t_1 + v_t 
  # v_t covariance matrix is Q
  gamma1 <- par[1]
  gamma2 <- par[2]
  gamma3 <- par[3]
  gamma4 <- par[4]
  
  psi11 <- par[5]
  psi12 <- par[6]
  psi13 <- par[7]
  psi14 <- par[8]
  
  psi21 <- par[9]
  psi22 <- par[10]
  psi23 <- par[11]
  psi24 <- par[12]
  
  phi1 <- par[13]
  phi2 <- par[14]
  
  var1 <- par[15]
  var2 <- par[16]
  var3 <- par[17]
  var4 <- par[18]
  Z <- rbind(c(gamma1,0,1,0,0,0,0,0,0,0),
             c(gamma2,0,0,0,1,0,0,0,0,0),
             c(gamma3,0,0,0,0,0,1,0,0,0),
             c(gamma4,0,0,0,0,0,0,0,1,0))
  
  T <- rbind(c(phi1,phi2,0,0,0,0,0,0,0,0),
             c(1,0,0,0,0,0,0,0,0,0),
             c(0,0,psi11,psi21,0,0,0,0,0,0),
             c(0,0,1,0,0,0,0,0,0,0),
             c(0,0,0,0,psi12,psi22,0,0,0,0),
             c(0,0,0,0,1,0,0,0,0,0),
             c(0,0,0,0,0,0,psi13,psi23,0,0),
             c(0,0,0,0,0,0,1,0,0,0),
             c(0,0,0,0,0,0,0,0,psi14,psi24),
             c(0,0,0,0,0,0,0,0,1,0)) 
  Q <- diag(c(1,0,var1,0,var2,0,var3,0,var4,0))
  return(list('Z'=Z,'T'=T,'Q'=Q))
}

# Kalman Filter: conventions
# a_t=T*a_{t-1}+eta_t. eta_t~N(0,Q)
# yt=Z*a_t+e_t e_t~N(0,H)

# Initial state: mi0 e sigma0
# Y: Matrix. Each column is an observation.
kalman_filter <- function(a0,P0,Z,T,Q,H,y) {
  nt <- ncol(y)
  at <- matrix(0,nrow = dim(T)[1],ncol=nt)
  atm <- at
  P <- array(0,dim = c(nt,length(a0),length(a0)))
  ptm <- P
  K <- array(0,dim = c(nt,dim(P)[2],dim(H)[2]))
  inov <- y
  cov_inov_matrix <- array(0,dim=c(nt,dim(y)[1],dim(y)[1]))
  for (i in 1:nt) {
    if (i==1) {
      atm[,i] <- T%*%a0
      ptm[i,,] <- T%*%P0%*%t(T)+Q
      ytm <- Z%*%atm[,i]
      F <- Z%*%ptm[i,,]%*%t(Z)+H
      K[i,,] <- ptm[i,,]%*%t(Z)%*%solve(F)
      at[,i] <- atm[,i]+K[i,,]%*%(y[,i]-Z%*%atm[,i])
      P[i,,] <- ptm[i,,]-K[i,,]%*%Z%*%ptm[i,,]
      inov[,i] <- y[,i]-ytm
      cov_inov_matrix[i,,] <- F
    }
    else{
      atm[,i] <- T%*%at[,i-1]
      ptm[i,,] <- T%*%P[i-1,,]%*%t(T)+Q
      ytm <- Z%*%atm[,i]
      F <- Z%*%ptm[i,,]%*%t(Z)+H
      K[i,,] <- ptm[i,,]%*%t(Z)%*%solve(F)
      at[,i] <- atm[,i]+K[i,,]%*%(y[,i]-Z%*%atm[,i])
      P[i,,] <- ptm[i,,]-K[i,,]%*%Z%*%ptm[i,,]
      inov[,i] <- y[,i]-ytm
      cov_inov_matrix[i,,] <- F
    }
  }
  return(list('gain'=K,'state'=at,'state_m1'=atm,'inov'=inov,'cov_inov_matrix'=cov_inov_matrix,'P'=P,'p_pred'=ptm))
}

recover_level_C <- function(deltac,delta_hat,C0) {
  C <- deltac
  for (i in 1:length(C)) {
    if (i==1) {
      C[i] <- C0+deltac[i]+delta_hat
    }
    else{
      C[i] <- C[i-1]+deltac[i]+delta_hat
    }
  }
  return(C)
}

log_likelihood <- function(a0,P0,Z,T,Q,H,y) {
  kalman <- kalman_filter(a0,P0,Z,T,Q,H,y)
  innov <- kalman$inov
  sigma_innov <- kalman$cov_inov_matrix
  loglike <- 0
  n <- ncol(innov)
  for (i in 1:n){
    loglike <- loglike-1/2*log(det(sigma_innov[i,,]))-1/2*as.numeric(t(innov[,i])%*%solve(sigma_innov[i,,])%*%innov[,i])
  }
  return(loglike)
}
obj_function_ML <- function(par,a0,y) {
  ss <- state_space_representation(par)
  Z <- ss$Z
  T <- ss$T
  Q <- ss$Q
  P0 <- T%*%Q%*%t(T)
  var_inov_y <- matrix(0,nrow=dim(y)[1],ncol=dim(y)[1])
  return(-log_likelihood(a0,P0,Z,T,Q,H = var_inov_y,y))
}
# load and clean data -----------------------------------------------------
data <- fread(here('Data','kim_data.csv'),drop = 6)
data$date <- as.character(data$date)
data$date <- paste0('01.',data$date)
data$date <- as.Date(format(as.Date(data$date,format = "%d.%m.%y"),"19%y-%m-%d"))
data <- data %>% arrange(ymd(data$date))
# apply log to data
data <- data %>% mutate_at(c(2:5),.funs=log)

# first difference
data <- data %>% mutate_at(c(2:5),.funs = function(x) x-lag(x))
dada_without_demean <- data
# parameter of the model

gamma1 <- 0.717
gamma2 <- 0.521
gamma3 <- 0.47
gamma4 <- 0.602

psi11 <- -0.04
psi12 <- -0.087
psi13 <- -0.414
psi14 <- 0.108

psi21 <- -0.137
psi22 <- 0.154
psi23 <- -0.206
psi24 <- 0.448

phi1 <- 0.545
phi2 <- 0.032
var1 <- (0.488*1e-2)^2
var2 <- (0.769*1e-2)^2
var3 <- (0.735*1e-2)^2
var4 <- (0.540*1e-2)^2

data_filter <- data[date<='1987-12-01' & date>='1959-02-01',]
#demean
data_filter[,2:5] <- sweep(data_filter[,2:5] ,2,colMeans(data_filter[,2:5],na.rm = T))

dada_without_demean_filter <- dada_without_demean[date<='1987-12-01' & date>='1959-02-01',]
y <- t(as.matrix(data_filter[,2:5]))
Y <- t(as.matrix(dada_without_demean_filter[,2:5]))

write.table(y,file=here('Mlab','y.csv'),row.names = F)
write.table(Y,file=here('Mlab','Y_with_mean.csv'),row.names = F)

par <- c(gamma1,gamma2,gamma3,gamma4,
         psi11,psi12,psi13,psi14,
         psi21,psi22,psi23,psi24,
         phi1,phi2,
         var1,
         var2,
         var3,
         var4)
ss <- state_space_representation(par)
Q <- ss$Q
T <- ss$T
Z <- ss$Z

# Kalman Filter: conventions
# a_t=T*a_{t-1}+R_T*eta_t. eta_t~N(0,Q)
# yt=Z*a_t+e_t e_t~N(0,H)
# Initial state: mi0 e sigma0

a0 <- matrix(0,nrow=nrow(T),ncol=1)
var_inov_y <- matrix(0,nrow=dim(y)[1],ncol=dim(y)[1])

P0 <- T%*%Q%*%t(T)


# solve part a ------------------------------------------------------------
results <- kalman_filter(a0,P0,Z,T,Q,H = var_inov_y,y = y)
state <- results$state
gain_ss<- results$gain[dim(results$gain)[1],,]
m <- dim(gain_ss)[1]
W_1 <-(solve(I(m)-(I(m)-gain_ss%*%Z)%*%T)%*%gain_ss)[1,]
delta_hat <- W_1%*%rowMeans(Y)      
delta_c <- results$state[1,]

C <- recover_level_C(delta_c,C0 = 1, delta_hat = delta_hat)
#  adjust C
level_c_results <- data.frame('C'=C,'date'=as.Date(data_filter$date))
normalization <- as.numeric( filter(level_c_results,date=='1967-01-01') %>% select(C))
level_c_results$C <- level_c_results$C/normalization*100
level_c_results <- level_c_results %>% filter(date>='1961-01-01')
plt <- ggplot(level_c_results,aes(x=date,y=C))+
  geom_line()+
  scale_x_date(date_labels="%y",date_breaks  ="2 year")+
  theme_classic()+
  xlab('Year')+
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_line())

plt
ggsave(plot=plt,here('C_replicate_paper.png'),width = 12,height = 8,dpi=200)
  
# treat data for part b ------------------------------------------------------------

# loading data from fred
emp <- fread(here('Data','PAYEMS.csv'))
ip <- fread(here('Data','INDPRO.csv'))
per_income <- fread(here('Data','W875RX1.csv'))
mtp <- fread(here('Data','CMRMTSPL.csv'))
data_fred <- list(ip, per_income,mtp,emp) %>% reduce(full_join, by = "DATE")
colnames(data_fred) <- c('date','ip','personal_income','mtp','emp')
data_fred <- data_fred %>% arrange(ymd(data_fred$date))

# apply log to data
data_fred <- data_fred %>% mutate_at(c(2:5),.funs=log)

# first difference
data_fred <- data_fred %>% mutate_at(c(2:5),.funs = function(x) x-lag(x))
data_fred_filter <- data_fred[ date<='2019-12-01'& date>='1967-02-01',]
# remove mean
data_fred_filter[,2:5] <- sweep(data_fred_filter[,2:5] ,2,colMeans(data_fred_filter[,2:5],na.rm = T))

y_fred <- t(as.matrix(data_fred_filter[,2:5]))
write.table(y_fred,file=here('Mlab','y_fred.csv'),row.names = F)

