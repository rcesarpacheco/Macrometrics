library(here)
library(doParallel)
library(doRNG)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(stargazer)
library(clipr)
################################
#                              #
#          Problem 2           #
#                              #
################################
set.seed(6413)
# function to draw from inverse gamma distribution
inv_gamma_draw <- function(a,b,ndraws) {
  shape <- a
  scale <- b
  draws <- rgamma(ndraws, shape=shape, scale = scale)
  return(1/draws)
}

# function to draw discrete r.v. assuming values in vec_values with specified probs
# must be order from small to largest

discrete_draw <- function(values,probs) {
  cum_prob <- cumsum(probs)
  u <- runif(n=1)
  idx <- 1
  while (u>cum_prob[idx]) {
    idx <- idx+1
  }
  # return(sample(values, size=1, prob = probs))
  return(values[idx])
}

# function to generate the data
generate_data_q2 <- function(theta_1,theta_2,beta_1,beta_2,var_1,var_2,lambda,T) {
  x <- rnorm(T)
  y <- numeric(T)
  for (i in 1:lambda) {
    y[i] <- rnorm(1,theta_1+theta_2*x[i],sqrt(var_1))
  }
  for (i in (lambda+1):T) {
    y[i] <- rnorm(1,beta_1+beta_2*x[i],sqrt(var_2))
  }
  return(list('y'=y,'x'=x))
}

# generate random numbers from a bivariate normal dist
rand_biv_normal <- function(mu,cov_matrix,n_draws,burn_in) {
  mu_1 <- mu[1]
  mu_2 <- mu[2]
  sigma_1 <- sqrt(cov_matrix[1,1])
  sigma_2 <- sqrt(cov_matrix[2,2])
  rho <- cov_matrix[1,2]/sigma_1/sigma_2
  n_gibbs <- n_draws+burn_in
  x <- matrix(0,n_gibbs,2)
  for (i in 2:n_gibbs) {
    x[i,1] <- rnorm(1,mean = mu_1+sigma_1/sigma_2*rho*(x[i-1,2]-mu_2),sigma_1*sqrt(1-rho^2))
    x[i,2] <- rnorm(1,mean = mu_2+sigma_2/sigma_1*rho*(x[i,1]-mu_1),sigma_2*sqrt(1-rho^2))
  }
  return(x[(burn_in+1):n_gibbs,])
}

# generate data
true_theta <- c(2,1)
true_beta <- c(1.5,0.8)
true_lambda <- 85
true_var_1 <- 0.2
true_var_2 <- 0.5

data <- generate_data_q2(theta_1 = true_theta[1],
                         theta_2 = true_theta[2],
                         beta_1 = true_beta[1],
                         beta_2=true_beta[2],
                         var_1 = true_var_1,
                         var_2 = true_var_2,lambda = true_lambda,T = 500)

# priors
V_theta <- diag(100,nrow = 2)
V_beta <- diag(100,nrow = 2)
mu_theta <- rep(0,2)
mu_beta <- rep(0,2)
a1 <- 3
b1 <- 3
a2 <- 0.5
b2 <- 0.5

N <- length(data$y)
S <- 10000
B <- 1000
N_gibbs <- S+B

theta <- matrix(nrow = N_gibbs,ncol=2)
beta <- matrix(nrow = N_gibbs,ncol=2)
lambda <- numeric(N_gibbs)
var_1 <- numeric(N_gibbs)
var_2 <- numeric(N_gibbs)
lambda_values <- 1:(N-1)

# initial values. Prior mean
lambda[1] <- 250
var_1[1] <- 1/((a1-1)*a2)
var_2[1] <- 1/((b1-1)*b2)
for (i in 2:N_gibbs) {
  X_theta <- cbind(1,data$x[1:lambda[i-1]])
  X_beta <- cbind(1,data$x[(lambda[i-1]+1):N])
  Y_theta <- data$y[1:lambda[i-1]]
  Y_beta <- data$y[(lambda[i-1]+1):N]
  
  theta_hat <- solve(t(X_theta)%*%X_theta)%*%t(X_theta)%*%Y_theta
  vcov_theta_posterior <- solve(t(X_theta)%*%X_theta/var_1[i-1]+solve(V_theta))
  mu_theta_posterior <- vcov_theta_posterior%*%(t(X_theta)%*%X_theta%*%theta_hat/var_1[i-1]+solve(V_theta)%*%mu_theta)
  theta[i,] <- rand_biv_normal(mu_theta_posterior,vcov_theta_posterior,1,100)
  
  beta_hat <- solve(t(X_beta)%*%X_beta)%*%t(X_beta)%*%Y_beta
  vcov_beta_posterior <- solve(t(X_beta)%*%X_beta/var_2[i-1]+solve(V_beta))
  mu_beta_posterior <- vcov_beta_posterior%*%(t(X_beta)%*%X_beta%*%beta_hat/var_2[i-1]+solve(V_beta)%*%mu_beta)
  beta[i,] <- rand_biv_normal(mu_beta_posterior,vcov_beta_posterior,1,100)
  
  a1_posterior <- a1+lambda[i-1]/2
  a2_posterior <- 1/(1/a2+1/2*sum((Y_theta-X_theta%*%theta[i,])^2))
  b1_posterior <- b1+(N-lambda[i-1])/2
  b2_posterior <- 1/(1/b2+1/2*sum((Y_beta-X_beta%*%beta[i,])^2)) 
  var_1[i] <- inv_gamma_draw(a1_posterior,a2_posterior,1)
  var_2[i] <- inv_gamma_draw(b1_posterior,b2_posterior,1)
  
  # update the probabilities of lambda
  probs <- numeric(length(lambda_values))
  for (idx_lambda in 1:length(lambda_values)) {
    X_theta_temp <- cbind(1,data$x[1:lambda_values[idx_lambda]])
    X_beta_temp <- cbind(1,data$x[(lambda_values[idx_lambda]+1):N])
    Y_theta_temp <- data$y[1:lambda_values[idx_lambda]]
    Y_beta_temp <- data$y[(lambda_values[idx_lambda]+1):N]
    probs[idx_lambda] <- var_1[i]^(-idx_lambda/2)*exp(-1/(2*var_1[i])*sum((Y_theta_temp-X_theta_temp%*%theta[i,])^2))*var_2[i]^(-(T-idx_lambda)/2)*exp(-1/(2*var_2[i])*sum((Y_beta_temp-X_beta_temp%*%beta[i,])^2))
  }
  probs <- probs/sum(probs)
  lambda[i] <- discrete_draw(lambda_values,probs)
}

# results discarding burn in
colMeans(theta[(B+1):N,])
colMeans(beta[(B+1):N,])
mean(var_1[(B+1):N])
mean(var_2[(B+1):N])
mean(lambda[(B+1):N])

################################
#                              #
#          Problem 3           #
#                              #
################################

n <- 100
e <- rnorm(n)
x <- rnorm(n)
X <- cbind(rep(1,n),x)
beta_true <- 0.5
y_star <- beta_true*x+e
y <- 1*(y_star>0)

# generate random numbers from a normal distribution truncated at a and b, a<b
rtrun_norm <- function(lower_bound,upper_bound,mu,sigma,n_draws) {
    u <- runif(n_draws)
    Fa <- pnorm(lower_bound,mean = mu,sd = sigma)
    Fb <- pnorm(upper_bound,mean = mu,sd = sigma)
    tn <- qnorm(u*(Fb-Fa)+Fa,mean = mu,sd = sigma)
  return(tn)
}

gibbs_sampling <- function(X,y,burn_in,S) {
  y_star <- y
  n <- length(y)
  n_gibbs <- burn_in+S
  xtx_inv <-solve(t(X)%*%X)
  beta <- matrix(NA,nrow = 2,ncol = n_gibbs)
  prob <- numeric(n_gibbs)
  for (iter_gibs in 1:n_gibbs) {
    beta_hat <- xtx_inv%*%t(X)%*%y_star
    beta[,iter_gibs] <- rand_biv_normal(mu = beta_hat,cov_matrix = xtx_inv,n_draws = 1,burn_in = 100)
    for (i in 1:n) {
      if (y[i]==0) {
        y_star[i] <- rtrun_norm(-Inf,0,X[i,]%*%beta[,iter_gibs],1,1)
      }
      else{
        y_star[i] <- rtrun_norm(0,Inf,X[i,]%*%beta[,iter_gibs],1,1)
      }
    }
  prob[iter_gibs] <- pnorm(c(1,2)%*%beta[,iter_gibs])
  }
  return(list('betas'= beta[,(burn_in+1):n_gibbs],'prob'=prob[(burn_in+1):n_gibbs]))
}

results_gibbs <- gibbs_sampling(X,y,1000,10000)
betas_gibbs <- rowMeans(results_gibbs$betas)
prob_gibbs <- mean(results_gibbs$prob)
sd_prob <- sd(results_gibbs$prob)
sd_betas <-apply(results_gibbs$betas,  MARGIN = 1,sd)


df_results <- data.frame('Variable' = c("beta_0", "beta_{1}","Phi(beta_0+2\beta_1)"),
                         'Posterior Mean' = c(betas_gibbs,prob_gibbs),
                         'Posterior St. Dev.' = c(sd_betas,sd_prob),
                         check.names=FALSE)
clipr::write_clip(df_results)       


# # Monte Carlo to check if the mean is centered properly
B <- 1000
registerDoParallel(cores=8)
results_monte_carlo <- foreach(icount(B),.options.RNG=123, .combine=rbind) %dorng% {
  n <- 100
  e <- rnorm(n)
  x <- rnorm(n)
  X <- cbind(rep(1,n),x)
  beta_true <- 0.5
  y_star <- beta_true*x+e
  y <- 1*(y_star>0)
  results_gibbs <- gibbs_sampling(X,y,1000,10000)
  betas_monte_carlo <- rowMeans(results_gibbs$betas)
  prob_gibbs <- mean(results_gibbs$prob)
  c(betas_monte_carlo,prob_gibbs)
}
save.image(here('r_workspace'))


plt1 <- ggplot(data.frame('Estimates'=results_monte_carlo[,1]),aes(x=Estimates,y=..density..))+
  geom_histogram()+
  geom_vline(aes(xintercept=0),
             color="blue", linetype="dashed")+
  ggtitle(TeX('Histogram $\\hat{\\beta_0}$'))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

plt2 <- ggplot(data.frame('Estimates'=results_monte_carlo[,2]),aes(x=Estimates,y=..density..))+
  geom_histogram()+
  geom_vline(aes(xintercept=0.5),
             color="blue", linetype="dashed")+
  ggtitle(TeX('Histogram $\\hat{\\beta_1}$'))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

plt3 <- ggplot(data.frame('Estimates'=results_monte_carlo[,3]),aes(x=Estimates,y=..density..))+
  geom_histogram()+
  geom_vline(aes(xintercept=pnorm(1)),
             color="blue", linetype="dashed")+
  ggtitle(TeX('Histogram $\\Phi(\\hat{\\beta_0}+2\\hat{\\beta_1})$'))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
plt4 <- ggarrange(plt1,plt2,plt3,ncol=3 )
plt4
ggsave(plot=plt4,here('histogram.png'),width = 14,height = 8,dpi=200)

load(here('r_workspace'))