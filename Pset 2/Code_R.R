library(data.table)
library(ggplot2)
library(latex2exp)
library(here)
library(quantmod)
library(expm)
library(lubridate)
library(bayesforecast)
library(ggpubr)
library(forecast)
################################
#                              #
#          Problem 1           #
#                              #
################################

alpha <- 0.5
theta <- 0.8

psi <- c()
psi[1] <- 1
n_psi <- 20
sigma2 <- 1

for (i in 2:n_psi) {
  psi[i] <- alpha^(i-2)*(alpha+theta)
}


plt <- ggplot(data.frame(x=0:(length(psi)-1),psi),aes(x=x,y=psi))+
  geom_point()+
  xlab('j')+
  ylab(TeX('$ \\psi_j$ '))+
  scale_x_continuous(breaks = 0:(length(psi)-1))+
  theme_classic()+
  theme(panel.grid.major = element_line())
ggsave(here('psi_j.png'),width = 10,height = 7,dpi=200)


# part e
autocov <- function(lag,alpha,theta,sigma2) {
  if (lag==0) {
    return (sigma2*(1+(alpha+theta)^2/(1-alpha^2)))
  }
  else
    return (sigma2*alpha^(lag-1)*(alpha+theta)*(1+alpha*(alpha+theta)/(1-alpha^2)))
}

# part f
T <- 10
e <- rep(0,T)
e[2] <- 1
y <- rep(0,T)
for (i in 2:T) {
  y[i] <- alpha*y[i-1]+e[i]+theta*e[i-1]
}
t <- 
plt <- ggplot(data.frame(t=-1:(length(y)-2),y),aes(x=t,y=y))+
  geom_point()+
  xlab('t')+
  ylab(TeX('$ \\y_t$ '))+
  scale_x_continuous(breaks = -1:(length(y)-2))+
  theme_classic()+
  theme(panel.grid.major = element_line())
plt
ggsave(here('y_t.png'),width = 10,height = 7,dpi=200)


# Problem 3 ---------------------------------------------------------------


################################
#                              #
#          Problem 3           #
#                              #
################################


gdp <- getSymbols('GDPC1',src='FRED',env = NULL)
gdp_growth <- xts(gdp/lag(gdp)-1)
gdp_growth <- gdp_growth['1961/']
gdp_growth_filtered <- gdp_growth['1961/2015']

x <- as.matrix(cbind(lag(gdp_growth_filtered,1),lag(gdp_growth_filtered,2),lag(gdp_growth_filtered,3),lag(gdp_growth_filtered,4)))
x <- na.omit(x)
gdp_growth_filtered <- last(gdp_growth_filtered,length(gdp_growth_filtered)-4)
betas <- solve(t(x)%*%x)%*%t(x)%*%gdp_growth_filtered

data <- data.frame(gdp_growth_filtered,x)
colnames(data) <- c('y','y1','y2','y3','y4')
ols <- lm(y~ -1 + y1+y2+y3+y4,data)
summary(ols)

# create companion matrix

A <- rbind(c(betas[1],betas[2],betas[3],betas[4]),c(1,0,0,0),c(0,1,0,0),c(0,0,1,0))

# forecast function
forecast_arp <- function(A,y,h,sigma2,level_ic) {
  p <- dim(A)[1]
  Y_t <- as.numeric(last(y,p))
  Y_hat <- (A%^%h)%*%Y_t
  mse <- 0
  for (j in 0:(h-1)) {
    mse <- mse+(A%^%j)%*%t(A%^%j)
  }
  mse <- sigma2*mse
  
  ic_sup <- Y_hat[1]+qnorm(1-(1-level_ic)/2)*sqrt(mse[1,1])
  ic_inf <- Y_hat[1]-qnorm(1-(1-level_ic)/2)*sqrt(mse[1,1])
  
  return(list('forecast'=Y_hat[1],'mse'=mse[1,1],ic=c(ic_inf,ic_sup)))
}

# part i
level <- 0.95
forecast_arp(A,gdp_growth_filtered,1,var(residuals(ols)),level)
forecast_arp(A,gdp_growth_filtered,2,var(residuals(ols)),level)

#part ii)
forecasts_df <- data.frame(date=Date(0),
                           forecast_date1=Date(),
                           forecast_date2=Date(0),
                           forecast1=numeric(0),
                           forecast1_ic_inf=numeric(0),
                           forecast1_ic_sup=numeric(0),
                           forecast2=numeric(0),
                           forecast2_ic_inf=numeric(0),
                           forecast2_ic_sup=numeric(0),
                           realized_y=numeric(0))
end_dates <- seq.Date(as.Date("2016/1/1"), as.Date("2019/10/1"), by = "quarter")
for (i in 1:length(end_dates)) {
  gdp_growth_filtered <- window(gdp_growth,end = end_dates[i])
  x <- as.matrix(cbind(lag(gdp_growth_filtered,1),lag(gdp_growth_filtered,2),lag(gdp_growth_filtered,3),lag(gdp_growth_filtered,4)))
  x <- na.omit(x)
  gdp_growth_filtered <- last(gdp_growth_filtered,length(gdp_growth_filtered)-4)
  data <- data.frame(gdp_growth_filtered,x)
  colnames(data) <- c('y','y1','y2','y3','y4')
  ols <- lm(y~ -1 + y1+y2+y3+y4,data)
  betas <- coefficients(ols)
  A <- rbind(c(betas[1],betas[2],betas[3],betas[4]),c(1,0,0,0),c(0,1,0,0),c(0,0,1,0))
  date_forecast_1 <- end_dates[i]+months(3)
  date_forecast_2 <- end_dates[i]+months(6)
  forecast_res1 <- forecast_arp(A,gdp_growth_filtered,1,var(residuals(ols)),level)
  forecast_res2 <- forecast_arp(A,gdp_growth_filtered,2,var(residuals(ols)),level)
  
  results <- list('date'=as.Date(end_dates[i]),
               'forecast_date1'=as.Date(date_forecast_1),
               'forecast_date2'=as.Date(date_forecast_2),
               'forecast1'=forecast_res1$forecast,
               'forecast1_ic_inf'=forecast_res1$ic[1],
               'forecast1_ic_sup'=forecast_res1$ic[2],
               'forecast2'=forecast_res2$forecast,
               'forecast2_ic_inf'=forecast_res2$ic[1],
               'forecast2_ic_sup'=forecast_res2$ic[2],
               'realized_y1'=gdp_growth[date_forecast_1],
               'realized_y2'=gdp_growth[date_forecast_2])
  forecasts_df <-  rbind(forecasts_df,results)
}
forecasts_df[1:3] <- lapply(forecasts_df[1:3],FUN = as.Date)
forecasts_df$err1 <- forecasts_df$forecast1 - forecasts_df$realized_y1
forecasts_df$err2 <- forecasts_df$forecast2 - forecasts_df$realized_y2
e1 <- forecasts_df$err1
plt1 <- ggplot(forecasts_df,aes(x=forecast_date1))+
  geom_line(aes(y=forecast1,color='1 Step Forecast'))+
  geom_point(aes(y=forecast1,color='1 Step Forecast',shape='1 Step Forecast'))+
  geom_ribbon(aes(ymin = forecast1_ic_inf, ymax = forecast1_ic_sup), alpha = 0.2,fill='red')+
  geom_line(aes(y=realized_y1,color='Actual'))+
  geom_point(aes(y=realized_y1,color='Actual',shape='Actual'))+
  labs(x="Date", y="GDP growth", color = "Legend")+
  scale_color_manual(name='',labels=c('1 Step Forecast','Actual'),values=c('red','black'))+
  scale_shape_manual(name='',labels=c('1 Step Forecast','Actual'),values=c(20,15))+
  scale_x_date(breaks=forecasts_df$forecast_date1)+  
  theme_classic()+
  ggtitle('1 step Forecast')+
  theme(panel.grid.major = element_line(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5))

plt2 <- ggplot(forecasts_df,aes(x=forecast_date2))+
  geom_line(aes(y=forecast2,color='2 Step Forecast'))+
  geom_point(aes(y=forecast2,color='2 Step Forecast',shape='2 Step Forecast'))+
  geom_ribbon(aes(ymin = forecast2_ic_inf, ymax = forecast2_ic_sup), alpha = 0.2,fill='red')+
  geom_line(aes(y=realized_y2,color='Actual'))+
  geom_point(aes(y=realized_y2,color='Actual',shape='Actual'))+
  labs(x="Date", y="GDP growth", color = "Legend")+
  scale_color_manual(name='',labels=c('2 Step Forecast','Actual'),values=c('red','black'))+
  scale_shape_manual(name='',labels=c('2 Step Forecast','Actual'),values=c(20,15))+  
  scale_x_date(breaks=forecasts_df$forecast_date2)+  
  theme_classic()+
  ggtitle('2 step Forecast')+
  theme(panel.grid.major = element_line(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5))

plt <- ggarrange(plt1,plt2,nrow=2)
ggsave(plot=plt,here('forecasts.png'),width = 10,height = 8,dpi=200)

plt3 <- ggacf(forecasts_df$err1)+
  theme_classic()+
  ggtitle('Correlogram 1 step forecast')+
  theme(plot.title = element_text(hjust = 0.5))


plt4 <- ggacf(forecasts_df$err2)+
  ggtitle('Correlogram 2 step forecast')+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

  
plt5 <- ggarrange(plt3,plt4,nrow=2)

ggsave(plot=plt5,here('correlogram.png'),width = 10,height = 8,dpi=200)

# iv)
# compare the model with an AR(5)


forecasts_model2_df <- data.frame(date=Date(0),
                          forecast_date1=Date(),
                          forecast1=numeric(0),
                          forecast1_ic_inf=numeric(0),
                          forecast1_ic_sup=numeric(0),
                          realized_y=numeric(0))
for (i in 1:length(end_dates)) {
  gdp_growth_filtered <- window(gdp_growth,end = end_dates[i])
  x <- as.matrix(cbind(lag(gdp_growth_filtered,1),lag(gdp_growth_filtered,2),lag(gdp_growth_filtered,3),lag(gdp_growth_filtered,4),lag(gdp_growth_filtered,5)))
  x <- na.omit(x)
  gdp_growth_filtered <- last(gdp_growth_filtered,length(gdp_growth_filtered)-5)
  data <- data.frame(gdp_growth_filtered,x)
  colnames(data) <- c('y','y1','y2','y3','y4','y5')
  ols <- lm(y~ -1 + y1+y2+y3+y4+y5,data)
  betas <- coefficients(ols)
  A <- rbind(c(betas[1],betas[2],betas[3],betas[4],betas[5]),c(1,0,0,0,0),c(0,1,0,0,0),c(0,0,1,0,0),c(0,0,0,1,0))
  date_forecast_1 <- end_dates[i]+months(3)
  forecast_res1 <- forecast_arp(A,gdp_growth_filtered,1,var(residuals(ols)),level)
  
  results <- list('date'=as.Date(end_dates[i]),
                  'forecast_date1'=as.Date(date_forecast_1),
                  'forecast1'=forecast_res1$forecast,
                  'forecast1_ic_inf'=forecast_res1$ic[1],
                  'forecast1_ic_sup'=forecast_res1$ic[2],
                  'realized_y1'=gdp_growth[date_forecast_1])
  forecasts_model2_df <-  rbind(forecasts_model2_df,results)
}

forecasts_model2_df[1:2] <- lapply(forecasts_model2_df[1:2],FUN = as.Date)
e2 <- forecasts_model2_df$err1 <- forecasts_model2_df$forecast1 - forecasts_model2_df$realized_y1
# dm.test from package for comparison
dm.test(e1, e2, alternative ="two.sided", h = 1,
        power = 2)
# manual calculation
d <- (forecasts_df$err1^2-forecasts_model2_df$err1^2)
P <- length(d)
acf_estimation <- Acf(d,plot = F,lag.max = 1,type = 'covariance')
acf_estimation <- acf_estimation$acf
long_run_variance <- sum(c(acf_estimation[1],2*acf_estimation[-1]))

S <- sqrt(P)*mean(d)/sqrt(long_run_variance)
2*(1-pnorm(S))
