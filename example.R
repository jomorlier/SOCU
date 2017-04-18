
library("DiceKriging")
library("DiceOptim")

 x <- c(0, 0.4, 0.6, 0.8, 1)
 x<-c(-0.8879654,  0.9625596, 0.2956500)
 y <- 10 * c(-0.6, 0, -2, 0.5, 0.9)
 y<-c(8.519378, 9.485646, 3.611862)
 theta <- 0.1
 sigma <- 10
 trend <- 5 * c(-2, 1)
 #model <- km(~x, design = data.frame(x = x), response = y,
 #               coef.trend = trend, covtype = "gauss", coef.cov = theta,
 #               coef.var = sigma^2)
 
 model <- km(~x, design = data.frame(x = x), response = y,
              covtype = "gauss")
 
 
  t <- seq(from = 0, to = 1, by = 0.005)
  p <- predict(model, newdata = data.frame(x = t), type = "UK")
  EI_values <- apply(as.matrix(t), 1, EI, model, type = "UK")
plot(t,p$mean,type="l",col="blue",xlim=c(0,1),ylim=c(-40,30))
points(x,y,col="red")
lines(t,p$lower95,lty=2)
lines(t,p$upper95,lty=2)
abline(h=0)

EI_values <- apply(as.matrix(t), 1, EI, model, type = "UK")
plot(t,EI_values,type="l")

x_star <- max_EI(model, lower = 0, upper = 1)