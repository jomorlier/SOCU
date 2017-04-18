xvals <- 1:100
yvals <- rnorm(100)
data <- data.frame(xvals,yvals)

library(manipulate)
library(ggplot2)
manipulate({
  #define plotting function 
  ggplot(data,aes(xvals,yvals)) +
    geom_smooth(method="loess",span=span.val) +
    browser()
    geom_point()},
  #define variable that will be changed in plot
  span.val=slider(0.1,1)
)