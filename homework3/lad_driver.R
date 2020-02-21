###
#
# Jimmy Hickey
#
###

library(jhickeyST790)
library(ggplot2)

## Number of International Calls from Belgium,
## taken from the Belgian Statistical Survey,
## published by the Ministry of Economy,
##
## 73 subjects, 2 variables:
##  Year(x[i])
##  Number of Calls (y[i], in tens of millions)
##
## http://www.uni-koeln.de/themen/statistik/data/rousseeuw/
## Datasets used in Robust Regression and Outlier Detection (Rousseeuw and Leroy, 1986).
## Provided on-line at the University of Cologne.

x <- c(50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68,
       69, 70, 71, 72, 73)

y <- c(0.44, 0.47, 0.47, 0.59, 0.66, 0.73, 0.81, 0.88, 1.06, 1.20, 1.35, 1.49, 1.61,
       2.12, 11.90, 12.40, 14.20, 15.90, 18.20, 21.20, 4.30, 2.40, 2.70, 2.90)

x = as.matrix(x)

beta0 = c(0.1)
epsilon = 0.25

output_smLAD = smLAD(y=y, X=x, beta=beta0, epsilon=epsilon, max_iter=1e2, tol=1e-3)

telephone_output = data.frame(cbind(
                                    "smooth" = x %*% output_smLAD$final_iterate,
                                    "x" = x[,1],
                                    "y" = y
                                    ))

plot(output_smLAD$objective_history)

ggplot(telephone_output) +
  geom_point(size=2, shape=23, aes(x = x, y = y)) +
  geom_line(aes(x = x, y = V1))


####
# Newton
###

set.seed(12345)

## Data set 1
n <- 200
p <- 300
X1 <- matrix(rnorm(n*p),n,p)
beta0 <- matrix(rnorm(p),p,1)
y1 <- X1%*%beta0 + rnorm(n)
lambda1 = 10
epsilon1 = 0.25

naive1 = system.time(
  lad_newton(y = y1,
             X=X1,
             beta = beta0,
             epsilon=epsilon,
             lambda=lambda,
             naive=TRUE,
             max_iter=1e2,
             tol=1e-3))[[3]]

smw1 = system.time(
  lad_newton(y = y1,
             X=X1,
             beta = beta0,
             epsilon=epsilon,
             lambda=lambda,
             naive=FALSE,
             max_iter=1e2,
             tol=1e-3))[[3]]

## Data set 2
p <- 600
X2 <- matrix(rnorm(n*p),n,p)
beta0 <- matrix(rnorm(p),p,1)
y2 <- X2%*%beta0 + rnorm(n)

naive2 = system.time(
  lad_newton(y = y2,
             X=X2,
             beta = beta0,
             epsilon=epsilon,
             lambda=lambda,
             naive=TRUE,
             max_iter=1e2,
             tol=1e-3))[[3]]

smw2 = system.time(
  lad_newton(y = y2,
             X=X2,
             beta = beta0,
             epsilon=epsilon,
             lambda=lambda,
             naive=FALSE,
             max_iter=1e2,
             tol=1e-3))[[3]]


## Data set 3
p <- 1200
X3 <- matrix(rnorm(n*p),n,p)
beta0 <- matrix(rnorm(p),p,1)
y3 <- X3%*%beta0 + rnorm(n)

naive3 = system.time(
  lad_newton(y = y3,
             X=X3,
             beta = beta0,
             epsilon=epsilon,
             lambda=lambda,
             naive=TRUE,
             max_iter=1e2,
             tol=1e-3))[[3]]

smw3 = system.time(
  lad_newton(y = y3,
             X=X3,
             beta = beta0,
             epsilon=epsilon,
             lambda=lambda,
             naive=FALSE,
             max_iter=1e2,
             tol=1e-3))[[3]]


timing = data.frame(
    "naive" = rbind(naive1, naive2, naive3),
    "smw" = rbind(smw1, smw2, smw3),
    "p" = rbind(300, 600, 1200)
)

ggplot(timing) +
  geom_line(aes(x = p, y = naive, color = "naive")) +
  geom_line(aes(x = p, y = smw, color = "smw")) +
  scale_color_manual(values=c("naive"="#FF0000", "smw"="#0000FF")) +
  ylab("time (s)") +
  xlab("number of columns") +
  ggtitle("Comparison of naive and SMW Newton Method")

