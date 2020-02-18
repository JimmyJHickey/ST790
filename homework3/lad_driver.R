library(jhickeyST790)

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

fx_lad(y, x, beta0, epsilon = 0.25)

gradf_lad(y, x, beta0, epsilon = 0.25)

mm_update(y, x, beta0, epsilon = 0.25)
