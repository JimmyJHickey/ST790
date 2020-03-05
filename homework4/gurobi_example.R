###
#
# Jimmy Hickey
# 2020-03-03
#
###

# install gurobi
install.packages('/Library/gurobi901/mac64/R/gurobi_9.0-1_R_3.6.1.tgz', repos=NULL)
install.packages("slam", repos = "https://cloud.r-project.org")
#install.packages('/Users/jimmy/Documents/gurobi_9.0-1_R_3.6.1.tgz', repos=NULL)

install.packages('/Library/gurobi811/mac64/R/gurobi_8.1-1_R_3.5.0.tgz', repos=NULL)



# install.packages("slam", repos = "https://cloud.r-project.org")


library(gurobi)


# Copyright 2019, Gurobi Optimization, LLC
#
# This example formulates and solves the following simple MIP model:
#  maximize
#        x +   y + 2 z
#  subject to
#        x + 2 y + 3 z <= 4
#        x +   y       >= 1
#        x, y, z binary

model <- list()

model$A          <- matrix(c(1,2,3,1,1,0), nrow=2, ncol=3, byrow=T)
model$obj        <- c(1,1,2)
model$modelsense <- 'max'
model$rhs        <- c(4,1)
model$sense      <- c('<', '>')
model$vtype      <- 'B'

params <- list(OutputFlag=0)

result <- gurobi(model, params)

print('Solution:')
print(result$objval)
print(result$x)
``
# Clear space
rm(model, result, params)
