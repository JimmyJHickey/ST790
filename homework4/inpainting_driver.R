###
#
# Jimmy Hickey
# 2020-03-05
#
# Run inpainting code
#
###
library(jhickeyST790)

ruscha = read.csv("~/git/jhickeyST790/homework4/data/Ruscha_in.csv", header=FALSE)
cox = read.csv("~/git/jhickeyST790/homework4/data/gertrude_cox_in.csv", header=FALSE)



r_mat = as.matrix(ruscha)
image(r_mat)

r_final = myInpaint(r_mat)
image(matrix(r_final$x[1:(nrow(r_mat) * ncol(r_mat))],
                nrow = nrow(r_mat),
                ncol = ncol(r_mat)))

c_mat = as.matrix(cox)
image(c_mat)

c_final = myInpaint(c_mat)
image(matrix(c_final$x[1:(nrow(c_mat) * ncol(c_mat))],
                nrow = nrow(c_mat),
                ncol = ncol(c_mat)))
