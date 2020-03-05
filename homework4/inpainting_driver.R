###
#
# Jimmy Hickey
# 2020-03-05
#
# Run inpainting code
#
###


ruscha = read.csv("~/git/jhickeyST790/homework4/data/Ruscha_in.csv")
cox = read.csv("~/git/jhickeyST790/homework4/data/gertrude_cox_in.csv")

r_mat = as.matrix(ruscha)
image(r_mat)

c_mat = as.matrix(cox)
image(c_mat)
