library(jhickeyST790)

# call plot restrictions on the logdet function
plot_restrictions( fx = logdet,
                   rx = restrictions_positive_definite_matrices,
                   nRow = 3,
                   nCol = 3)
