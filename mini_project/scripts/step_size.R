library(ggplot2)

# look at x in [-5, 5]
x = seq(-5,5,0.01)

x0 = 2.5

alpha1 = 1
alpha2 = 0.1

square = function(x)
{
 return(x^2)
}

gradient = function(x)
{
 return(2*x)
}


new_iterate = function(iterate, gradient_function, step_size)
{
 return( iterate - step_size * gradient_function(iterate) )
}


x1_1 = new_iterate(x0, gradient, alpha1)
x1_1

x1_2 = new_iterate(x0, gradient, alpha2)
x1_2

fx = data.frame(x = x, x2 = square(x))
test_df = data.frame(cbind(x0, square(x0)))

ggplot(df, aes(x, V2)) +
 geom_line() +
 geom_point(data=test_df, aes(x0, V2, color="blue")) +
 scale_color_manual(values=c("blue"="#0000FF"))


