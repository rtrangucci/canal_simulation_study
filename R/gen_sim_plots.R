source("R/fun_sim_data.R")

dat <- gen_data(S=100, G=3, seed = 123, int = 0.05, n=564, l=0.1,house_dist = 'clustered_uniform')
grid_x <- seq(0,dat$M,length.out = 100)
grid_y <- seq(0,dat$C,length.out = 100)
grid_df <- expand.grid(x = grid_x, y = grid_y)
plt_grid <- matrix(0, length(grid_x), length(grid_y))
for (i in 1:length(grid_x))
  for (j in 1:length(grid_y))
    plt_grid[i,j] <- 1 - exp(- (dat$int + dat$total_risk_bind(grid_x[i], grid_y[j]))*1)
contour(grid_x, grid_y, plt_grid, main = expression(paste('True P(',Y[ij]==1,'|s) of disease')), xlab = 'x-location (kilometers)', ylab = 'y-location (kilometers)')
abline(h = 0, col='red')
abline(h = 2*dat$C/3, col='red')
abline(v = M/2, col='red')
points(5,4.1,pch='x',col='blue')
points(10.1,8/3,pch='x',col='blue')
points(10.1,0,pch='x',col='blue')
arrows(4.3,4,4.3,3.2,col='blue',length=0.1)
arrows(4,3.2,0.5,3.2,col='blue',length=0.1)
arrows(10,3.2,6,3.2,col='blue',length=0.1)
arrows(4.3,2,4.3,0.7,col='blue',length=0.1)
arrows(4,0.5,0.5,0.5,col='blue',length=0.1)
arrows(10,0.5,6,0.5,col='blue',length=0.1)

plot(dat$data$coords[,1],dat$data$coords[,2], main = 'Household locations', ylab = 'y-location (kilometers)', xlab = 'x-location (kilometers)')
abline(h = 0, col='red')
abline(h = 2 * dat$C/3, col='red')
abline(v = dat$M/2, col='red')

boundary_points <- matrix(
  c(0,0,
    5,0,
    10,0,
    0,8/3,
    5,8/3,
    10,8/3,
    5,4
  ), byrow=T, 
  7, 2)

df <- data.frame(
  x = boundary_points[,1],
  y = boundary_points[,2],
  lab = 1:7
)

plot(df$x, df$y, main = 'Gaussian process prior construction', ylab = 'y-location (kilometers)', xlab = 'x-location (kilometers)',xlim=c(0,10),ylim=c(0,4),
     pch = 'x', col = 'black')
with(df[c(1,3,4,6),], text(x, y, lab, pos = 3))
with(df[2,], text(x+0.4, y, bquote(p[x%*%y]), pos = 3))
with(df[5,], text(x+0.4, y, bquote(p[x[2]%*%y]), pos = 3))
with(df[2,], text(x+0.4, y+0.4, bquote(eta[x%*%y]), pos = 3))
with(df[5,], text(x+0.4, y+0.4, bquote(eta[x[2]%*%y]), pos = 3))
with(df[7,], text(x, y, lab, pos = 4))
abline(h = 0, col='red')
abline(h = 2 * dat$C/3, col='red')
abline(v = dat$M/2, col='red')
arrows(4.3,4,4.3,3.2,col='blue',length=0.1)
arrows(4,3.2,0.5,3.2,col='blue',length=0.1)
arrows(10,3.2,6,3.2,col='blue',length=0.1)
arrows(4.3,2,4.3,0.7,col='blue',length=0.1)
arrows(4,0.5,0.5,0.5,col='blue',length=0.1)
arrows(10,0.5,6,0.5,col='blue',length=0.1)
text(2,3.4,expression(eta[x[2]]))
text(8,3.4,expression(eta[x[2]]))
text(2,0.7,expression(eta[x]))
text(8,0.7,expression(eta[x]))
text(3.9,1.35,expression(eta[y]))
text(3.9,3.5,expression(eta[y]))
