# depending on where your code is
setwd("~/projects/EPS236-Env.-Modeling-/")

# function to create eigenvalue diagonal matrix Lambda
exp_Lambda = function(lambda_values, t) {
  return( diag(exp(lambda_values * t)) )
}

# diffusion
add_diffusion_term = function(K, k.h.dt, n_box) {
  for (i in seq(n_box)) {
    K[i, i] = K[i, i] + -2 * k.h.dt
    if (i == 1) {
      K[i, n_box] = K[i, n_box] + k.h.dt
    } else {
      K[i, (n_box+i-1)%%n_box] = K[i, (n_box+i-1)%%n_box] + k.h.dt
    }
    if (i == n_box - 1) {
      K[i, n_box] = K[i, n_box] + k.h.dt
    }
    K[i, (n_box+i+1)%%n_box] = K[i, (n_box+i+1)%%n_box] + k.h.dt
  }
  return( K )
}

# advection - call if v_a != 0
add_advection_term = function(K, v_a, n_box) {
  for (i in seq(n_box)) {
    K[i, i] = K[i, i] + -1 * v_a
    if (v_a > 0) {
      if (i == 1) {
        K[i, n_box] = K[i, n_box] + v_a
      } else {
        K[i, (n_box+i-1)%%n_box] = K[i, (n_box+i-1)%%n_box] + v_a
      }
    } else {
      if (i == n_box - 1) {
        K[i, n_box] = K[i, n_box] + v_a
      }
      K[i, (n_box+i+1)%%n_box] = K[i, (n_box+i+1)%%n_box] + v_a
    }
  }
  return( K )
}

circular_boxes = function(
                          # number of boxes
                          n_box,
                          # advection velocity
                          v_a,
                          # diffusion coeff
                          k.h.dt,
                          # initial box number
                          i_box,
                          # initial concentration
                          C_i,
                          # time series
                          t_max
                          ) {
  t = 0:t_max
  K = matrix(0, nrow = n_box, ncol = n_box)
  # construct transition matrix
  if (v_a != 0) {
    K = add_advection_term(K, v_a, n_box)
  }
  K = add_diffusion_term(K, k.h.dt, n_box)
  # eigenvalues and eigenvectors
  eigens = eigen(K)
  e_vals = eigens$values
  Xe = eigens$vectors
  
  # initial condition
  C_0 = matrix(0, n_box, 1)
  C_0[i_box] = C_i
  
  # matrix to store the transition of concentration
  Conc = matrix(nrow = n_box, ncol = length(t))
  Conc[, 1] = C_0
  # loop over time to calculate C(t) for t
  for (t_i in t) {
    if (t_i != 0) {
      Lambda_t = exp_Lambda(e_vals, t_i)
      conc_t_i = Re(Xe %*% Lambda_t %*% solve(Xe, C_0))
      Conc[, t_i+1] = conc_t_i
      filename = formatC(t_i, width = 4, format = "d", flag = "0")
      png(paste("./part1_movie/", filename, ".png", sep=""))
      plot(x = 1:n_box, y = conc_t_i, xlim = c(1, n_box), ylim = c(0, C_i), xlab = "#box", ylab = "concentration")
      dev.off()
    }
  }
  return( Conc )
}


# initial conditions 1
# number of boxes
n_box = 10
# advection velocity
v_a = .02
# diffusion coeff
k.h.dt = .05
# initial box number
i_box = 1
# initial concentration
C_i = 1
# time series
t_max = 1000
# produce plots
Conc = circular_boxes(n_box, v_a, k.h.dt, i_box, C_i, t_max)

# initial conditions 2
# number of boxes
n_box = 10
# advection velocity
v_a = 0
# diffusion coeff
k.h.dt = .05
# initial box number
i_box = 5
# initial concentration
C_i = 1
# time series
t_max = 1000
# produce plots
Conc = circular_boxes(n_box, v_a, k.h.dt, i_box, C_i, t_max)

# part (2)
# advection diffusion eqn

# initial conditions
# number of boxes
n_box = 50
# advection velocity
v_a = .014
# diffusion coeff
k.h.dt = .01
# initial box number
i_box = 25
# initial concentration
C_i = 1
# time series
t_max = 100

Conc = circular_boxes(n_box, v_a, k.h.dt, i_box, C_i, t_max)

advection_diffusion_eqn = function(D, x, t_i) {
  exp_term = exp( -(x - v_a*t_i)^2 / (4*D*t_i) )
  return(1 / sqrt(4*pi*D*t_i) * exp_term)
}

t_1 = 1:100
d_list = seq(0.01, 0.05, 0.0001)
didx = 1
x = seq(1, n_box) - i_box

ds = rep(NA, length(t_1))

for (t_i in t_1) {
  idx = 1
  sqrt_err_list = rep(NA, length(d_list))
  for (d_guess in d_list) {
    diff = Conc[, t_i] - advection_diffusion_eqn(d_guess, x, t_i)
    sqrt_err = sqrt(sum(diff^2))
    sqrt_err_list[idx] = sqrt_err
    
    idx = idx + 1
  }
  min_sqrt_err_idx = which(sqrt_err_list == min(sqrt_err_list))
  d_best = d_list[min_sqrt_err_idx]
  ds[didx] = d_best
  #plot(x, Conc[, t_i], ylim = c(0, C_i), ylab = paste("conc at time", t_i))
  didx = didx + 1
}

# best D for different time t_i
plot(t_1, ds, xlab = "time", ylab = "best D")
ds
abline(h = ds[length(ds)], col = "red")

# diffusion coeff D converges to 0.0160
d_best = ds[length(ds)]
plot(x, Conc[, 60], xlab = "x", ylab = "conc at time 60")
lines(x, advection_diffusion_eqn(d_best, x, t_1[60]))

# part (3)
# number of boxes
n_box = 10
# advection velocity
v_a = .014
# diffusion coeff
k.h.dt = .01
# initial box number
i_box = 5
# initial concentration
C_i = 1
# time series
t_max = 100
# produce plots
Conc = circular_boxes(n_box, v_a, k.h.dt, i_box, C_i, t_max)

# number of boxes
n_box = 100
# advection velocity
v_a = .014 * (100/10)
# diffusion coeff
k.h.dt = .01 * (100/10)^2
# initial box number
i_box = 50
# initial concentration
C_i = 10
# time series
t_max = 100
# produce plots
Conc = circular_boxes(n_box, v_a, k.h.dt, i_box, C_i, t_max)
