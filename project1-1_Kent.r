setwd("~/projects/EPS236-Env.-Modeling-/")

# function to create eigenvalue diagonal matrix Lambda
exp_Lambda = function(lambda_values, t) {
  return( diag(exp(lambda_values * t)) )
}

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
t = 0:100

K = matrix(0, nrow = n_box, ncol = n_box)
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

# construct transition matrix
K = add_advection_term(K, v_a, n_box)
K = add_diffusion_term(K, k.h.dt, n_box)

# eigenvalues and eigenvectors
eigens = eigen(K)
e_vals = eigens$values
Xe = eigens$vectors
Xe.inv = solve(Xe)

# initial condition
C_0 = matrix(0, n_box, 1)
C_0[i_box] = C_i

# matrix to store the transition of concentration
Conc = matrix(nrow = n_box, ncol = length(t))

# loop over time to calculate C(t) for t
for (t_i in t) {
  Lambda = exp_Lambda(e_vals, t_i)
  conc_t_i = Re(Xe %*% Lambda %*% Xe.inv %*% C_0)
  Conc[, t_i+1] = conc_t_i
  png(paste("./part1_movie/", t_i, ".png", sep=""))
  plot(x = 1:n_box, y = conc_t_i, xlim = c(1, n_box), ylim = c(0, C_i), xlab = "#box", ylab = "concentration")
  dev.off()
}

advection_diffusion_eqn = function(D, x, t_i) {
  exp_term = exp( -(x - v_a*t_i)^2 / (4*D*t_i) )
  return(1 / sqrt(4*pi*D*t_i) * exp_term)
}

t_1 = 1:100
box_num = 26
plot(t, Conc[box_num,], ylim = c(0, C_i))
d_list = seq(0.01, 0.05, 0.001)
sqrt_err_list = rep(NA, length(d_list))
idx = 1
for (d_guess in d_list) {
  diff = Conc[box_num,] - advection_diffusion_eqn(d_guess, abs(box_num - i_box), t_1)
  sqrt_err = sqrt(sum(diff^2))
  sqrt_err_list[idx] = sqrt_err
  idx = idx + 1
}
min_sqrt_err_idx = which(sqrt_err_list == min(sqrt_err_list))
d_best = d_list[min_sqrt_err_idx]
d_best
lines(t_1, advection_diffusion_eqn(d_best, abs(box_num - i_box), t_1))
