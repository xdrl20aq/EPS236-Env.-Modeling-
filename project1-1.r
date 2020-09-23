# function to create eigenvalue diagonal matrix Lambda
exp_Lambda = function(lambda_values, t) {
  return( diag(exp(lambda_values * t)) )
}

# initial conditions
# number of boxes
n_box = 10
# advection velocity
v_a = .5
# diffusion coeff
k.h.dt = .1
# initial box number
i_box = 5
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
      K[i, (10+i-1)%%n_box] = K[i, (10+i-1)%%n_box] + k.h.dt
    }
    if (i == 9) {
      K[i, n_box] = K[i, n_box] + k.h.dt
    }
    K[i, (10+i+1)%%n_box] = K[i, (10+i+1)%%n_box] + k.h.dt
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
        K[i, (10+i-1)%%n_box] = K[i, (10+i-1)%%n_box] + v_a
      }
    } else {
      if (i == 9) {
        K[i, n_box] = K[i, n_box] + v_a
      }
      K[i, (10+i+1)%%n_box] = K[i, (10+i+1)%%n_box] + v_a
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
C_0

# matrix to store the transition of concentration
Conc = matrix(nrow = n_box, ncol = length(t))

# loop over time to calculate C(t) for t
for (t_i in t) {
  Lambda = exp_Lambda(e_vals, t_i)
  conc_t_i = Xe %*% Lambda %*% Xe.inv %*% C_0
  Conc[, t_i] = conc_t_i
}
conc_sol = Re(Conc)

# plot box 1
plot(conc_sol[1, ])

