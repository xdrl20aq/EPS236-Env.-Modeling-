setwd("~/projects/EPS236-Env.-Modeling-/")

# function to create eigenvalue diagonal matrix Lambda
exp_Lambda = function(lambda_values, t) {
  return( diag(exp(lambda_values * t)) )
}

# initial conditions
# number of layers
n_layer = 4
boxes = rep(NA, n_layer)
for (i in 1:n_layer) {
  if (i == 1) {
    boxes[i] = 1
  } else {
    edge = i * 2 - 1
    boxes[i] = (edge - 1) * 4
  }
}
boxes
n_box = sum(boxes)

# advection velocity
v_a = .014
# diffusion coeff
kh1 = .05
kh2 = .03
kh3 = .03
# initial box number
i_box = 10
# initial concentration
C_i = 10
# time series
t = 0:100

# flows
corner1 = c(kh3, kh1, kh1, kh3)
edge1 = c(kh3, kh1, kh2, kh1)
corner2 = c(kh3, kh3, kh1, kh1)
edge2 = c(kh1, kh3, kh1, kh2)
corner3 = c(kh1, kh3, kh3, kh1)
edge3 = c(kh2, kh1, kh3, kh1)
corner4 = c(kh1, kh1, kh3, kh3)
edge4 = c(kh1, kh2, kh1, kh3)


K = matrix(0.0, nrow = n_box, ncol = n_box)

# diffusion
add_circular_diffusion_term = function(K, kh1, kh2, kh3, n_layer, boxes) {
  for (layer_i in 1:n_layer) {
    if (layer_i == 1) {
      K[1, c(3, 5, 7, 9)] = rep(kh3, 4)
      K[1, 1] = kh3 * -4
    } else if (layer_i == 2) {
      K[2, c(11, 3, 9, 25)] = corner1
      K[2, 2] = K[2, 2] - sum(corner1)
      K[3, c(12, 4, 1, 2)] = edge1
      K[3, 3] = K[3, 3] - sum(edge1)
      K[4, c(13, 15, 5, 3)] = corner2
      K[4, 4] = K[4, 4] - sum(corner2)
      K[5, c(4, 16, 6, 1)] = edge2
      K[5, 5] = K[5, 5] - sum(edge2)
      K[6, c(5, 17, 19, 7)] = corner3
      K[6, 6] = K[6, 6] - sum(corner3)
      K[7, c(1, 6, 20, 8)] = edge3
      K[7, 7] = K[7, 7] - sum(edge3)
      K[8, c(9, 7, 21, 23)] = corner4
      K[8, 8] = K[8, 8] - sum(corner4)
      K[9, c(2, 1, 8, 24)] = edge4
      K[9, 9] = K[9, 9] - sum(edge4)
    } else {
      outer = layer_i * 2 - 1
      inner = (layer_i - 1) * 2 - 1
      startbox = inner^2+1
      endbox = outer^2
      for (i in startbox:endbox) {
        flow = rep(NA, 4)
        top = 0; right = 0; bottom = 0; left = 0
        cornerNum = (i - startbox) / ((layer_i - 1) * 2)
        if (i == startbox) {
          top = i + (9 + 8 * (layer_i - 2))
          right = i + 1
          bottom = endbox
          left = i + (39 + 16 * (layer_i - 3))
          flow = corner1
        } else if (0 < cornerNum & cornerNum < 1) {
          top = i + (9 + 8 * (layer_i - 2))
          right = i + 1
          bottom = i - (9 + 8 * (layer_i - 3))
          left = i - 1
          flow = edge1
        } else if (cornerNum == 1) {
          top = i + (9 + 8 * (layer_i - 2))
          right = i + (9 + 8 * (layer_i - 2) + 2)
          bottom = i + 1
          left = i - 1
          flow = corner2
        } else if (1 < cornerNum & cornerNum < 2) {
          top = i - 1
          right = i + (9 + 8 * (layer_i - 2) + 2)
          bottom = i + 1
          left = i - (9 + 8 * (layer_i - 3) + 2)
          flow = edge2
        } else if (cornerNum == 2) {
          top = i - 1
          right = i + (9 + 8 * (layer_i - 2) + 2)
          bottom = i + (9 + 8 * (layer_i - 2) + 4)
          left = i + 1
          flow = corner3
        } else if (2 < cornerNum & cornerNum < 3) {
          top = i - (9 + 8 * (layer_i - 3) + 4)
          right = i - 1
          bottom = i + (9 + 8 * (layer_i - 2) + 4)
          left = i + 1
          flow = edge3
        } else if (cornerNum == 3) {
          top = i + 1
          right = i - 1
          bottom = i + (9 + 8 * (layer_i - 2) + 4)
          left = i + (9 + 8 * (layer_i - 2) + 6)
          flow = corner4
        } else if (i == endbox) {
          top = startbox
          right = i - (23 + 16 * (layer_i - 3))
          bottom = i - 1
          left = i + (9 + 8 * (layer_i - 2) + 6)
          flow = edge4
        } else if (3 < cornerNum & cornerNum < 4) {
          top = i + 1
          right = i - (9 + 8 * (layer_i - 3) + 6)
          bottom = i - 1
          left = i + (9 + 8 * (layer_i - 2) + 6)
          flow = edge4
        }
        adjacent = c(top, right, bottom, left)
        flow_fix = flow[adjacent <= n_box]
        adjacent_fix = adjacent[adjacent <= n_box]
        K[i, adjacent_fix] = K[i, adjacent_fix] + flow_fix
        K[i, i] =　K[i, i] - sum(flow_fix)
      }
    }
  }
  return( K )
}


# advection - call if v_a != 0
add_circular_advection_term = function(K, v_a, n_box, n_layer) {
  for (layer_i in 2:n_layer) {
    outer = layer_i * 2 - 1
    inner = (layer_i - 1) * 2 - 1
    startbox = inner^2+1
    endbox = outer^2
    for (i in startbox:endbox) {
      K[i, i] = K[i, i] - v_a
      if (v_a > 0) {
        if (i == startbox) {
          K[i, endbox] = K[i, endbox] + v_a
        } else {
          K[i, i-1] = K[i, i-1] + v_a
        }
      } else {
        if (i == endbox) {
          K[i, startbox] = K[i, startbox] + v_a
        } else {
          K[i, i+1] = K[i, i+1] + v_a
        }
      }
    }
  }
  return( K )
}

# construct transition matrix
K = add_circular_diffusion_term(K, kh1, kh2, kh3, n_layer, n_box)

# making sure mass is conserved
sum(K)
which(abs(rowSums(K)) > 0.001)
which(abs(colSums(K)) > 0.001)

K = add_circular_advection_term(K, v_a, n_box, n_layer)

sum(K)
which(abs(rowSums(K)) > 0.01)

library(Matrix)
K_0 = bdiag(K, K)

add_connection_diffusion = function(K, kh1, n_layer) {
  nexttop = (n_layer * 2 - 1)^2 + 1
  K[1, nexttop] = K[1, nexttop] + kh1
  K[1, 1] = K[1, 1] - kh1
  K[nexttop, 1] = K[nexttop, 1] + kh1
  K[nexttop, nexttop] = K[nexttop, nexttop] - kh1
  return( K )
}

K_0 = add_connection_diffusion(K_0, kh1, n_layer)

# eigenvalues and eigenvectors
eigens = eigen(K_0)
e_vals = eigens$values
Xe = eigens$vectors

# initial condition
C_0 = matrix(0, nrow = n_box * 2, ncol = 1)
C_0[i_box] = C_i


# matrix to store the transition of concentration
Conc = matrix(nrow = n_box * 2, ncol = length(t))
Conc[, 1] = C_0
# loop over time to calculate C(t) for t
for (t_i in t) {
  if (t_i != 0) {
    Lambda_t = exp_Lambda(e_vals, t_i)
    conc_t_i = Re(Xe %*% Lambda_t %*% solve(Xe, C_0))
    Conc[, t_i+1] = conc_t_i
    # plotting!
    filename = formatC(t_i, width = 3, format = "d", flag = "0")
    png(paste("./part2_movie/", filename, ".png", sep=""))
    plot(x = 1:(n_box*2), y = conc_t_i, xlim = c(1, n_box*2), ylim = c(0, C_i), xlab = "#box", ylab = "concentration")
    dev.off()
  }
}

