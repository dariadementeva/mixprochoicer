#' Create a random initial design with PV for MNL
#'
#' \code{MNL_MPV_random_design} returns an initial random design for mixtures
#' and process variables when MNL model is implied
#'
#' @export
#'

MNL_MPV_random_design = function(q, J, S, m = 0,  m_bounds = c(-1,1), seed = NULL){

   if(!is.null(seed)) set.seed(seed)

  random_design = array(rep(NA_real_, (q+m)*J*S), dim = c((q+m), J, S))

  for(j in 1:J){
    for(s in 1:S){
      rands = runif(q)
      random_design[1:q,j, s] = rands/sum(rands)
    }
  }
  if(m > 0){
    for(z in 1:m){
      random_design[z + q, , ] = runif(n = J*S, min = m_bounds[1], max = m_bounds[2])
    }
  }
  return(random_design)
}

#' Create a moments matrix for I-optimality for MNL-MPV model
#'
#' \code{MNL_MPV_create_moment_matrix} returns a moments matrix for I-optimality
#' in the MNL context framework
#'
#' @export
#'

MNL_MPV_create_moment_matrix = function(q, m, order, m_bounds=c(-1,1)) {

  stopifnot(order %in% 1:4)

  if (order == 1) {
    parameters = q - 1
  } else {
    if (order == 2) {
      parameters = q - 1 + q * (q - 1)/2
    } else {
      if (order == 3) {
        parameters = q - 1 + q*(q - 1)/2 + q * (q - 1) * (q - 2)/6
      }
      else {
        parameters = q - 1 + q*(q-1)/2 + q*m + m*(m-1)/2 + m
      }
    }
  }

  auxiliary_matrix_f = matrix(rep(0L, parameters*(q + m)), ncol = q + m)

  counter = 0
  for(i in 1:(q-1)){
    counter = counter + 1
    auxiliary_matrix_f[counter, i] = 1
  }

  if(order >= 2){
    for(i in 1:(q-1)){
      for(j in (i+1):q){
        counter = counter + 1
        auxiliary_matrix_f[counter, i] = 1
        auxiliary_matrix_f[counter, j] = 1
      }
    }
  }

  if(order == 3){
    for(i in 1:(q-2)){
      for(j in (i+1):(q-1)){
        for(k in (j+1):q){
          counter = counter + 1
          auxiliary_matrix_f[counter, i] = 1
          auxiliary_matrix_f[counter, j] = 1
          auxiliary_matrix_f[counter, k] = 1
        }
      }
    }
  }
  if(order == 4){

    for(i in 1:m){
      for(k in 1:q){
        counter = counter + 1
        auxiliary_matrix_f[counter, q + i] = 1
        auxiliary_matrix_f[counter, k] = 1
      }
    }

    if(m > 1){
      for(i in 1:(m-1)){
        for(j in (i+1):m){
          counter = counter + 1
          auxiliary_matrix_f[counter, q + i] = 1
          auxiliary_matrix_f[counter, q + j] = 1
        }
      }
    }

    for(i in 1:m){
      counter = counter + 1
      auxiliary_matrix_f[counter, q + i] = 2
    }
  }

  MNL_MPV_moments_matrix = matrix(rep(NA_real_, parameters^2), ncol = parameters)

  for(i in 1:(parameters)){
    for(j in 1:(parameters)){

      auxiliary_vector = auxiliary_matrix_f[i, 1:q] + auxiliary_matrix_f[j, 1:q]
      numerator = prod(factorial(auxiliary_vector))
      denominator = factorial(q - 1 + sum(auxiliary_vector))

      if(m > 0){
        auxiliary_vector2 = auxiliary_matrix_f[i, (q+1):(q+m)] + auxiliary_matrix_f[j, (q+1):(q+m)]
        numerator_2 = 1
        for(l in 1:m){
          numerator_2 = numerator_2 * (m_bounds[2]^(auxiliary_vector2[l] + 1) - m_bounds[1]^(auxiliary_vector2[l] + 1))
        }

      } else{
        auxiliary_vector2 = 0
        numerator_2 = 1
      }
      denominator_2 = prod(1 +auxiliary_vector2)
      MNL_MPV_moments_matrix[i,j] = (numerator * numerator_2)/(denominator * denominator_2)
    }
  }

  return(MNL_MPV_moments_matrix)
}

#' Create a random initial design with PV for MNL
#'
#' \code{MNL_MPV_model_matrix_order4} returns model expansion for order 4 for MNL with PV
#'
#' @export
#'

# MNL_MPV_model_matrix_order4 = function(design, m) {
#   #
#   q = dim(design)[1] - m
#  # m = dim(design)[1] - q
#   J = dim(design)[2]
#   S = dim(design)[3]
#
#   m_array_MNL = design[1:(q-1), ,]
#   m_array_RS = design[1:q, ,]
#   pv_array = design[(q+1):m, ,]
#   n = nrow(design)
#
#   mixture_parameters = q - 1 + q*(q-1)/2
#   pv_parameter_int = m*(m-1)/2
#   pv_parameter_quad = m
#   mpv_parameters = q*m
#
#   mixture_array= array(rep(NA_real_, mixture_parameters*J*S), dim = c(mixture_parameters, J, S))
#   mixture_array[1:(q-1), , ] =  m_array_MNL
#   counter = q-1
#   for(i in 1:(q-1)){
#     for(k in (i+1):q){
#       counter = counter + 1
#       mixture_array[counter, , ] = m_array_RS[i, , ]*m_array_RS[k, , ]
#     }
#   }
#
#   mpv_array = array(rep(NA_real_, mpv_parameters*J*S), dim = c(mpv_parameters, J, S))
#   s = 0
#   for (i in 1:m) {
#     for (k in 1:q){
#       s = s+1
#       mpv_array[s, , ] = pv_array[i, , ]*m_array_RS[k, , ]
#     }
#   }
#
#   pv_array_int =  array(rep(NA_real_, pv_parameter_int*J*S), dim = c(pv_parameter_int, J, S))
#   w = 0
#   for (i in 1:(m-1)) {
#     for (j in (i+1):m){
#       w = w+1
#       pv_array_int[w, , ] = pv_array[i, , ]*pv_array[j, , ]
#     }
#   }
#
#   pv_array_quad = array(rep(NA_real_, pv_parameter_quad*J*S), dim = c(pv_parameter_quad, J, S))
#   r = 0
#   for (i in 1:m) {
#     r = r+1
#     pv_array_quad[r, , ] =  pv_array[i, , ]*pv_array[i, , ]
#   }
#
#   model_array = abind(mixture_array, mpv_array, pv_array_int, pv_array_quad, along=1)
#   return(model_array)
# }


# DISCUSS

#' Compute a model matrix for MPV MNL case
#'
#' \code{MNL_MPV_compute_model_matrix} returns model expansion for order 1,2,3,4 for MNL-MPV case
#'
#' @export
#'

# MNL_MPV_compute_model_matrix = function(design, m, order){
#   if(order == 4){
#     MNL_MPV_model_matrix_order4(design, m)
#   } else{
#     MNL_compute_model_matrix(design, order)
#   }
# }


MNL_MPV_compute_model_matrix <- function(design, m, order) {
  q = dim(design)[1] - m
  J = dim(design)[2]
  S = dim(design)[3]

  p1 = q - 1
  p2 = q * (q - 1)/2
  p3 = q * (q - 1) * (q - 2)/6
  p4 = q*m + m*(m-1)/2 + m


  if (order == 1) {
    n_parameters = p1
  } else {
    if (order == 2) {
      n_parameters = p1 + p2
    } else {
      if (order == 3) {
        n_parameters = p1 + p2 + p3
      } else {
        if (order == 4) {
          n_parameters = p1 + p2 + p4
        }
      }
    }
  }

  model_array = array(rep(NA_real_, n_parameters * J * S), dim = c(n_parameters, J, S))

  if (order >= 1) {

    model_array[1:(q - 1), , ] = design[1:(q - 1), , ]
  }

  if (order >= 2) {
    counter = q - 1
    for (i in 1:(q - 1)) {
      for (j in (i + 1):q) {
        counter = counter + 1
        model_array[counter, , ] = design[i, , ] * design[j, , ]
      }
    }
  }
  if (order == 3) {
    # counter = 0
    for (i in 1:(q - 2)) {
      for (j in (i + 1):(q - 1)) {
        for (k in (j + 1):q) {
          counter = counter + 1
          model_array[counter, , ] = design[i, , ] * design[j, , ] * design[k, , ]
        }
      }
    }
  }
  if (order == 4) {
    #counter = p1 + p2
    for (z in 1:m) {
      for (p in 1:q){
        counter = counter + 1
        model_array[counter, , ] = design[p, , ] * design[z+q, , ]
      }
    }

    for (i in 1:(m-1)) {
      for (j in (i+1):m){
        counter = counter + 1
        model_array[counter, , ] = design[i+q, , ]*design[j+q, , ]
      }
    }

    for (i in 1:m) {
      counter = counter + 1
      model_array[counter, , ] =  design[i+q, , ]*design[i+q, , ]
    }
  }
  return(model_array)
}

#' Compute a model matrix for MPV MNL case
#'
#' \code{MNL_MPV_compute_choice_probabilites_order4} returns choice probabilities matrix
#' @export
#'

MNL_MPV_compute_choice_probabilites4 = function(design, m, beta, order=4) {
  q = dim(design)[1] - m
  J = dim(design)[2]
  S = dim(design)[3]

  U = matrix(rep(NA_real_, J * S), ncol = S)

  model_matrix = MNL_MPV_compute_model_matrix(design, m, order=4)

  S = dim(model_matrix)[3]

  for(s in 1:S){
    Xs = model_matrix[, , s]
    Us = t(Xs) %*% matrix(beta, ncol = 1)
    U[, s] = Us
  }
  # define choice probability of alternative j in choice s
  P = matrix(rep(NA_real_, J * S), ncol = S)
  for (j in 1:J) {
    for (s in 1:S) {
      u_js = U[j, s]
      p_js = exp(u_js)/sum(exp(U[, s]))
      P[j, s] = p_js  # choice probabilities form matrix of model expansion
    }
  }
  return(P)
}


#' Compute a model matrix for MPV MNL case
#'
#' \code{MNL_MPV_compute_choice_probabilites} returns choice probabilities matrix
#' @export
#'

# MNL_MPV_compute_choice_probabilites = function(design, m, beta, order){
#   if(order == 4){
#     MNL_MPV_compute_choice_probabilites_order4(design, m, beta)
#   } else{
#     MNL_compute_choice_probabilites(design, beta, order)
#   }
# }

#' Computeinformation matrix order for MPV MNL case
#'
#' \code{MNL_MPV_compute_information_matrix_order4} returns information matrix for order 4
#' @export
#'

MNL_MPV_compute_information_matrix4 = function(design, m, beta, order=4) {
  q = dim(design)[1] - m
  J = dim(design)[2]
  S = dim(design)[3]

  p1 = q - 1
  p2 = q * (q - 1)/2
  p3 = q * (q - 1) * (q - 2)/6
  p4 = q*m + m*(m-1)/2 + m

  n_parameters = p1 + p2 + p4

  model_array = MNL_MPV_compute_model_matrix(design, m, order=4)
  P = MNL_MPV_compute_choice_probabilites4(design, m, beta, order=4)

  information_matrix = matrix(rep(0, n_parameters * n_parameters), ncol = n_parameters)
  for (s in 1:S) {
    p_s = P[, s]
    I_s = model_array[1:n_parameters, , s] %*% (diag(p_s) - (t(t(p_s)) %*% t(p_s))) %*%
      t(model_array[1:n_parameters, , s])
    information_matrix = information_matrix + I_s
  }
  return(information_matrix)
}

#' Compute information matrix order for MPV MNL case, all orders
#'
#' \code{MNL_MPV_compute_information_matrix} returns information for orders 1,2,3,4
#' @export
#'

# MNL_MPV_compute_information_matrix = function(design, m, beta, order){
#   if(order == 4){
#     MNL_MPV_compute_information_matrix_order4(design, m, beta)
#   } else{
#     MNL_compute_information_matrix(design, beta, order)
#   }
# }

#' Compute D-optimality for MNL MPV case
#'
#' \code{MNL_MPV_compute_d_optimality} computes D-optimality for MNL MPV case for all orders
#'

MNL_MPV_compute_d_optimality = function(design, beta, m, order=4) {
  mnl_mpv_information_matrix = MNL_MPV_compute_information_matrix4(design, m, beta, order=4)
  log_determinant = determinant(t(mnl_mpv_information_matrix) %*% mnl_mpv_information_matrix)$modulus
  return(as.numeric(log_determinant))
}


#' Compute I-optimality value for MNL MPV case
#' \code{MNL_MPV_compute_i_optimality} computes I-optimality value for MNL-MPV case for all orders
#'
#' @export

MNL_MPV_compute_i_optimality = function(design, m, order, m_bounds=c(1,1)){
  q = dim(design)[1] - m
  printed_model_matrix = MNL_MPV_compute_model_matrix(design, m, order)
  printed_moments_matrix = MNL_MPV_create_moment_matrix(q, m, order, m_bounds=c(-1,1))

  i_opt = sum(diag(solve(printed_information_matrix, printed_moments_matrix)))
  return(i_opt)
}

# ASK MARIO ______________________________________________________________
#' Perform Coordinate-Exchange for D-optimality, MNL-MPV Case
#' \code{MNL_MPV_coordinate_exchange_d_optimal} performs Coordinate-Exchange for D-optimality
#'
#' @export

MNL_MPV_coordinate_exchange_d_optimal <- function(design, m, order, beta, n_points, n_iter = 10) {
  q = dim(design)[1] - m
  J = dim(design)[2]
  S = dim(design)[3]

  d_optimal_design = design

  d_value_opt = MNL_MPV_compute_d_optimality(design, beta, m, order)
  # n = nrow(design)

  iter = 0
  while (iter < n_iter) {
    iter = iter + 1

    for (j in 1:J) {
      for (s in 1:S) {
        for (ing in 1:q) {
          cox_dir = compute_cox_direction(d_optimal_design[, j, s], ing,
            n_points)
          for (entry in 1:nrow(cox_dir)) {
            new_design = d_optimal_design
            # entry = 1
            new_design[, j, s] = cox_dir[entry, ]
            d_value_opt = MNL_MPV_compute_d_optimality(d_optimal_design, beta, m, order)
            d_value_new = MNL_MPV_compute_d_optimality(new_design, beta, m, order)

            if (d_value_opt >= d_value_new) {
              d_optimal_design = new_design
              d_value_opt = d_value_new
            }
          }

          if(m > 0){
            for (z in 1:m){
              points_to_evaluate = seq(-1, 1, length.out = n_points)
              for(point in points_to_evaluate) {
                new_design = d_optimal_design
                # entry=1
                new_design[q+z, j, s] = point #ASK!

                d_value_opt = MNL_MPV_compute_d_optimality(d_optimal_design, beta, m, order)
                d_value_new = MNL_MPV_compute_d_optimality(new_design, beta, m, order)

                if (d_value_opt >= d_value_new){
                  d_optimal_design = new_design
                  d_value_opt = d_value_new
                }
              }
            }
          } # end if

        }
      }
    }
  }
  return(d_optimal_design)
}


# ASK MARIO ______________________________________________________________
#' Perform Coordinate-Exchange for I-optimality, MNL-MPV Case
#' \code{MNL_MPV_coordinate_exchange_i_optimal} performs Coordinate-Exchange for I-optimality
#'
#' @export

MNL_MPV_coordinate_exchange_i_optimal <- function(design, m, beta, n_points, n_iter = 10) {
  q = dim(design)[1] - m
  J = dim(design)[2]
  S = dim(design)[3]

  i_optimal_design = design

  i_value_opt = MNL_MPV_compute_i_optimality(design, m, order, m_bounds=c(1,1))
  n = nrow(design)

  iter = 0
  while (iter < n_iter) {
    iter = iter + 1

    for (j in 1:J) {
      for (s in 1:S) {
        for (ing in 1:q) {
          cox_dir = compute_cox_direction(i_optimal_design[, j, s], ing,
            n_points)
          for (entry in 1:nrow(cox_dir)) {
            new_design = i_optimal_design
            # entry = 1
            new_design[, j, s] = cox_dir[entry, ]
            i_value_opt = MNL_MPV_compute_i_optimality(i_optimal_design, m, order, m_bounds=c(1,1))
            i_value_new = MNL_MPV_compute_i_optimality(new_design,m, order, m_bounds=c(1,1))

            if (i_value_opt >= i_value_new) {
              i_optimal_design = new_design
              i_value_opt = i_value_new
            }
          }

          if(m > 0){
            for (z in 1:m){
              points_to_evaluate = seq(-1, 1, length.out = n_points)
              for(point in points_to_evaluate) {
                new_design = i_optimal_design
                # entry=1
                new_design[q+z, j, s] = point #ASK!

                i_value_opt = MNL_MPV_compute_i_optimality(i_optimal_design, beta, m, order)
                i_value_new = MNL_MPV_compute_i_optimality(new_design, beta, m, order)

                if (i_value_opt <= i_value_new){
                  i_optimal_design = new_design
                  i_value_opt = i_value_new
                }
              }
            }
          } # end if
        }
      }
    }
  }
  return(i_optimal_design)
}




