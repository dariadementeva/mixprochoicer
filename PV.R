#' Create a random initial design with PV
#'
#' \code{MPV_random_design} returns an initial random design for mixtures
#' and process variables
#'
#' @param q number of ingredients in a mixture
#' @param n number of runs required
#' @param z the number of process variables
#' @param m the number of process variables
#' @param seed set seed for reproducibility
#' @examples
#' MPV_random_design(3,30,3)
#' @export
#'

# In a Cartesian join, a mixture design is assigned to every treatment combination
# of a process-variable design, or equivalently,
# a process-variable design is assigned to every mixture blend of a mixture design.
# MPV = Mdesign x PVdesign
# Mario: pv bounds should be enclosed in [-1,1] interval

MPV_random_design <- function(q, n, m, m_bounds = c(-1, 1), seed = NULL) {

  if(!is.null(seed)) set.seed(seed)

  random_matrix = matrix(rep(NA_real_, q*n), nrow = n)

  for(i in 1:nrow(random_matrix)){
    random_values = runif(q)
    random_matrix[i,] = random_values / sum(random_values)
  }
  if(m > 0) {
    pv_matrix = matrix(rep(NA_real_, n*m), ncol = m)

    for(p in 1:m){
      pv_matrix[, p] = runif(n, min = m_bounds[1], max = m_bounds[2])
    }
    random_design = cbind(random_matrix, pv_matrix)
     }
  return(random_design)
  print(random_design)
  }

#' Computation of design matrix MPV model
#'
#'\code{MPV_model_matrix} returns a model matrix for the model expansion
#' of Kowalski et al., 2000, as cited in Goos & Jones, 2011
#' @param design the MPV design
#' @param q number of ingredients
#' @param m the number of process variables
#' @examples
#' MPV_model_matrix(design, 3,3)
#' @export
#'

MPV_model_matrix <- function(design, q, m){

  m_matrix = design[,1:q]
  pv_matrix = design[,(q+1):ncol(design)]
  n = nrow(design)

  mixture_parameters = q + q*(q-1)/2
  pv_parameter_int = m*(m-1)/2
  pv_parameter_quad = m
  mpv_parameters = q*m
  # total_parameters = mixture_parameters +  mpv_parameters + pv_parameter_int + pv_parameter_quad

  mixture_matrix= matrix(rep(NA_real_, mixture_parameters*n), nrow = n)
  mixture_matrix[,1:q] = m_matrix
  p = q
  for(k in 1:(q-1)){
    for(l in (k+1):q){
      p = p+1
      mixture_matrix[,p] = m_matrix[,k]*m_matrix[,l]
    }
  }

  mpv_matrix = matrix(rep(NA_real_, mpv_parameters*n), nrow = n)
  s = 0
  for (i in 1:m) {
    for (k in 1:q){
      s = s+1
      mpv_matrix[,s] = pv_matrix[,i]*m_matrix[,k]
    }
  }

  pv_matrix_int = matrix(rep(NA_real_, pv_parameter_int*n), nrow = n)
  w = 0
  for (i in 1:(m-1)) {
    for (j in (i+1):m){
      w = w+1
      pv_matrix_int[,w] =  pv_matrix[,i]*pv_matrix[,j]
    }
  }

  pv_matrix_quad = matrix(rep(NA_real_, pv_parameter_quad*n), nrow = n)
  r = 0
  for (i in 1:m) {
    r = r+1
    pv_matrix_quad[,r] =  pv_matrix[,i]*pv_matrix[,i]
  }

  model_matrix = cbind(mixture_matrix, mpv_matrix, pv_matrix_int, pv_matrix_quad)
  return(model_matrix)
}


#' Coordinate-exchange for MPV design
#'
#'\code{MPV_compute_d_optimality} returns a model matrix for the model expansion


MPV_compute_d_optimality = function(design, q, m){
  printed_model_matrix = MPV_model_matrix(design, q, m)
  transposed_mm = t(printed_model_matrix)
  d_opt = -as.numeric(determinant(transposed_mm%*%printed_model_matrix)$modulus)

  # if(d_opt < 1e-16){
  #   stop("Matrix is not invertible")
  # }

  return(d_opt)
}


#' Coordinate-Exchange for MPV design based on D-optimality
#'
#'\code{MPV_coordinate_exchange_d_optimal} returns a model matrix for the model expansion
#'
#'
#'
#

MPV_coordinate_exchange_d_optimal <- function(design, q, m, n_points, n_iter = 10) {

  d_optimal_design = design

  d_value_opt = MPV_compute_d_optimality(d_optimal_design, q, m)
  n = nrow(design)

  iter = 0
  while(iter < n_iter){
    iter = iter + 1

    for (row in 1:n) {

      for (ing in 1:q) {

        # cat("row: ", row, "ing: ", ing, "\n")
        cox_dir = compute_cox_direction(d_optimal_design[row,1:q], ing, n_points)
        for (entry in 1:nrow(cox_dir)) {
          new_design = d_optimal_design
          # entry=1
          new_design[row, 1:q] = cox_dir[entry,]
          d_value_opt = MPV_compute_d_optimality(d_optimal_design, q, m)
          d_value_new = MPV_compute_d_optimality(new_design, q, m)

          if (d_value_opt >= d_value_new){
            d_optimal_design = new_design
            d_value_opt = d_value_new
          }

        }

        for (z in 1:m){
          points_to_evaluate = seq(-1, 1, length.out = n_points)
          for(point in points_to_evaluate) {
            new_design = d_optimal_design
            # entry=1
            new_design[row, q+z] = point

            d_value_opt = MPV_compute_d_optimality(d_optimal_design, q, m)
            d_value_new = MPV_compute_d_optimality(new_design, q, m)

            if (d_value_opt >= d_value_new){
              d_optimal_design = new_design
              d_value_opt = d_value_new
            }

          }
        }

      }
    }
  }
  return(d_optimal_design)
}

#' Moments Matrix for I-optimality
#'
#'\code{MPV_coordinate_exchange_d_optimal} returns a model matrix for the model expansion
#'
#'
#'
#

MPV_create_moment_matrix = function(q, m, order, m_bounds=c(-1,1)) {

  stopifnot(order %in% 1:4)

  if (order == 1) {
    parameters = q
  } else {
    if (order == 2) {
      parameters = q + q * (q - 1)/2
    } else {
      if (order == 3) {
        parameters = q + q*(q - 1)/2 + q * (q - 1) * (q - 2)/6
      }
      else {
        parameters = q + q*(q-1)/2 + q*m + m*(m-1)/2 + m
      }
    }
  }

  auxiliary_matrix_f = matrix(rep(0L, parameters*(q + m)), ncol = q + m)

  counter = 0
  for(i in 1:q){
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

  MPV_moments_matrix = matrix(rep(NA_real_, parameters^2), ncol = parameters)

  for(i in 1:(parameters)){
    for(j in 1:(parameters)){

      auxiliary_vector = auxiliary_matrix_f[i, 1:q] + auxiliary_matrix_f[j, 1:q]
      numerator = prod(factorial(auxiliary_vector))
      denominator = factorial(q + sum(auxiliary_vector))

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
      denominator_2 = prod(1 + auxiliary_vector2)
      MPV_moments_matrix[i,j] = (numerator * numerator_2)/(denominator * denominator_2)
    }
  }

  return(MPV_moments_matrix)
}


#' Compute I-optimality value for the pre-specified design and order
#'
#' \code{MNL_MPV_compute_i_optimality} computes I-optimality value for the
#' pre-specified design and order
#'
#' @param design specifies the design matrix
#' @param order specifies the order of the model (1, 2 or 3)
#' @param q specifies the number of ingredients
#' @param m specifies the number of process variables
#'
#'
#'
#' @export

MPV_compute_i_optimality = function(design, q, m, order=4, m_bounds=c(1,1)){

  printed_model_matrix = MPV_model_matrix(design, q, m)
  printed_moments_matrix = MPV_create_moment_matrix(q, m, order, m_bounds=c(-1,1))

  i_opt = tr(solve(t(printed_model_matrix)%*%printed_model_matrix)%*%printed_moments_matrix)
  return(i_opt)
}


#' Coordinate-Exchange for MPV design based on I-optimality
#'
#'\code{MPV_coordinate_exchange_d_optimal} returns a model matrix for the model expansion
#'
#'
#'
#

MPV_coordinate_exchange_i_optimal <- function(design, q, m, n_points, n_iter = 10) {

  i_optimal_design = design

  i_value_opt = MPV_compute_i_optimality(i_optimal_design, q, m)
  n = nrow(design)

  iter = 0
  while(iter < n_iter){
    iter = iter + 1

    for (row in 1:n) {

      for (ing in 1:q) {

        # cat("row: ", row, "ing: ", ing, "\n")
        cox_dir = compute_cox_direction(i_optimal_design[row,1:q], ing, n_points)
        for (entry in 1:nrow(cox_dir)) {
          new_design = i_optimal_design
          # entry=1
          new_design[row, 1:q] = cox_dir[entry,]
          i_value_opt = MPV_compute_i_optimality(i_optimal_design, q, m)
          i_value_new = MPV_compute_i_optimality(new_design, q, m)

          if (i_value_opt >= i_value_new){
            i_optimal_design = new_design
            i_value_opt = i_value_new
          }

        }

        for (z in 1:m){
          points_to_evaluate = seq(-1, 1, length.out = n_points)
          for(point in points_to_evaluate) {
            new_design = i_optimal_design
            # entry=1
            new_design[row, q+z] = point

            i_value_opt = MPV_compute_i_optimality(i_optimal_design, q, m)
            i_value_new = MPV_compute_i_optimality(new_design, q, m)

            if (i_value_opt >= i_value_new){
              i_optimal_design = new_design
              i_value_opt = i_value_new
            }

          }
        }

      }
    }
  }
  return(i_optimal_design)
}





#' Computation of design matrix for Scheffe Model with PV: Order 1
#'
#' \code{PV_scheffe_model_order1} returns the model matrix for Scheffe Model (1)
#'
#' @param design specifies the model matrix
#' @examples
#'
#'
#' @export
#'


# coordinate-exchange for MPV design
# Model matrix: Equation 6.12, p. 127
# Data Part: Transportation Example - ask Peter OR Theoretical Example
#'
#'
#' PV_scheffe_model_order1 <- function(design){
#'   m_matrix = design[,1:q]
#'   pv_matrix = design[,(q+1):ncol(design)]
#'   n = nrow(design)
#'   q = ncol(m_matrix)
#'   z = ncol(pv_matrix)
#'
#'   parameters = q*(1+ z + (z*(z-1)/2))
#'
#'   model_matrix = matrix(rep(NA_real_, parameters*n), nrow = n)
#'   model_matrix[,1:q] = m_matrix
#'   k = q
#'   for (j in 1:q){
#'     for (i in 1:z){
#'       k = k+1
#'       model_matrix[,k] = m_matrix[,j]*pv_matrix[,i]
#'     }
#'
#'
#'   }
#' return(model_matrix)
#' }
#'
#' #' Computation of model matrix for Scheffe Model with PV: Order 2
#' #'
#' #' \code{PV_scheffe_model_order2} returns the model matrix for Scheffe Model (2)
#' #'
#' #' @param design specifies the design matrix
#' #' @examples
#' #'
#' #'
#' #' @export
#'
#' PV_scheffe_model_order2 <- function(design){
#'   m_matrix = design[,1:q]
#'   pv_matrix = design[,(q+1):ncol(design)]
#'   n = nrow(design)
#'   q = ncol(m_matrix)
#'   z = ncol(pv_matrix)
#'   parameters = ((q + q*(q-1)/2)) * (1+ z + (z*(z-1)/2))
#'
#'   model_matrix = matrix(rep(NA_real_, parameters*n), nrow = n)
#'   model_matrix[,1:q] = m_matrix
#'
#'   # k = q
#'   # for(i in 1:(q-1)){
#'   #   for(j in (i+1):q){
#'   #     k = k+1
#'   #     model_matrix[,k] = design[,i]*design[,j]
#'   #   }
#'   # }
#'   # return(model_matrix)
#' }
#'
#' #' Computation of model matrix for Scheffe Model with PV: Order 3. Built on SM-O2
#' #'
#' #' \code{PV_scheffe_model_order3} returns the model matrix for Scheffe Model (3)
#' #'
#' #' @param design specifies the design matrix
#' #' @examples
#' #' scheffe_model_order3(my_random_design)
#' #'
#' #'@export
#'
#' PV_scheffe_model_order3 <- function(design) {
#'   m_matrix = design[,1:q]
#'   pv_matrix = design[,(q+1):ncol(design)]
#'   n = nrow(design)
#'   q = ncol(m_matrix)
#'   z = ncol(pv_matrix)
#'
#'   parameters = ((q * (q - 1) * (q - 2)/6) * (1+ z + (z*(z-1)/2)))
#'
#'   # for(i in 1:(q-2)){
#'   #   for(j in (i+1):(q-1)){
#'   #     for(k in (j+1):q){
#'   #       # cat("i = ", i, ", j = ", j, "\n") # uncomment if needs to be printed
#'   #       model_matrix = cbind(model_matrix, design[,i]*design[,j]*design[,k])
#'   #     }
#'   #   }
#'   # }
#'   return(model_matrix)
#' }
