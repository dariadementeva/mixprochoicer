install.packages(c("devtools", "roxygen2", "testthat", "knitr"))
install.packages("rstudioapi")
install.packages("formatR")
formatR::tidy_dir("R")
rstudioapi::isAvailable("0.99.149")
library(devtools)
library(roxygen2)
library(testthat)
library(knit)
library(dplyr)
library(ggplot2)
library(ggtern)
library(psych)
has_devel()
devtools::load_all()
devtools::document()
roxygen2::roxygenise()

#' Create a random initial design
#'
#' \code{create_random_initial_design} returns an initial random design
#'
#' @param q number of ingredients in a mixture
#' @param n number of runs required
#' @param seed set seed for reproducibility
#' @examples
#' create_random_initial_design(4,12,1000)
#' create_random_initial_design(5,36)
#'
#' @export

create_random_initial_design <- function(q, n, seed = NULL) {
  if(!is.null(seed)) set.seed(seed)
  random_matrix = matrix(rep(NA_real_, q*n), nrow = n)
  for(i in 1:nrow(random_matrix)){
    random_values = runif(q)
    random_matrix[i,] = random_values / sum(random_values)
  }
  return(random_matrix)
  print(random_matrix)
}

#' Computation of Cox-effect direction (backbone function)
#'
#' \code{compute_cox_direction} computates Cox-effect direction for the each row
#'
#' @param s vector of the proportions of the q components (sums to 1)
#' @param i the  proportion ofcomponent i
#' @param n_points number of cox-effect points
#' @examples
#' compute_cox_direction(my_random_design[2,],2,30)
#'
#' @export

compute_cox_direction = function(s, i, n_points){
  q = length(s)
  my_sequence=seq(0,1,length.out=n_points)
  cox_direction = matrix(rep(NA_real_, q*n_points), ncol = q)
  for (j in (1:length(my_sequence))) {
    cox_direction[j,i] = my_sequence[j]
    delta = cox_direction[j,i] - s[i]
    for(k in setdiff(1:q, i)){ #setdiff leaves 2nd and 3rd elements of the row to
      # be transformed via cox-effect direction in our for loop below
      # if i=1 and q=3
      if(class(all.equal(s[i], 1)) != "character"){
        # cox_direction[j,k] = (1 - cox_direction[j,i])/(q-1)
        cox_direction[j,k] = (1 - my_sequence[j])/(q-1)
      } else{
        cox_direction[j,k] = s[k]*(1 - delta/(1 - s[i]))
      }

    }
  }
  return(cox_direction)
}

#' Computation of design matrix for Scheffe Model: Order 1
#'
#' \code{scheffe_model_order1} returns the model matrix for Scheffe Model (1)
#'
#' @param design specifies the model matrix
#' @examples
#' scheffe_model_order1(my_random_design)
#'
#' @export

scheffe_model_order1 <- function(design){
  q = ncol(design)
  n = nrow(design)
  parameters = q
  model_matrix = matrix(rep(NA_real_, parameters*n), nrow = n)
  model_matrix[,1:q] = design
  return(model_matrix)
}

#' Computation of model matrix for Scheffe Model: Order 2
#'
#' \code{scheffe_model_order2} returns the model matrix for Scheffe Model (2)
#'
#' @param design specifies the design matrix
#' @examples
#' scheffe_model_order2(my_random_design)
#'
#' @export

scheffe_model_order2 <- function(design){
  q = ncol(design)
  n = nrow(design)
  parameters = q + q*(q-1) / 2
  model_matrix = matrix(rep(NA_real_, parameters*n), nrow = n)
  model_matrix[,1:q] = design

  k = q
  for(i in 1:(q-1)){
    for(j in (i+1):q){
      # cat("i = ", i, ", j = ", j, "\n") # uncomment if needs to be printed
      k = k+1
      model_matrix[,k] = design[,i]*design[,j]
    }
  }
  return(model_matrix)
}

#' Computation of model matrix for Scheffe Model: Order 3. Built on SM-O2
#'
#' \code{scheffe_model_order3} returns the model matrix for Scheffe Model (3)
#'
#' @param design specifies the design matrix
#' @examples
#' scheffe_model_order3(my_random_design)
#'
#'@export

scheffe_model_order3 <- function(design) {
  q = ncol(design)
  model_matrix = scheffe_model_order2(design)
  for(i in 1:(q-2)){
    for(j in (i+1):(q-1)){
      for(k in (j+1):q){
        # cat("i = ", i, ", j = ", j, "\n") # uncomment if needs to be printed
        model_matrix = cbind(model_matrix, design[,i]*design[,j]*design[,k])
      }
    }
  }
  return(model_matrix)
}

#' Computation of model matrix for Scheffe Models: All Orders.
#'
#' \code{compute_model_matrix} returns the model matrix for Scheffe Model of any
#' order between 1 and 3
#'
#' @param design specifies the design matrix
#' @param order specifies the order of the model (1, 2 or 3)
#' @examples
#' compute_model_matrix(my_random_design,3)
#'
#' @export

compute_model_matrix <- function(design, order) {
  if(order == 1)
    printed_model_matrix = scheffe_model_order1(design)
  else{
    if(order == 2){
      printed_model_matrix = scheffe_model_order2(design)
    } else {
      printed_model_matrix = scheffe_model_order3(design)
    }
  }
  return(printed_model_matrix)
}

#' Compute D-optimality value for the pre-specified design and order
#'
#' \code{compute_d_optimality} computes D-optimality value for the
#' pre-specified design and order
#'
#' @param design specifies the design matrix
#' @param order specifies the order of the model (1, 2 or 3)
#' @examples
#' compute_d_optimality(my_random_design,3)
#'
#' @export


compute_d_optimality = function(design, order){
  if(order == 1){
    printed_model_matrix = compute_model_matrix(design, order = 1)

  } else{
    if(order == 2){
      printed_model_matrix = compute_model_matrix(design, order = 2)

    } else{
      if(order == 3){
        printed_model_matrix = compute_model_matrix(design, order = 3)
      }
    }
  }
  transposed_mm = t(printed_model_matrix)
  d_opt = - as.numeric(determinant(transposed_mm%*%printed_model_matrix)$modulus)
  # if(d_opt < 1e-16){
  #   stop("Matrix is not invertible")
  # }
  return(d_opt)
}


# compute_d_optimality = function(design, order){
#   if(order == 1){
#     printed_model_matrix = compute_model_matrix(design, order = 1)
#
#   } else{
#     if(order == 2){
#       printed_model_matrix = compute_model_matrix(design, order = 2)
#
#     } else{
#       if(order == 3){
#         printed_model_matrix = compute_model_matrix(design, order = 3)
#
#       }
#     }
#   }
#   d_opt = as.numeric(determinant(t(printed_model_matrix)%*%printed_model_matrix)$modulus)
#   if(d_opt < 1e-16){
#     stop("Matrix is not invertible")
#   }
#   return(d_opt)
# }




#' Perform coordinate exchange algorithm for the pre-specified design, order and
#' D-optimality criterion imposed
#'
#' \code{coordinate_exchange_d_optimal} computes I-optimality value for the
#' pre-specified design and order
#'
#' @param design specifies the design matrix
#' @param order specifies the order of the model (1, 2 or 3)
#' @param n_points specifies the number of Cox direction points
#' @param n_iter specifies the number of iterations
#'
#' @examples
#' coordinate_exchange_d_optimal(my_random_design,3,1000,10)
#'
#' @export

coordinate_exchange_d_optimal <- function(design, order, n_points, n_iter = 10) {
  d_optimal_design = design

  d_value_opt = compute_d_optimality(d_optimal_design, order)

  n=nrow(design)
  q=ncol(design)

  iter = 0
  while(iter < n_iter){
    iter = iter + 1

    for (row in 1:n) {

      for (ing in 1:q) {
        # cat("row: ", row, "ing: ", ing, "\n")
        cox_dir = compute_cox_direction(d_optimal_design[row,], ing, n_points)
        for (entry in 1:nrow(cox_dir)) {
          new_design = d_optimal_design
          entry=1
          new_design[row,] = cox_dir[entry,]
          d_value_opt = compute_d_optimality(d_optimal_design, order)
          d_value_new = compute_d_optimality(new_design, order)

          if (d_value_opt >= d_value_new){
            d_optimal_design = new_design
            d_value_opt = d_value_new
          }
          # else {
          #   new_design=new_design
          # }
        }
        # cat("d_value_opt: ", d_value_opt, "\n\n")
        # d_optimal_design
        # cat("d_value_new: ", d_value_new, "\nd_value_opt: ", d_value_opt, "\n\n")
      }
    }
  }
  return(d_optimal_design)
}

#' Create a moment matrix for I-optimality based design
#'
#' \code{create_moment_matrix} Create a moment matrix for I-optimality based design
#' for model of order 1, 2, or 3 with q ingredients.
#'
#' @param q the number of ingredients in a mixture
#' @param order specifies the order of the model (1, 2 or 3)
#'
#' @export

create_moment_matrix = function(q, order){

  # q=ncol(design) # disregard, make it simplier
  # define the number of parameters for model order \in 1, 2, 3

  if(order == 1){
    parameters = q
  } else{
    if(order == 2){
      parameters = q + q*(q-1) / 2
    } else{
      if(order == 3){
        parameters = q + q*(q-1) / 2 + q*(q-1)*(q-2) / 6
      }
    }
  }

  # Auxiliary matrix F is of size parameters with q columns;
  # We fill it in with 1 if the parameter is assumed to be in the moment matrix

  auxiliary_matrix_f = matrix(0, nrow = parameters, ncol = q)
  counter = 0
  for(i in 1:q){
    counter = counter + 1
    auxiliary_matrix_f[counter, i] = 1
  }

  # Indeces for order 2:
  # \sum_{i}^{q-1}\sum_{j=i+1}^{q}
  # j=i+1 and i=1

  if(order>= 2){
    for(i in 1:(q-1)){
      for(j in (i+1):q){
        counter = counter + 1
        auxiliary_matrix_f[counter, i] = 1
        auxiliary_matrix_f[counter, j] = 1
      }
    }
  }

  # Indeces for order 3:
  # \sum_{i}^{q-2}\sum_{j=i+1}^{q-1}\sum_{k=j+1}^{q}
  # i=1, j=i+1, k=j+1
  # each ith, jth, kth element are filled in with 1

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

  moment_matrix = matrix(NA_real_, nrow = parameters, ncol = parameters)

  for(i in 1:parameters){
    for(j in 1:parameters){

      # auxiliary_vector is the sum of ith and jth row of auxiliary matrix F

      auxiliary_vector = auxiliary_matrix_f[i, 1:q] + auxiliary_matrix_f[j, 1:q]
      numerator = prod(factorial(auxiliary_vector))
      denominator = factorial(q-1 + sum(auxiliary_vector))
      moment_matrix[i,j] = numerator / denominator
    }
  }
  return(moment_matrix)
}

#' Compute I-optimality value for the pre-specified design and order
#'
#' \code{compute_i_optimality} computes I-optimality value for the
#' pre-specified design and order
#'
#' @param design specifies the design matrix
#' @param order specifies the order of the model (1, 2 or 3)
#' @examples
#' compute_i_optimality(my_random_design,3)
#'
#' @export

compute_i_optimality = function(design, order){
  q = ncol(design)

  printed_model_matrix = compute_model_matrix(design, order)
  printed_moments_matrix = create_moment_matrix(q, order)

  i_opt = tr(solve(t(printed_model_matrix)%*%printed_model_matrix)%*%printed_moments_matrix)
  return(i_opt)
}

#' Perform coordinate exchange algorithm for the pre-specified design, order and
#' I-optimality criterion imposed
#'
#' \code{coordinate_exchange_i_optimal} computes I-optimality value for the
#' pre-specified design and order
#'
#' @param design specifies the design matrix
#' @param order specifies the order of the model (1, 2 or 3)
#' @param n_points specifies the number of Cox direction points
#' @param n_iter specifies the number of iterations
#'
#' @examples
#' coordinate_exchange_i_optimal(my_random_design,3,1000,10)
#'
#' @export

coordinate_exchange_i_optimal <- function(design, order, n_points, n_iter = 10) {
  i_optimal_design = design

  i_value_opt = compute_i_optimality(i_optimal_design, order)

  n=nrow(design)
  q=ncol(design)

  iter = 0
  while(iter < n_iter){
    iter = iter + 1

    for (row in 1:n) {

      for (ing in 1:q) {
        # cat("row: ", row, "ing: ", ing, "\n")
        cox_dir = compute_cox_direction(i_optimal_design[row,], ing, n_points)
        for (entry in 1:nrow(cox_dir)) {
          new_design = i_optimal_design
          entry=1
          new_design[row,] = cox_dir[entry,]
          i_value_opt = compute_i_optimality(i_optimal_design, order)
          i_value_new = compute_i_optimality(new_design, order)

          if (i_value_opt >= i_value_new){
            i_optimal_design = new_design
            i_value_opt = i_value_new
          }
          # else {
          #   new_design=new_design
          # }
        }
        # cat("d_value_opt: ", d_value_opt, "\n\n")
        # d_optimal_design
        # cat("d_value_new: ", d_value_new, "\nd_value_opt: ", d_value_opt, "\n\n")
      }
    }
  }

  return(i_optimal_design)
}

#' Perform coordinate exchange algorithm for the pre-specified design, order and
#' optimality criterion
#'
#' \code{coordinate_exchange_general} computes I-optimality value for the
#' pre-specified design and order
#'
#' @param design specifies the design matrix
#' @param order specifies the order of the model (1, 2 or 3)
#' @param n_points specifies the number of Cox direction points
#' @param optimality_criterion specifies the optimality criterion ("D" or "I")
#' @param n_iter specifies the number of iterations
#'
#' @examples
#' coordinate_exchange_general(my_random_design,3,1000,I,100)
#'
#' @export


coordinate_exchange_general<-function(design, order, n_points, optimality_criterion, n_iter = 100) {
  if (optimality_criterion =="D"){
    coordinate_exchange_d_optimal(design, order, n_points, n_iter = 100)
  } else {
    coordinate_exchange_i_optimal(design, order, n_points, n_iter = 100)
  }
}















