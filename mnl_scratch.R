#' Create a random initial design (multinomial logit model)
#'
#' \code{create_random_initial_design_discrete} creates a random initial design
#' for multinomial logit model with q ingredients, j alternatives, and s choice
#' sets
#' @param q the number of ingredients in a mixture
#' @param alternatives the number of alternatives within each choice set
#' @param n_choice_sets the number of choice sets
#' @param seed set seed for reproducibility
#'
#' @export

create_random_initial_design_discrete = function(n_ingredients, n_alternatives, n_choice_sets,
    seed) {
    # alternatives=j # does not work n_choice_sets=s # does not work, so use letter
    # instead
    if (!is.null(seed))
        set.seed(seed)
    my_array = array(rep(NA_real_, n_ingredients * n_alternatives * n_choice_sets),
        dim = c(n_ingredients, n_alternatives, n_choice_sets))
    for (j in 1:n_alternatives) {
        for (k in 1:n_choice_sets) {
            for (m in 1:n_ingredients) {
                random_values = runif(n_ingredients)
                my_array[, j, k] = random_values/sum(random_values)
            }
        }
    }
    return(my_array)
}


MNL_compute_choice_probabilites = function(design, beta, order) {
    #
    q = dim(design)[1]
    J = dim(design)[2]
    S = dim(design)[3]


    p1 = q - 1
    p2 = q * (q - 1)/2
    p3 = q * (q - 1) * (q - 2)/6

    if (order == 1) {
        n_parameters = p1
    } else {
        if (order == 2) {
            n_parameters = p1 + p2
        } else {
            if (order == 3) {
                n_parameters = p1 + p2 + p3
            }
        }
    }



    # beta = rep(1,)
    beta_star = beta[1:p1]  # modify
    beta_2FI = beta[(p1 + 1):p2]
    beta_3FI = beta[(p2 + 1):p3]

    stopifnot(length(beta) == n_parameters)
    stopifnot(order %in% 1:3)

    # define utility (matrix) of alternative j in choice set s for (i in 1:q){
    U = matrix(rep(NA_real_, J * S), ncol = S)
    for (j in 1:J) {
        for (s in 1:S) {

            # ORDER-1 X_JS

            x_js = design[, j, s]
            U_js_term1 = sum(x_js[1:p1] * beta_star)
            u_js = U_js_term1


            if (order >= 2) {
                # ORDER-2 X_IJS*X_KJS - 2FI terms
                x_second_order = rep(NA_real_, p2)
                counter = 0
                for (i in 1:(q - 1)) {
                  for (k in (i + 1):q) {
                    counter = counter + 1
                    x_ijs_x_kjs = design[i, j, s] * design[k, j, s]
                    x_second_order[counter] = x_ijs_x_kjs
                  }
                }
                U_js_term2 = sum(x_second_order * beta_2FI)
                u_js = u_js + U_js_term2
            }
            # x_ijs_x_kjs = design[i,j,s]*design[k,j,s]
            if (order >= 3) {
                # 0RDER-3 X_IJS*X_KJS*XLJS - 3FI terms
                x_third_order = rep(NA_real_, p3)
                counter = 0
                # i = 1 k = i+1 l = k+1
                for (i in 1:(q - 2)) {
                  for (k in (i + 1):(q - 1)) {
                    for (l in (k + 1):q) {
                      counter = counter + 1
                      x_ijs_x_kjs_xljs = design[i, j, s] * design[k, j, s] * design[l,
                        j, s]
                      x_third_order[counter] = x_ijs_x_kjs_xljs
                    }
                  }
                }
                U_js_term3 = sum(x_third_order * beta_3FI)
                u_js = u_js + U_js_term3
            }
            # x_ijs_x_kjs_xljs = design[i,j,s]*design[k,j,s]*design[l,j,s]


            # U_js_term2 = sum(x_ijs_x_kjs[(p1+1):p2]*beta_2FI) U_js_term2 =
            # sum(x_second_order*beta_2FI) U_js_term3 = sum(x_ijs_x_kjs[(p2+1):p3]*beta_3FI)
            # U_js_term3 = sum(x_third_order*beta_3FI) u_js = U_js_term1 + U_js_term2 +
            # U_js_term3
            U[j, s] = u_js  # matrix of model expansion
        }
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

MNL_compute_model_matrix <- function(design, order) {
    q = dim(design)[1]
    J = dim(design)[2]
    S = dim(design)[3]


    p1 = q - 1
    p2 = q * (q - 1)/2
    p3 = q * (q - 1) * (q - 2)/6

    if (order == 1) {
        n_parameters = p1
    } else {
        if (order == 2) {
            n_parameters = p1 + p2
        } else {
            if (order == 3) {
                n_parameters = p1 + p2 + p3
            }
        }
    }



    model_array = array(rep(NA_real_, n_parameters * J * S), dim = c(n_parameters,
        J, S))

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
    if (order >= 3) {
        # counter = 0
        for (i in 1:(q - 2)) {
            for (j in (i + 1):(q - 1)) {
                for (k in (j + 1):q) {
                  counter = counter + 1
                  model_array[counter, , ] = design[i, , ] * design[j, , ] * design[k,
                    , ]
                }
            }
        }
    }

    return(model_array)

}


MNL_compute_information_matrix = function(design, beta, order) {
    q = dim(design)[1]
    J = dim(design)[2]
    S = dim(design)[3]

    p1 = q - 1
    p2 = q * (q - 1)/2
    p3 = q * (q - 1) * (q - 2)/6

    if (order == 1) {
        n_parameters = p1
    } else {
        if (order == 2) {
            n_parameters = p1 + p2
        } else {
            if (order == 3) {
                n_parameters = p1 + p2 + p3
            }
        }
    }

    model_array = MNL_compute_model_matrix(design, order)
    P = MNL_compute_choice_probabilites(design, beta, order)

    information_matrix = matrix(rep(0, n_parameters * n_parameters), ncol = n_parameters)  # generalize this
    for (s in 1:S) {
        p_s = P[, s]
        I_s = model_array[1:n_parameters, , s] %*% (diag(p_s) - (t(t(p_s)) %*% t(p_s))) %*%
            t(model_array[1:n_parameters, , s])
        information_matrix = information_matrix + I_s
    }
    return(information_matrix)
}

# _______________________________________________________________________________
# Check with Mario

MNL_compute_d_optimality = function(design, beta, order) {
    mnl_information_matrix = MNL_compute_information_matrix(design, beta, order)
    log_determinant = determinant(t(mnl_information_matrix) %*% mnl_information_matrix)$modulus
    return(as.numeric(log_determinant))
}


MNL_coordinate_exchange_d_optimal <- function(design, order, beta, n_points, n_iter = 10) {
    q = dim(design)[1]
    J = dim(design)[2]
    S = dim(design)[3]

    n_parameters = q - 1 + q * (q - 1)/2 + q * (q - 1) * (q - 2)/6

    d_optimal_design = design
    d_value_opt = MNL_compute_d_optimality(design, beta, order)

    stopifnot(length(beta) == n_parameters)
    stopifnot(order %in% 1:3)

    iter = 0
    while (iter < n_iter) {
        iter = iter + 1

        for (j in 1:J) {
            for (s in 1:S) {
                for (ing in 1:q) {
                  # cat('row: ', row, 'ing: ', ing, '\n')
                  cox_dir = compute_cox_direction(d_optimal_design[, j, s], ing,
                    n_points)
                  for (entry in 1:nrow(cox_dir)) {
                    new_design = d_optimal_design
                    entry = 1
                    new_design[, j, s] = cox_dir[entry, ]
                    d_value_opt = MNL_compute_d_optimality(d_optimal_design, beta,
                      order)
                    d_value_new = MNL_compute_d_optimality(new_design, beta, order)

                    if (d_value_opt <= d_value_new) {
                      d_optimal_design = new_design
                      d_value_opt = d_value_new
                    }
                    # else { new_design=new_design }
                  }
                  # cat('d_value_opt: ', d_value_opt, '\n\n') d_optimal_design cat('d_value_new:
                  # ', d_value_new, '\nd_value_opt: ', d_value_opt, '\n\n')

                }
            }
        }
    }
    return(d_optimal_design)
}

MNL_create_moment_matrix(3, 3)

MNL_create_moment_matrix = function(q, order) {

    # q=ncol(design) # disregard, make it simplier define the number of parameters
    # for model order \in 1, 2, 3

    if (order == 1) {
        parameters = q - 1
    } else {
        if (order == 2) {
            parameters = q - 1 + q * (q - 1)/2
        } else {
            if (order == 3) {
                parameters = q - 1 + q * (q - 1)/2 + q * (q - 1) * (q - 2)/6
            }
        }
    }

    # Auxiliary matrix F is of size parameters with q columns; We fill it in with 1
    # if the parameter is assumed to be in the moment matrix

    auxiliary_matrix_f = matrix(rep(0, parameters * q), ncol = q)
    counter = 0
    for (i in 1:(q - 1)) {
        counter = counter + 1
        auxiliary_matrix_f[counter, i] = 1
    }

    # Indeces for order 2: \sum_{i}^{q-1}\sum_{j=i+1}^{q} j=i+1 and i=1

    if (order >= 2) {
        for (i in 1:(q - 1)) {
            for (j in (i + 1):q) {
                counter = counter + 1
                auxiliary_matrix_f[counter, i] = 1
                auxiliary_matrix_f[counter, j] = 1
            }
        }
    }

    # Indeces for order 3: \sum_{i}^{q-2}\sum_{j=i+1}^{q-1}\sum_{k=j+1}^{q} i=1,
    # j=i+1, k=j+1 each ith, jth, kth element are filled in with 1

    if (order >= 3) {
        for (i in 1:(q - 2)) {
            for (j in (i + 1):(q - 1)) {
                for (k in (j + 1):q) {
                  counter = counter + 1
                  auxiliary_matrix_f[counter, i] = 1
                  auxiliary_matrix_f[counter, j] = 1
                  auxiliary_matrix_f[counter, k] = 1
                }
            }
        }
    }

    MNL_moment_matrix = matrix(rep(NA_real_, parameters * parameters), ncol = parameters)

    for (i in 1:parameters) {
        for (j in 1:parameters) {

            # auxiliary_vector is the sum of ith and jth row of auxiliary matrix F

            auxiliary_vector = auxiliary_matrix_f[i, ] + auxiliary_matrix_f[j, ]
            numerator = prod(factorial(auxiliary_vector))
            denominator = factorial(q - 1 + sum(auxiliary_vector))
            MNL_moment_matrix[i, j] = numerator/denominator
        }
    }
    return(MNL_moment_matrix)
}


MNL_compute_i_optimality = function(design, order, beta) {
    q = dim(design)[1]
    printed_information_matrix = MNL_compute_information_matrix(design, beta, order)
    printed_moments_matrix = MNL_create_moment_matrix(q, order)
    # MNL_i_opt = inv(printed_information_matrix) %*% printed_moments_matrix
    MNL_i_opt = sum(diag(solve(printed_information_matrix, printed_moments_matrix)))  # trace is the sum of the diagonal elements
    return(MNL_i_opt)
}


MNL_coordinate_exchange_i_optimal <- function(design, order, beta, n_points, n_iter = 10) {
    q = dim(design)[1]
    J = dim(design)[2]
    S = dim(design)[3]

    n_parameters = q - 1 + q * (q - 1)/2 + q * (q - 1) * (q - 2)/6

    i_optimal_design = design
    i_value_opt = MNL_compute_i_optimality(design, beta, order)

    stopifnot(length(beta) == n_parameters)
    stopifnot(order %in% 1:3)

    iter = 0
    while (iter < n_iter) {
        iter = iter + 1

        for (j in 1:J) {
            for (s in 1:S) {
                for (ing in 1:q) {
                  # cat('row: ', row, 'ing: ', ing, '\n')
                  cox_dir = compute_cox_direction(i_optimal_design[, j, s], ing,
                    n_points)
                  for (entry in 1:nrow(cox_dir)) {
                    new_design = i_optimal_design
                    entry = 1
                    new_design[, j, s] = cox_dir[entry, ]
                    i_value_opt = MNL_compute_i_optimality(i_optimal_design, beta,
                      order)
                    i_value_new = MNL_compute_i_optimality(new_design, beta, order)

                    if (i_value_opt >= i_value_new) {
                      i_optimal_design = new_design
                      i_value_opt = i_value_new
                    }
                    # else { new_design=new_design }
                  }
                  # cat('d_value_opt: ', d_value_opt, '\n\n') d_optimal_design cat('d_value_new:
                  # ', d_value_new, '\nd_value_opt: ', d_value_opt, '\n\n')

                }
            }
        }
    }
    return(i_optimal_design)
}


MNL_coordinate_exchange_general <- function(design, order, n_points, beta, optimality_criterion,
    n_iter = 100) {
    if (optimality_criterion == "D") {
        MNL_coordinate_exchange_d_optimal(design, order, n_points, n_iter = 100)
    } else {
        MNL_coordinate_exchange_i_optimal(design, order, n_points, n_iter = 100)
    }
}



