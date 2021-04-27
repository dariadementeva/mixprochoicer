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




MNL_compute_choice_probabilites = function(design, m, beta, order) {
  #
  q = dim(design)[1] - m
  J = dim(design)[2]
  S = dim(design)[3]


  p1 = q - 1
  p2 = q * (q - 1)/2
  p3 = q * (q - 1) * (q - 2)/6
  p4 = q*(q-1)/2 + q*m + m*(m-1)/2 + m

  if (order == 1) {
    n_parameters = p1
  } else {
    if (order == 2) {
      n_parameters = p1 + p2
    } else {
      if (order == 3) {
        n_parameters = p1 + p2 + p3
      }
      else {
        n_parameters = p1 + p4
      }
    }
  }



  # beta = rep(1,)
  beta_star = beta[1:p1]  # modify
  beta_2FI = beta[(p1 + 1):p2]
  beta_3FI = beta[(p2 + 1):p3]
  beta_PV = beta[(p3 + 1):p4]


  # stopifnot(length(beta) == n_parameters) stopifnot(order %in% 1:3)

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

      if (order >= 3) {
        # 0RDER-3 X_IJS*X_KJS*XLJS - 3FI terms
        x_third_order = rep(NA_real_, p3)
        counter = 0
        for (i in 1:(q - 2)) {
          for (k in (i + 1):(q - 1)) {
            for (l in (k + 1):q) {
              counter = counter + 1
              x_ijs_x_kjs_xljs = design[i, j, s] * design[k, j, s] * design[l, j, s]
              x_third_order[counter] = x_ijs_x_kjs_xljs
            }
          }
        }
        U_js_term3 = sum(x_third_order * beta_3FI)
        u_js = u_js + U_js_term3
      }

      if (order >= 4) {
        x_order4 = rep(NA_real_, p4)





        u_js = u_js + U_js_term4
      }

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

MNL_MPV_compute_model_matrix <- function(design, m) {
  q = dim(design)[1] - m
  J = dim(design)[2]
  S = dim(design)[3]

  n_parameters = q - 1 + q*(q-1)/2 + q*m + m*(m-1)/2 + m

  model_array = array(rep(NA_real_, n_parameters * J * S), dim = c(n_parameters, J, S))

  # Main effects
  model_array[1:(q - 1), , ] = design[1:(q - 1), , ]

  # 2FI
    for (i in 1:(q - 1)) {
      for (j in (i + 1):q) {
        counter = q - 1 + 1
        model_array[counter, , ] = design[i, , ] * design[j, , ]
      }
    }

  #MPV interactions

    for (l in 1:m) {
      for (k in 1:q){
        counter =  ### SET COUNTER FROM 2FI
        model_array[counter, ,] = design[l, , ] * design[k, , ]
      }
    }
    return(model_array)
}

    #PV interactions
    w = 0
    for (e in 1:(m-1)) {
      for (r in (i+1):m){
        w = w+1
        model_array[w, , ] =  model_array[e, , ]* model_array[r, , ]
      }
    }

    #PV quadratic effects
    t = 0
    for (z in 1:m) {
      t = t+1
      model_array[r, , ] =  model_array[z, , ]*model_array[z, , ]
    }

  }
return(model_array)
  }


MNL_MPV_compute_information_matrix = function(design, m, beta, order=4) {
  q = dim(design)[1] - m
  J = dim(design)[2]
  S = dim(design)[3]


  p1 = q - 1
  p2 = q * (q - 1)/2
  p3 = q * (q - 1) * (q - 2)/6
  p4 = q*(q-1)/2 + q*m + m*(m-1)/2 + m

  if (order == 1) {
    n_parameters = p1
  } else {
    if (order == 2) {
      n_parameters = p1 + p2
    } else {
      if (order == 3) {
        n_parameters = p1 + p2 + p3
      }
    } else {
      if (order==4)
      n_parameters = p1 + p4
    }
  }

  model_array = MNL_compute_model_matrix(design, m, order)
  P = MNL_compute_choice_probabilites(design, m, beta, order)

  information_matrix = matrix(rep(0, n_parameters * n_parameters), ncol = n_parameters)  # generalize this
  for (s in 1:S) {
    p_s = P[, s]
    I_s = model_array[1:n_parameters, , s] %*% (diag(p_s) - (t(t(p_s)) %*% t(p_s))) %*% t(model_array[1:n_parameters, , s])
    information_matrix = information_matrix + I_s
  }
  return(information_matrix)
}

