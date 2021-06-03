## quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' MM Algorithm for M-HPCA
#'
#' @param MHPCA MHPCA object
#' @param group group
#' @param maxiter maximum iterations
#' @param epsilon epsilon value
#' @param j_tot total repetitions
#'
#'
mm = function(MHPCA, group, maxiter = 1000, epsilon = 1e-6, j_tot) {

  subjects        = unique(MHPCA$data[[group]]$Subject)
  names(subjects) = subjects

  # determining whether there are missing visits/obs
  complete = purrr::map(
    subjects,
    ~ MHPCA$data[[group]][MHPCA$data[[group]]$Subject == subjects[.x], ] %>%
      dplyr::pull(Repetition) %>% unique()
  ) %>%
    purrr::map_dbl(~ length(unique(.x)) == j_tot) %>%
    sum() == length(subjects)

  if (complete) {
    mm_complete(MHPCA, group, maxiter, epsilon)
  } else {
    mm_incomplete(MHPCA, group, maxiter, epsilon)
  }
}


# MM Complete Data --------------------------------------------------------

#' MM Algorithm for M-HPCA with complete data
#'
#' @param MHPCA MHPCA object
#' @param group group
#' @param maxiter maximum iterations
#' @param epsilon epsilon value
#'
#'
mm_complete = function(MHPCA, group, maxiter, epsilon) {

  # inputs ----
  subjects        = unique(MHPCA$data[[group]]$Subject)
  names(subjects) = subjects

  # store the number of obs, regions, and timepoints
  obs        = unique(MHPCA$data[[group]]$Repetition)
  names(obs) = obs
  J          = length(obs)
  r_tot      = length(unique(MHPCA$data[[group]]$reg))
  f_tot      = length(unique(MHPCA$data[[group]]$func))
  TR         = f_tot * r_tot
  TRJ        = TR * J

  n = length(subjects)
  G = MHPCA$data[[group]] %>%
    dplyr::select(dplyr::starts_with("between")) %>%
    ncol()
  H = MHPCA$data[[group]] %>%
    dplyr::select(dplyr::starts_with("within")) %>%
    ncol()
  HJ = J * H

  # create identity and ones matrices for kronecker products
  one_J  = rep(1, J)
  ones_J = matrix(1, J, J)
  I_J    = Matrix::Diagonal(J, 1)
  I_TR   = Matrix::Diagonal(TR, 1)
  I_HJ   = Matrix::Diagonal(HJ, 1)
  I_G    = Matrix::Diagonal(G, 1)

  y_i = MHPCA$data[[group]] %>%
    split(.$Subject) %>%
    map(~ pull(.x, y_centered2))

  Phi1 = MHPCA$data[[group]] %>%
    dplyr::select(dplyr::starts_with("between")) %>%
    dplyr::slice(1:TR) %>%
    as.matrix(dimnames = list(names(.), NULL))

  level1_names = colnames(Phi1)

  z_i = Matrix::kronecker(one_J, Phi1)

  Phi2 = MHPCA$data[[group]] %>%
    dplyr::select(dplyr::starts_with("within")) %>%
    dplyr::slice(1:TR) %>%
    as.matrix(dimnames = list(names(.), NULL))

  level2_names = colnames(Phi2)

  Phi2t = Matrix::t(Phi2)

  w_i = Matrix::kronecker(I_J, Phi2)

  wt = Matrix::t(w_i)
  wtw = Matrix::kronecker(I_J, crossprod(Phi2))
  ztz = J * Matrix::crossprod(Phi1)
  wtz = Matrix::kronecker(one_J, Matrix::crossprod(Phi2, Phi1))

  # storage ----
  sigma2      = rep(NA_real_, maxiter + 1)
  loglik_vec  = rep(NA_real_, maxiter + 1)
  Lambda1     = vector(mode = "list", length = maxiter + 1)
  Lambda2     = vector(mode = "list", length = maxiter + 1)
  zeta_i      = vector(mode = "list", length = maxiter + 1)
  xi_i        = vector(mode = "list", length = maxiter + 1)

  # initialize variances ----
  Lambda1[[1]] = Matrix::Diagonal(G, 1)
  Lambda2[[1]] = Matrix::Diagonal(H, 1)
  sigma2[[1]]  = MHPCA$marg$within$functional[[group]]$sigma_2d

  E_mat = Lambda1[[1]]
  F_mat = Matrix::kronecker(I_J, Lambda2[[1]])

  E_sqrt = sqrt(E_mat)
  F_sqrt = sqrt(F_mat)

  # estimated level 1 covariance
  K_dB = Matrix::tcrossprod(Phi1, Phi1 %*% Lambda1[[1]])

  # Qinv = sigma2_d * I_TR + K_dW
  Q = 1 / sigma2[1] * I_TR - 1 / sigma2[1] ^ 2 * Phi2 %*% Matrix::solve(Matrix::solve(Lambda2[[1]]) + 1 / sigma2[1] * Matrix::crossprod(Phi2), Phi2t)

  # find inverse of JK_dB + Qinv
  Q_Phi1 = Q %*% Phi1
  JK_Qinv_inv = Q - J * Q_Phi1 %*% Matrix::solve(Matrix::solve(Lambda1[[1]]) + J * Matrix::crossprod(Phi1, Q_Phi1), Matrix::crossprod(Phi1, Q))

  temp =  Q %*% K_dB %*% JK_Qinv_inv
  omega_inv = Matrix::kronecker(I_J, Q) - Matrix::kronecker(ones_J, temp)

  omegainv_y    = vector(mode = "list", length = n)
  yt_omegainv_y = vector(length = n)

  for (i in 1:n) {
    omegainv_y[[i]] = omega_inv %*% y_i[[i]]
    yt_omegainv_y[i] = as.numeric(Matrix::crossprod(y_i[[i]], omegainv_y[[i]]))
  }

  loglik_vec[1] = loglik_calc_complete(
    sigma2[1], TRJ, n, z_i, w_i, E_sqrt, F_mat, F_sqrt, I_G, I_HJ,
    omegainv_y, yt_omegainv_y, wtw, ztz, wtz
  )

  # iterate for maxiter ----
  for (k in 1:maxiter) {

    # calculate scores
    zeta_i[[k]]        = vector(mode = "list", length = n)
    names(zeta_i[[k]]) = subjects
    xi_i[[k]]          = vector(mode = "list", length = n)
    names(xi_i[[k]])   = subjects

    y_sum = 0

    Lambda1_Zt    = Matrix::tcrossprod(Lambda1[[k]], z_i)
    Lambda2_Phi2t = Matrix::kronecker(I_J, Matrix::tcrossprod(Lambda2[[k]], Phi2))
    for (i in 1:n) {
      zeta_i[[k]][[i]] = Lambda1_Zt %*% omegainv_y[[i]]
      xi_i[[k]][[i]]   = Lambda2_Phi2t %*% omegainv_y[[i]]
      rownames(zeta_i[[k]][[i]]) = level1_names
      rownames(xi_i[[k]][[i]])   = rep(level2_names, each = J)
      y_sum = y_sum + Matrix::crossprod(omegainv_y[[i]])
    }

    trace_sum = n * sum(Matrix::diag(omega_inv))

    sigma2[k + 1] = as.vector(sigma2[k] * sqrt(y_sum / trace_sum))

    Lambda1[[k + 1]] = Matrix::Diagonal(G, x = 0)
    Lambda2[[k + 1]] = Matrix::Diagonal(H, x = 0)

    A = Matrix::crossprod(z_i, omega_inv %*% z_i)
    for (j in 1:G) {
      gamma_sum = 0
      zeta_sum  = 0
      for (i in 1:n) {
        gamma_sum = gamma_sum + A[j, j]
        zeta_sum = zeta_sum + zeta_i[[k]][[i]][j]^2
      }
      Lambda1[[k + 1]][j, j] = sqrt(zeta_sum / gamma_sum)
    }

    B = Matrix::crossprod(w_i, omega_inv %*% w_i)
    for (h in 1:H) {
      gamma_sum = 0
      xi_sum    = 0
      for (j in 1:J) {
        index = h + (j - 1) * H
        gamma_sum = gamma_sum + n * B[index, index]
        for (i in 1:n) {
          xi_sum = xi_sum + xi_i[[k]][[i]][index] ^ 2
        }
      }
      Lambda2[[k + 1]][h, h] = sqrt(xi_sum / gamma_sum)
    }

    E_mat = Lambda1[[k + 1]]
    F_mat = Matrix::kronecker(I_J, Lambda2[[k + 1]])
    E_sqrt = sqrt(E_mat)
    F_sqrt = sqrt(F_mat)

    # estimated level 1 covariance
    phi_lambda = Phi1 %*% Lambda1[[k + 1]]
    K_dB = Matrix::tcrossprod(phi_lambda, Phi1)

    # Qinv = sigma2_d * I_TR + K_dW
    Q = 1 / sigma2[k + 1] * I_TR - 1 / sigma2[k + 1] ^ 2 * Phi2 %*% Matrix::solve(Matrix::solve(Lambda2[[k + 1]]) + 1 / sigma2[k + 1] * Matrix::crossprod(Phi2), Phi2t)

    # find inverse of JK_dB + Qinv
    Q_Phi1 = Q %*% Phi1
    JK_Qinv_inv = Q - J * Q_Phi1 %*% Matrix::solve(Matrix::solve(Lambda1[[k + 1]]) + J * Matrix::crossprod(Phi1, Q_Phi1), Matrix::crossprod(Phi1, Q))

    temp = Q %*% K_dB %*% JK_Qinv_inv
    omega_inv = Matrix::kronecker(I_J, Q) - Matrix::kronecker(ones_J, temp)

    omegainv_y    = vector(mode = "list", length = n)
    yt_omegainv_y = vector(length = n)
    for (i in 1:n) {
      omegainv_y[[i]] = omega_inv %*% y_i[[i]]
      yt_omegainv_y[i] = as.numeric(Matrix::crossprod(y_i[[i]], omegainv_y[[i]]))
    }

    loglik_vec[k + 1] = loglik_calc_complete(
      sigma2[k + 1], TRJ, n, z_i, w_i, E_sqrt, F_mat, F_sqrt, I_G, I_HJ,
      omegainv_y, yt_omegainv_y, wtw, ztz, wtz
    )

    if (abs(loglik_vec[k + 1] - loglik_vec[k]) <
        epsilon * abs(loglik_vec[k + 1] + 1)) {
      break
    }
  }

  return(list(
    "level1_components" = names(dplyr::select(MHPCA$data[[group]], dplyr::starts_with("between"))),
    "level2_components" = names(dplyr::select(MHPCA$data[[group]], dplyr::starts_with("within"))),
    # "loglik_vec" = loglik_vec[!is.na(sigma2)],
    # "Lambda1"    = Lambda1[!sapply(Lambda1, is.null)],
    # "Lambda2"    = Lambda2[!sapply(Lambda2, is.null)],
    # "sigma2"     = sigma2[!is.na(sigma2)],
    # "zeta_i"     = zeta_i[!sapply(zeta_i, is.null)],
    # "xi_i"       = xi_i[!sapply(xi_i, is.null)]
    "loglik_vec" = loglik_vec[(k + 1)],
    "Lambda1"    = Lambda1[[(k + 1)]],
    "Lambda2"    = Lambda2[[(k + 1)]],
    "sigma2"     = sigma2[(k + 1)],
    "zeta_i"     = zeta_i[[k]],
    "xi_i"       = xi_i[[k]]
  ))
}

#' Log-likelihood calculation for MM algorithm with complete data
#'
#' @param sigma2 sigma2
#' @param TRJ TRJ
#' @param n n
#' @param z_i z_i
#' @param w_i w_i
#' @param E_sqrt E_sqrt
#' @param F_mat F_mat
#' @param F_sqrt F_sqrt
#' @param I_G I_G
#' @param I_HJ I_HJ
#' @param omegainv_y omegainv_y
#' @param yt_omegainv_y yt_omegainv_y
#' @param wtw wtw
#' @param ztz ztz
#' @param wtz wtz
#'
loglik_calc_complete = function(
  sigma2, TRJ, n, z_i, w_i, E_sqrt, F_mat, F_sqrt, I_G, I_HJ,
  omegainv_y, yt_omegainv_y, wtw, ztz, wtz
) {

  FWtWF = F_sqrt %*% wtw %*% F_sqrt
  EZtZE = E_sqrt %*% ztz %*% E_sqrt
  WtZE  = wtz %*% E_sqrt
  F_inv = Matrix::solve(F_mat)

  middle_inv = I_G + 1 / sigma2 * (EZtZE - Matrix::crossprod(WtZE, Matrix::solve(sigma2 * F_inv + wtw, WtZE)))

  logdet_level2_arg = I_HJ + 1 / sigma2 * FWtWF
  logdet_level2 = Matrix::determinant(logdet_level2_arg, logarithm = TRUE)$modulus
  logdet_level1 = Matrix::determinant(middle_inv, logarithm = TRUE)$modulus
  sum_yt_omegainv_y = sum(yt_omegainv_y)

  loglik = -n / 2 * (TRJ * log(sigma2) + logdet_level2 + logdet_level1 + sum_yt_omegainv_y)

  return(loglik)
}







# MM Incomplete Data ------------------------------------------------------


#' MM Algorithm for M-HPCA with incomplete data
#'
#' @param MHPCA MHPCA object
#' @param group group
#' @param maxiter maximum iterations
#' @param epsilon epsilon value
#'
#'
mm_incomplete = function(MHPCA, group, maxiter, epsilon) {

  # to avoid issues with non-standard evaluation in tidyeval, set "global
  # variables" to NULL and remove them. this won't cause an issue with the rest
  # of the code.
  Repetition <- NULL

  rm(list = c("Repetition"))


  # inputs ----
  subjects        = unique(MHPCA$data[[group]]$Subject)
  names(subjects) = subjects

  # store the number of obs, regions, and timepoints
  obs        = unique(MHPCA$data[[group]]$Repetition)
  names(obs) = obs
  J          = length(obs)
  r_tot      = length(unique(MHPCA$data[[group]]$reg))
  f_tot      = length(unique(MHPCA$data[[group]]$func))
  TR         = r_tot * f_tot
  TRJ        = TR * J

  n  = length(subjects)
  G  = MHPCA$data[[group]] %>%
    dplyr::select(dplyr::starts_with("between")) %>%
    ncol()
  H  = MHPCA$data[[group]] %>%
    dplyr::select(dplyr::starts_with("within")) %>%
    ncol()
  HJ = J * H

  # create identity and ones matrices for kronecker products
  one_J  = purrr::map(1:J, ~ rep(1, .x))
  ones_J = purrr::map(1:J, ~ matrix(1, .x, .x))
  I_J    = purrr::map(1:J, ~ Matrix::Diagonal(.x, 1))
  I_TR   = Matrix::Diagonal(TR, 1)
  I_HJ   = purrr::map(1:J, ~ Matrix::Diagonal(H * .x, 1))
  I_G    = Matrix::Diagonal(G, 1)

  y_i = purrr::map(
    subjects,
    ~ MHPCA$data[[group]][MHPCA$data[[group]]$Subject == subjects[.x], ] %>%
      dplyr::pull(y_centered2)
  )

  j_i = purrr::map(
    subjects,
    ~ MHPCA$data[[group]][MHPCA$data[[group]]$Subject == subjects[.x], ] %>%
      dplyr::pull(Repetition) %>% unique()
  )

  length_j_i = purrr::map(
    subjects,
    ~ MHPCA$data[[group]][MHPCA$data[[group]]$Subject == subjects[.x], ] %>%
      dplyr::pull(Repetition) %>% unique() %>% length()
  )

  Phi1 = MHPCA$data[[group]] %>%
    dplyr::select(dplyr::starts_with("between")) %>%
    dplyr::slice(1:TR) %>%
    as.matrix(dimnames = list(names(.), NULL))

  z = Matrix::kronecker(rep(1, J), Phi1) %>%
    data.frame() %>%
    dplyr::mutate(Repetition = rep(obs, each = TR))

  Phi2 = MHPCA$data[[group]] %>%
    dplyr::select(dplyr::starts_with("within")) %>%
    dplyr::slice(1:TR) %>%
    as.matrix(dimnames = list(names(.), NULL))

  Phi2t = Matrix::t(Phi2)

  w = Matrix::kronecker(I_J[[J]], Phi2) %>%
    as.matrix() %>%
    data.frame() %>%
    dplyr::mutate(Repetition = rep(obs, each = TR))

  # storage ----
  sigma2     = rep(NA_real_, maxiter + 1)
  loglik_vec = rep(NA_real_, maxiter + 1)
  Lambda1    = vector(mode = "list", length = maxiter + 1)
  Lambda2    = vector(mode = "list", length = maxiter + 1)
  zeta_i     = vector(mode = "list", length = maxiter + 1)
  xi_i       = vector(mode = "list", length = maxiter + 1)

  # initialize variances ----
  Lambda1[[1]] = Matrix::Diagonal(G, 1)
  Lambda2[[1]] = Matrix::Diagonal(H, 1)
  sigma2[[1]]  = MHPCA$marg$within$functional[[group]]$sigma_2d

  E_mat = Lambda1[[1]]
  E_sqrt = sqrt(E_mat)

  # estimated level 1 covariance
  K_dB = Matrix::tcrossprod(Phi1, Phi1 %*% Lambda1[[1]])

  # Qinv = sigma2_d * I_TR + K_dW
  Q = 1 / sigma2[1] * I_TR - 1 / sigma2[1] ^ 2 * Phi2 %*% Matrix::solve(Matrix::solve(Lambda2[[1]]) + 1 / sigma2[1] * Matrix::crossprod(Phi2), Phi2t)

  # find inverse of JK_dB + Qinv
  Q_Phi1 = Q %*% Phi1

  omega_inv = vector("list", J)
  JK_Qinv_inv = vector("list", J)

  # calculate omega_inv for each number of visits
  for (j in 1:J) {
    # find inverse of CK_dB + F
    JK_Qinv_inv[[j]] = Q - j * Q_Phi1 %*% Matrix::solve(Matrix::solve(Lambda1[[1]]) + j * Matrix::crossprod(Phi1, Q_Phi1), Matrix::crossprod(Phi1, Q))
    omega_inv[[j]] = kronecker(I_J[[j]], Q) -
      kronecker(ones_J[[j]], Q %*% K_dB %*% JK_Qinv_inv[[j]])
  }

  z_i = vector("list", n); names(z_i) = subjects
  w_i = vector("list", n); names(w_i) = subjects
  ztz = vector("list", n); names(ztz) = subjects
  wtw = vector("list", n); names(wtw) = subjects
  wtz = vector("list", n); names(wtz) = subjects
  F_mat  = vector("list", J)
  F_sqrt = vector("list", J)

  omegainv_y    = vector(mode = "list", length = n); names(omegainv_y) = subjects
  yt_omegainv_y = vector(mode = "list", length = n); names(yt_omegainv_y) = subjects

  for (j in 1:J) {
    F_mat[[j]] = Matrix::kronecker(I_J[[j]], Lambda2[[1]])
    F_sqrt[[j]] = sqrt(F_mat[[j]])
  }

  for (i in subjects) {
    z_i[[i]] = dplyr::filter(z, Repetition %in% j_i[[i]]) %>%
      dplyr::select(-Repetition) %>% as.matrix() %>% Matrix::Matrix()
    w_i[[i]] = dplyr::filter(w, Repetition %in% j_i[[i]]) %>%
      dplyr::select(-Repetition) %>% as.matrix() %>% Matrix::Matrix()

    ztz[[i]] = length_j_i[[i]] * Matrix::crossprod(Phi1)
    wtw[[i]] = Matrix::kronecker(I_J[[length_j_i[[i]]]], Matrix::crossprod(Phi2))
    wtz[[i]] = Matrix::kronecker(one_J[[length_j_i[[i]]]], Matrix::crossprod(Phi2, Phi1))

    omegainv_y[[i]] = omega_inv[[length_j_i[[i]]]] %*% y_i[[i]]
    yt_omegainv_y[[i]] = Matrix::crossprod(y_i[[i]], omegainv_y[[i]])
  }

  loglik_vec[1] = loglik_calc_incomplete(
    sigma2[1], TR, n, z_i, w_i, E_sqrt, F_mat, F_sqrt, I_G, I_HJ,
    omegainv_y, yt_omegainv_y, wtw, ztz, wtz, subjects, length_j_i, J
  )

  # iterate for maxiter ----
  for (k in 1:maxiter) {

    # calculate scores
    zeta_i[[k]] = vector(mode = "list", length = n); names(zeta_i[[k]]) = subjects
    xi_i[[k]]   = vector(mode = "list", length = n); names(xi_i[[k]])   = subjects

    y_sum = 0
    trace_sum = 0

    for (i in subjects) {
      zeta_i[[k]][[i]] = Lambda1[[k]] %*% Matrix::crossprod(z_i[[i]], omegainv_y[[i]])
      rownames(zeta_i[[k]][[i]]) = colnames(Phi1)
      xi_i[[k]][[i]]   = Matrix::kronecker(I_J[[J]], Lambda2[[k]]) %*%
        Matrix::crossprod(w_i[[i]], omegainv_y[[i]])
      rownames(xi_i[[k]][[i]]) = rep(colnames(Phi2), each = J)
      y_sum     = y_sum + Matrix::crossprod(omegainv_y[[i]])
      trace_sum = trace_sum + sum(Matrix::diag(omega_inv[[length_j_i[[i]]]]))
    }

    sigma2[k + 1] = as.vector(sigma2[k] * sqrt(y_sum / trace_sum))

    Lambda1[[k + 1]] = Matrix::Diagonal(G, x = 0)
    Lambda2[[k + 1]] = Matrix::Diagonal(H, x = 0)

    for (j in 1:G) {
      gamma_sum = 0
      zeta_sum  = 0
      for (i in subjects) {
        A = Matrix::crossprod(z_i[[i]], omega_inv[[length_j_i[[i]]]] %*% z_i[[i]])
        gamma_sum = gamma_sum + A[j, j]
        zeta_sum = zeta_sum + zeta_i[[k]][[i]][j] ^ 2
      }
      Lambda1[[k + 1]][j, j] = sqrt(zeta_sum / gamma_sum)
    }


    for (l in 1:H) {
      gamma_sum = 0
      xi_sum    = 0
      for (i in 1:n) {
        B = Matrix::crossprod(w_i[[i]], omega_inv[[length_j_i[[i]]]] %*% w_i[[i]])
        for (j in 1:J) {
          index = l + (j - 1) * H
          gamma_sum = gamma_sum + B[index, index]
          xi_sum = xi_sum + xi_i[[k]][[i]][index] ^ 2
        }
      }
      Lambda2[[k + 1]][l, l] = sqrt(xi_sum / gamma_sum)
    }

    E_mat = Lambda1[[k + 1]]
    E_sqrt = sqrt(E_mat)

    # estimated level 1 covariance
    K_dB = Matrix::tcrossprod(Phi1, Phi1 %*% Lambda1[[k + 1]])

    # Qinv = sigma2_d * I_TR + K_dW
    Q = 1 / sigma2[k + 1] * I_TR - 1 / sigma2[k + 1] ^ 2 * Phi2 %*% Matrix::solve(Matrix::solve(Lambda2[[k + 1]]) + 1 / sigma2[k + 1] * Matrix::crossprod(Phi2), Phi2t)

    # find inverse of JK_dB + Qinv
    Q_Phi1 = Q %*% Phi1

    omega_inv = vector("list", J)
    JK_Qinv_inv = vector("list", J)

    # calculate omega_inv for each number of visits
    for (j in 1:J) {
      # find inverse of CK_dB + F
      JK_Qinv_inv[[j]] = Q - j * Q_Phi1 %*% Matrix::solve(Matrix::solve(Lambda1[[k + 1]]) + j * Matrix::crossprod(Phi1, Q_Phi1), Matrix::crossprod(Phi1, Q))
      omega_inv[[j]] = Matrix::kronecker(I_J[[j]], Q) -
        Matrix::kronecker(ones_J[[j]], Q %*% K_dB %*% JK_Qinv_inv[[j]])
    }

    for (j in 1:J) {
      F_mat[[j]] = Matrix::kronecker(I_J[[j]], Lambda2[[1]])
      F_sqrt[[j]] = sqrt(F_mat[[j]])
    }

    for (i in subjects) {
      omegainv_y[[i]] = omega_inv[[length_j_i[[i]]]] %*% y_i[[i]]
      yt_omegainv_y[[i]] = Matrix::crossprod(y_i[[i]], omegainv_y[[i]])
    }

    loglik_vec[k + 1] = loglik_calc_incomplete(
      sigma2[k + 1], TR, n, z_i, w_i, E_sqrt, F_mat, F_sqrt, I_G, I_HJ,
      omegainv_y, yt_omegainv_y, wtw, ztz, wtz, subjects, length_j_i, J
    )

    if (abs(loglik_vec[k + 1] - loglik_vec[k]) < epsilon * abs(loglik_vec[k + 1] + 1)) {
      break
    }
  }

  return(list(
    "level1_components" = names(dplyr::select(MHPCA$data[[group]], dplyr::starts_with("between"))),
    "level2_components" = names(dplyr::select(MHPCA$data[[group]], dplyr::starts_with("within"))),
    # "loglik_vec" = loglik_vec[!is.na(sigma2)],
    # "Lambda1"    = Lambda1[!sapply(Lambda1, is.null)],
    # "Lambda2"    = Lambda2[!sapply(Lambda2, is.null)],
    # "sigma2"     = sigma2[!is.na(sigma2)],
    # "zeta_i"     = zeta_i[!sapply(zeta_i, is.null)],
    # "xi_i"       = xi_i[!sapply(xi_i, is.null)]
    "loglik_vec" = loglik_vec[(k + 1)],
    "Lambda1"    = Lambda1[[(k + 1)]],
    "Lambda2"    = Lambda2[[(k + 1)]],
    "sigma2"     = sigma2[(k + 1)],
    "zeta_i"     = zeta_i[[k]],
    "xi_i"       = xi_i[[k]]
  ))
}

#' Log-likelihood calculation for MM algorithm with incomplete data
#'
#' @param sigma2 sigma2
#' @param TR TR
#' @param n n
#' @param z_i z_i
#' @param w_i w_i
#' @param E_sqrt E_sqrt
#' @param F_mat F_mat
#' @param F_sqrt F_sqrt
#' @param I_G I_G
#' @param I_HJ I_HJ
#' @param omegainv_y omegainv_y
#' @param yt_omegainv_y yt_omegainv_y
#' @param wtw wtw
#' @param ztz ztz
#' @param wtz wtz
#' @param subjects subjects
#' @param length_j_i length_j_i
#' @param J J
#'
loglik_calc_incomplete = function(
  sigma2, TR, n, z_i, w_i, E_sqrt, F_mat, F_sqrt, I_G, I_HJ,
  omegainv_y, yt_omegainv_y, wtw, ztz, wtz, subjects, length_j_i, J
) {
  middle_inv    = vector("list", n); names(middle_inv) = subjects
  logdet_level2 = vector("list", n); names(logdet_level2) = subjects
  logdet_level1 = vector("list", n); names(logdet_level1) = subjects
  F_inv = vector("list", n); names(F_inv) = subjects
  FWtWF = vector("list", n); names(FWtWF) = subjects
  EZtZE = vector("list", n); names(EZtZE) = subjects
  WtZE  = vector("list", n); names(WtZE)  = subjects

  F_inv = vector("list", J)
  for (j in 1:J) {
    F_inv[[j]] =  Matrix::solve(F_mat[[j]])
  }

  loglik_i = rep(NA_real_, n); names(loglik_i) = subjects
  for (i in subjects) {
    FWtWF[[i]] = F_sqrt[[length_j_i[[i]]]] %*% wtw[[i]] %*% F_sqrt[[length_j_i[[i]]]]
    EZtZE[[i]] = E_sqrt %*% ztz[[i]] %*% E_sqrt
    WtZE[[i]]  = as.matrix(wtz[[i]] %*% E_sqrt)

    middle_inv[[i]] = I_G + 1 / sigma2 * (
      EZtZE[[i]] - Matrix::crossprod(WtZE[[i]], Matrix::solve(sigma2 * F_inv[[length_j_i[[i]]]] + wtw[[i]], WtZE[[i]]))
    )

    logdet_level2_arg = I_HJ[[length_j_i[[i]]]] + 1 / sigma2 * FWtWF[[i]]
    logdet_level2[[i]] = Matrix::determinant(logdet_level2_arg, logarithm = TRUE)$modulus
    logdet_level1[[i]] = Matrix::determinant(middle_inv[[i]], logarithm = TRUE)$modulus

    loglik_i[[i]] = as.numeric(-1 / 2 * (
      TR * length_j_i[[i]] * log(sigma2) +
        logdet_level2[[i]] + logdet_level1[[i]] +
        yt_omegainv_y[[i]]
    ))
  }

  loglik = sum(loglik_i)

  return(loglik)
}
