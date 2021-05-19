#' M-HPCA Simulation
#' @description Function for simulating data as described in the supplementary
#' material of Campos et al. (202?)
#'
#' @param sig_eps error standard deviation (scalar)
#' @param n_d group sample size (scalar)
#' @param J number of repetitions per subject (scalar)
#' @param D number of groups (scalar)
#' @param num_reg number of regions (scalar)
#' @param num_time number of functional points (scalar)
#' @param missing_level allow for missing repetitions (logical)
#' @param K number of marginal level 1 eigenvectors (scalar)
#' @param L number of marginal level 1 eigenfunctions (scalar)
#' @param P number of marginal level 2 eigenvectors (scalar)
#' @param M number of marginal level 2 eigenfunctions (scalar)
#'
#' @return A list with the simulated components.
#' \itemize{
#'   \item mu - The simulated mean function
#'   \item eta - The simulated group-region-repetition shift from the mean
#'   \item v_k - The simulated level 1 eigenvectors
#'   \item phi_l - the simulated level 1 eigenfunctions
#'   \item lambda_kl - The simulated level 1 eigenvalues
#'   \item v_p - The simulated level 2 eigenvectors
#'   \item phi_m - The simulated level 2 eigenfunctions
#'   \item lambda_pm - The simulated level 2 eigenvalues
#'   \item xi - The simulated subject-specific scores
#'   \item zeta - The simulated subject/repetition-specific scores
#'   \item data - The simulated dataframe with columns for
#'   \item sig_eps - The simulated measurement error variance
#' }
#' @export
#'
#' @importFrom dplyr arrange mutate filter pull rename select
#' @importFrom purrr map pmap_dfr
#' @importFrom stats rnorm var
#' @importFrom stringr str_pad str_length
MHPCA_simulation <- function(
  sig_eps, n_d, J, D, num_reg, num_time, missing_level = FALSE, K = 2, L = 2,
  P = 2, M = 2
) {

  # to avoid issues with non-standard evaluation in tidyeval, set "global
  # variables" to NULL and remove them. this won't cause an issue with the rest
  # of the code.
  cond       <- NULL
  group      <- NULL
  k          <- NULL
  l          <- NULL
  p          <- NULL
  m          <- NULL
  subject    <- NULL
  score      <- NULL
  repetition <- NULL
  y_nonoise  <- NULL
  noise      <- NULL
  y_noise    <- NULL


  rm(list = c(
    "cond", "group", "k", "l", "p", "m", "subject", "score", "repetition",
    "y_nonoise", "noise", "y_noise"
  ))

  # 1. Define Model Components ---------------------------------------------

  # _a. Domains ----

  # Regional domain
  reg <- c(1:num_reg)

  # Functional domain
  func <- seq(0, 1, length.out = num_time)


  # _b. Global Variables ----

  # Number of eigencomponents
  # K <- 2 # level 1, regional
  # L <- 2 # level 1, functional
  # P <- 2 # level 2, regional
  # M <- 2 # level 2, functional

  data <- expand.grid(
    cond  = 1:J,
    group = 1:D,
    reg   = reg,
    func  = func
  ) %>%
    dplyr::arrange(cond, group, reg, func)


  # _c. Fixed Effects ----

  # overall mean function
  mu_vec <- 5 * (func - 0.5) * (func - 0.75)

  # group-region-condition-specific shifts
  eta_vec <- purrr::map(1:D, function(d) {
    purrr::map(1:J, function(j) {
      purrr::map(reg, ~ 2 * (-1) ^ (j) * (func - (d - 1)) ^ 3)
    })
  })


  # _d. Eigenfunctions ----

  # level 1 regional marginal eigenvectors
  v_k <- purrr::map(1:K, ~ sin(.x * (reg - 1) * pi / 8) / 2)

  # level 2 regional marginal eigenvectors
  v_p <- purrr::map(1:P, ~ cos(.x * (reg - 1) * pi / 8) / sqrt(5))

  # level 1 functional marginal eigenfunctions
  phi_l <- purrr::map(1:L, ~ sqrt(2) * sin(.x * pi * func))

  # level 2 functional marginal eigenfunctions
  phi_m <- purrr::map(1:M, ~ sqrt(2) * cos(.x * pi * func))

  # _e. Variance Components ----

  # variance for subject-specific scores
  lambda_kl <- expand.grid(k = 1:K, l = 1:L) %>%
    dplyr::mutate(lambda = 0.5 ^ {(k - 1) + K * (l - 1)})

  # variance for subject-repetition-specific scores
  lambda_pm <- expand.grid(p = 1:P, m = 1:M) %>%
    dplyr::mutate(lambda = 0.5 * 0.5^{(p - 1) + P * (m - 1)})


  # Simulate Data ----------------------------------------------------------

  # _a. Simulate subject-specific scores ----
  zeta <- expand.grid(
    k = 1:K,
    l = 1:L,
    i = 1:n_d,
    d = 1:D
  ) %>%
    purrr::pmap_dfr(function(d, i, k, l) {
      data.frame(
        group = d,
        subject = i,
        K = k,
        L = l,
        score = stats::rnorm(
          1, 0, sqrt(lambda_kl[lambda_kl$k == k & lambda_kl$l == l, "lambda"])
        )
      )
    })

  # _b. Simulate subject-repetition-specific scores ----
  xi <- expand.grid(
    p = 1:P,
    m = 1:M,
    i = 1:n_d,
    d = 1:D,
    j = 1:J
  ) %>%
    purrr::pmap_dfr(function(j, d, i, p, m) {
      data.frame(
        repetition = j,
        group = d,
        subject = i,
        P = p,
        M = m,
        score = rnorm(
          1, 0, sqrt(lambda_pm[lambda_pm$p == p & lambda_pm$m == m, "lambda"])
        )
      )
    })

  # _c. Simulate subject-specific trajectories ----
  data <- expand.grid(
    d = 1:D,
    j = 1:J,
    i = 1:n_d,
    r = reg
  ) %>%
    purrr::pmap_dfr(function(d, j, i, r) {
      # the mean trajectory for subjects in group d, repetition j, region r
      subj_trajectory <- mu_vec + eta_vec[[d]][[j]][[r]]

      # set up level 1 and level 2 vectors
      level_1 <- vector("numeric", length(func))
      level_2 <- vector("numeric", length(func))

      # add level 1 components
      for (k in 1:K) {
        for (l in 1:L) {
          zeta_di_kl <- zeta %>%
            dplyr::filter(group == d, subject == i, K == k, L == l) %>%
            dplyr::pull(score)
          level_1 <- level_1 + zeta_di_kl * v_k[[k]][[r]] * phi_l[[l]]
        }
      }

      # add level 2 components
      for (m in 1:M) {
        for (p in 1:P) {
          xi_cdi_mp <- xi %>%
            dplyr::filter(
              repetition == j,
              group == d,
              subject == i,
              M == m,
              P == p
            ) %>%
            dplyr::pull(score)
          level_2 <- level_2 + xi_cdi_mp * v_p[[p]][[r]] * phi_m[[m]]
        }
      }

      # create subject trajectory + random noise
      data <- data.frame(
        Observation = paste("Repetition", rep(j, length(func))),
        Group = paste("Group", rep(d, length(func))),
        Subject = paste(
          "Subject",
          stringr::str_pad(
            rep(i + n_d * (d - 1), length(func)),
            width = stringr::str_length(n_d * D),
            pad = "0",
            side = "left"
          )
        ),
        reg = paste0(
          "E",
          rep(stringr::str_pad(r, width = 2, pad = "0"), length(func))
        ),
        func = func,
        y_nonoise = subj_trajectory + level_1 + level_2,
        noise = stats::rnorm(length(func), 0, sig_eps),
        stringsAsFactors = FALSE
      ) %>%
        dplyr::mutate(y_noise = y_nonoise + noise)

      SNR <- stats::var(data$y_nonoise) / sig_eps^2

      dplyr::mutate(data, SNR = SNR) %>%
        dplyr::rename(y = y_noise)
    })

  obs <- unique(data$Observation)

  SNR <- mean(data$SNR)
  data <- dplyr::select(data, -SNR)

  if (missing_level) {
    level_logical <- matrix(1, nrow = D * n_d, ncol = J)
    rownames(level_logical) <- unique(data$Subject)
    level_num <- sample(1:J, D * n_d, replace = TRUE, prob = c(0.3, 0.7))

    for (i in 1:length(level_num)) {
      if (level_num[i] < J) {
        level_logical[i, ] <- sample(
          c(rep(0, level_num[i]), rep(1, J - level_num[i]))
        )
      } else {
        level_logical[i, ] <- rep(0, J)
      }
    }

    for (i in unique(data$Subject)) {
      if (sum(level_logical[i, ]) > 0) {
        data = data[
          -which(
            data$Subject == i &
              data$Observation %in% obs[as.logical(level_logical[i, ])]),
        ]
      }
    }
  }


  return(list(
    "mu"        = mu_vec,
    "eta"       = eta_vec,
    "v_k"       = v_k,
    "phi_l"     = phi_l,
    "lambda_kl" = lambda_kl,
    "v_p"       = v_p,
    "phi_m"     = phi_m,
    "lambda_pm" = lambda_pm,
    "xi"        = xi,
    "zeta"      = zeta,
    "data"      = data,
    "sig_eps"   = sig_eps ^ 2
  ))
}


#' @title MHPCA Simulation for Group-Level Inference
#' @description Function for simulating data for group-level inference as
#' described in the supplementary materials.
#'
#' @param sig_eps error standard deviation (scalar)
#' @param n_d group sample size (scalar)
#' @param J number of repetitions per subject (scalar)
#' @param D number of groups (scalar)
#' @param num_reg number of regions (scalar)
#' @param num_time number of functional points (scalar)
#' @param missing_level allow for missing repetitions (logical)
#' @param test_group test group (scalar)
#' @param test_region test region (scalar)
#' @param delta tuning parameter (scalar)
#'
#' @return A list with the simulated components.
#' \itemize{
#'   \item mu - The simulated mean function
#'   \item eta - The simulated group-region-repetition shift from the mean
#'   \item v_k - The simulated level 1 eigenvectors
#'   \item phi_l - the simulated level 1 eigenfunctions
#'   \item lambda_kl - The simulated level 1 eigenvalues
#'   \item v_p - The simulated level 2 eigenvectors
#'   \item phi_m - The simulated level 2 eigenfunctions
#'   \item lambda_pm - The simulated level 2 eigenvalues
#'   \item xi - The simulated subject-specific scores
#'   \item zeta - The simulated subject/repetition-specific scores
#'   \item data - The simulated dataframe with columns for
#'   \item sig_eps - The simulated measurement error variance
#' }
#'
#' @importFrom dplyr arrange mutate filter pull rename select
#' @importFrom purrr map pmap_dfr
#' @importFrom stats rnorm var
#' @importFrom stringr str_pad str_length
MHPCA_simulation_within_group_test <- function(
  sig_eps, n_d, J, D, num_reg, num_time, missing_level = FALSE,
  test_group, test_region, delta
) {

  # to avoid issues with non-standard evaluation in tidyeval, set "global
  # variables" to NULL and remove them. this won't cause an issue with the rest
  # of the code.
  cond       <- NULL
  group      <- NULL
  k          <- NULL
  l          <- NULL
  p          <- NULL
  m          <- NULL
  subject    <- NULL
  score      <- NULL
  repetition <- NULL
  y_nonoise  <- NULL
  noise      <- NULL
  y_noise    <- NULL


  rm(list = c(
    "cond", "group", "k", "l", "p", "m", "subject", "score", "repetition",
    "y_nonoise", "noise", "y_noise"
  ))


  # 1. Define Model Components ---------------------------------------------

  # _a. Domain ----

  # Regional domain
  reg <- c(1:num_reg)

  # Functional domain
  func <- seq(0, 1, length.out = num_time)


  # _b. Global Variables ----

  # Number of eigencomponents
  K <- 2 # level 1, regional
  L <- 2 # level 1, functional
  P <- 2 # level 2, regional
  M <- 2 # level 2, functional

  data <- expand.grid(
    cond = 1:J,
    group = 1:D,
    reg = reg,
    func = func
  ) %>%
    dplyr::arrange(cond, group, reg, func)


  # _c. Fixed Effects ----

  # overall mean function
  mu_vec <- 5 * (func - 0.5) * (func - 0.75)

  # group-region-condition-specific shifts
  eta_vec <- purrr::map(1:D, function(d) {
    if (d == test_group) {
      purrr::map(1:J, function(j) {
        purrr::map(reg, function(r) {
          if (r == test_region) {
            rep((-1) ^ j * delta, num_time)
          } else {
            rep(0, num_time)
          }
        })
      })
    } else {
      purrr::map(1:J, function(j) {
        purrr::map(reg, ~ rep(0, num_time))
      })
    }
  })


  # _d. Eigenfunctions ----

  # level 1 regional marginal eigenvectors
  v_k <- purrr::map(1:K, ~ sin(.x * (reg - 1) * pi / 8) / 2)

  # level 2 regional marginal eigenvectors
  v_p <- purrr::map(1:P, ~ cos(.x * (reg - 1) * pi / 8) / sqrt(5))

  # level 1 functional marginal eigenfunctions
  phi_l <- purrr::map(1:L, ~ sqrt(2) * sin(.x * pi * func))

  # level 2 functional marginal eigenfunctions
  phi_m <- purrr::map(1:M, ~ sqrt(2) * cos(.x * pi * func))

  # _e. Variance Components ----

  # variance for subject-specific scores
  lambda_kl <- expand.grid(k = 1:K, l = 1:L) %>%
    dplyr::mutate(lambda = 0.5 ^ {(k - 1) + K * (l - 1)})

  # variance for subject-repetition-specific scores
  lambda_pm <- expand.grid(p = 1:P, m = 1:M) %>%
    dplyr::mutate(lambda = 0.5 * 0.5^{(p - 1) + P * (m - 1)})


  # Simulate Data ----------------------------------------------------------

  # _a. Simulate subject-specific scores ----
  zeta <- expand.grid(
    k = 1:K,
    l = 1:L,
    i = 1:n_d,
    d = 1:D
  ) %>%
    purrr::pmap_dfr(function(d, i, k, l) {
      data.frame(
        group = d,
        subject = i,
        K = k,
        L = l,
        score = stats::rnorm(
          1, 0, sqrt(lambda_kl[lambda_kl$k == k & lambda_kl$l == l, "lambda"])
        )
      )
    })

  # _b. Simulate subject-repetition-specific scores ----
  xi <- expand.grid(
    p = 1:P,
    m = 1:M,
    i = 1:n_d,
    d = 1:D,
    j = 1:J
  ) %>%
    purrr::pmap_dfr(function(j, d, i, p, m) {
      data.frame(
        repetition = j,
        group = d,
        subject = i,
        P = p,
        M = m,
        score = stats::rnorm(
          1, 0, sqrt(lambda_pm[lambda_pm$p == p & lambda_pm$m == m, "lambda"])
        )
      )
    })

  # _c. Simulate subject-specific trajectories ----
  data <- expand.grid(
    d = 1:D,
    j = 1:J,
    i = 1:n_d,
    r = reg
  ) %>%
    purrr::pmap_dfr(function(d, j, i, r) {
      # the mean trajectory for subjects in group d, repetition j, region r
      subj_trajectory <- mu_vec + eta_vec[[d]][[j]][[r]]

      # set up level 1 and level 2 vectors
      level_1 <- vector("numeric", length(func))
      level_2 <- vector("numeric", length(func))

      # add level 1 components
      for (k in 1:K) {
        for (l in 1:L) {
          zeta_di_kl <- zeta %>%
            filter(group == d, subject == i, K == k, L == l) %>%
            dplyr::filter(score)
          level_1 <- level_1 + zeta_di_kl * v_k[[k]][[r]] * phi_l[[l]]
        }
      }

      # add level 2 components
      for (m in 1:M) {
        for (p in 1:P) {
          xi_cdi_mp <- xi %>%
            dplyr::filter(
              repetition == j,
              group == d,
              subject == i,
              M == m,
              P == p
            ) %>%
            dplyr::pull(score)
          level_2 <- level_2 + xi_cdi_mp * v_p[[p]][[r]] * phi_m[[m]]
        }
      }

      # create subject trajectory + random noise
      data <- data.frame(
        Observation = paste("Repetition", rep(j, length(func))),
        Group = paste("Group", rep(d, length(func))),
        Subject = paste(
          "Subject",
          stringr::str_pad(
            rep(i + n_d * (d - 1), length(func)),
            width = stringr::str_length(n_d * D),
            pad = "0",
            side = "left"
          )
        ),
        reg = paste0(
          "E",
          rep(stringr::str_pad(r, width = 2, pad = "0"), length(func))
        ),
        func = func,
        y_nonoise = subj_trajectory + level_1 + level_2,
        noise = stats::rnorm(length(func), 0, sig_eps),
        stringsAsFactors = FALSE
      ) %>%
        dplyr::mutate(y_noise = y_nonoise + noise)

      SNR <- stats::var(data$y_nonoise) / sig_eps^2

      dplyr::mutate(data, SNR = SNR) %>%
        dplyr::rename(y = y_noise)
    })

  obs <- unique(data$Observation)

  SNR <- mean(data$SNR)
  data <- dplyr::select(data, -SNR)

  if (missing_level) {
    level_logical <- matrix(1, nrow = D * n_d, ncol = J)
    rownames(level_logical) <- unique(data$Subject)
    level_num <- sample(1:J, D * n_d, replace = TRUE, prob = c(0.3, 0.7))

    for (i in 1:length(level_num)) {
      if (level_num[i] < J) {
        level_logical[i, ] <- sample(
          c(rep(0, level_num[i]), rep(1, J - level_num[i]))
        )
      } else {
        level_logical[i, ] <- rep(0, J)
      }
    }

    for (i in unique(data$Subject)) {
      if (sum(level_logical[i, ]) > 0) {
        data = data[-which(
          data$Subject == i &
            data$Observation %in% obs[as.logical(level_logical[i, ])]
        ), ]
      }
    }
  }


  return(list(
    "mu"        = mu_vec,
    "eta"       = eta_vec,
    "v_k"       = v_k,
    "phi_l"     = phi_l,
    "lambda_kl" = lambda_kl,
    "v_p"       = v_p,
    "phi_m"     = phi_m,
    "lambda_pm" = lambda_pm,
    "xi"        = xi,
    "zeta"      = zeta,
    "data"      = data,
    "sig_eps"   = sig_eps ^ 2,
    "SNR"       = SNR
  ))
}
