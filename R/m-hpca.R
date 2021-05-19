#' @title Multilevel Hybrid Principal Component Analysis
#' @description Function for performing MHPCA decomposition described in
#'  "Multilevel Hybrid Principal Components Analysis For Region-Referenced
#'  Functional EEG Data" by Campos et al. (202?), including estimation of
#'  fixed effects and marginal covariance functions, marginal eigencompoents,
#'  subject-specific scores, variance components, and measurement error
#'  variance.
#' @param data dataframe in long format with six labeled columns
#'  (Observation: (character vector),
#'  Subject: subject IDs (character vector),
#'  Group: subject group (character vector),
#'  func: functional argument (numeric vector),
#'  reg: regional argument (character vector),
#'  y: region-referenced functional data (numeric vector))
#'  and row length equal to the length of the vectorized region-referenced
#'  observations across all subjects and groups
#' @param maxiter maximum number of iterations for MM algorithm (scalar)
#' @param epsilon epsilon value for determining log-likelihood convergence (scalar)
#' @param fve_cutoff fraction of variance cutoff for reducing the number of product components used in the mixed effects model (scalar in (0, 1))
#' @param nknots number of knots to use for smoothing splines
#' @param reduce should the number of product components be reduced and the mixed effects model re-estimated (logical)
#' @param quiet display messages for timing (logical)
#'
#' @return A list with
#' \itemize{
#'  \item mu: the overall mean function (vector),
#'  \item eta: group-region-level-specific shifts (dataframe),
#'  \item covar: list with total, between and within covariance matrices,
#'  \item marg: list with between and within marginal covariances,
#'  \item model: list of models for each group, including scores and variances,
#'  \item data: data with all of the different pieces from the estimation,
#'  \item FVE: fraction of variance explained.
#' }
#'
#' @export
#'
#' @import data.table
#' @importFrom dplyr mutate select all_of pull group_by ungroup filter full_join summarise inner_join starts_with
#' @importFrom Matrix rowMeans sparseMatrix crossprod
#' @importFrom mgcv gam te
#' @importFrom purrr map_dfc pmap_dfr map map_dfr
#' @importFrom stats smooth.spline predict complete.cases
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom tictoc tic toc
#' @importFrom tidyr spread
#' @importFrom data.table data.table
#' @importFrom dplyr filter mutate select all_of pull group_by ungroup full_join summarise inner_join starts_with
#' @importFrom fANCOVA loess.as
#' @importFrom Matrix rowMeans sparseMatrix crossprod
#' @importFrom mgcv gam te
#' @importFrom pracma trapz
#' @importFrom purrr map_dfc pmap_dfr map map_dfr
#' @importFrom stats smooth.spline predict complete.cases
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom tictoc tic toc
#' @importFrom tidyr spread
#' @importFrom tidyselect peek_vars
MHPCA_decomp <- function(
  data, maxiter, epsilon, fve_cutoff, nknots, reduce = TRUE, quiet = FALSE
) {

  # to avoid issues with non-standard evaluation in tidyeval, set "global
  # variables" to NULL and remove them. this won't cause an issue with the rest
  # of the code.
  Observation  <- NULL
  Group        <- NULL
  Subject      <- NULL
  reg          <- NULL
  func         <- NULL
  func_ind     <- NULL
  reg_ind      <- NULL
  index        <- NULL
  overall      <- NULL
  Overall_Mean <- NULL
  y_centered   <- NULL
  eta          <- NULL
  row_ind      <- NULL
  rt_ind       <- NULL
  y_centered2  <- NULL
  r            <- NULL
  r_prime      <- NULL
  t_prime      <- NULL
  covar        <- NULL
  y            <- NULL

  rm(list = c(
    "Observation", "Group", "Subject", "reg", "func", "func_ind", "reg_ind",
    "index", "overall", "Overall_Mean", "y_centered", "eta", "row_ind",
    "rt_ind", "y_centered2", "r", "r_prime", "t_prime", "covar", "y"
  ))

  tictoc::tic("MHPCA Decomposition and Estimation")


  # 0. Format data and create return list ----------------------------------

  tictoc::tic("0. Data formatting")

  # convert to data.table and arrange
  data <- data.table::data.table(data)
  data <- data[order(Observation, Group, Subject, reg, func)]


  # define global variables
  id                <- unique(data$Subject)    # subject IDs
  groups            <- unique(data$Group)      # groups
  names(groups)     <- groups                  # name groups
  d_tot             <- length(groups)          # total number of groups
  regional          <- unique(data$reg)
  names(regional)   <- regional
  functional        <- unique(data$func)
  f_tot             <- length(functional)      # total grid pts functional dom
  r_tot             <- length(regional)        # total grid pts regional dom
  obs               <- unique(data$Observation)  # obs
  names(obs)        <- obs              # name obs
  j_tot             <- length(obs)      # total number of obs

  # create dimension indices for merging
  data[, func_ind := (match(func, functional))]    # functional index
  data[, reg_ind  := (match(reg, regional))]       # regional index
  # functional, regional index
  data[, index    := ((reg_ind - 1) + r_tot * (func_ind - 1) + 1)]

  # create list for output
  MHPCA <- matrix(list(), 7, 1)
  names(MHPCA) <- c("mu", "eta", "covar", "marg", "model", "data", "FVE")

  tictoc::toc(quiet = quiet)


  # 1. Estimation of Fixed Effects -----------------------------------------

  tictoc::tic("1. Estimation of Fixed Effects")

  # _a. calculate overall mean function ----

  # calculate overall mean function by fitting a spline in each group then
  # calculating the point-wise average of those group means
  MHPCA$mu <- purrr::map_dfc(
    groups,
    function(d) {
      dat <- dplyr::filter(data, Group == d)
      spline_fit <- stats::smooth.spline(x = dat$func, y = dat$y, nknots = nknots)
      group_mean <- stats::predict(spline_fit, functional)$y
    }
  ) %>%
    dplyr::mutate(overall = Matrix::rowMeans(dplyr::select(., dplyr::all_of(groups)))) %>%
    dplyr::pull(overall)

  # center data by overall mean function
  data <- data %>%
    dplyr::group_by(Subject, Group, reg, Observation) %>%
    dplyr::mutate(
      Overall_Mean = MHPCA[["mu"]],
      y_centered = y - Overall_Mean
    ) %>%
    dplyr::ungroup()

  # _b. calculate group-region-level-specific shifts ----

  # obtain group-condition-region-specific shifts by fitting a spline on the
  # centered data within each group-region-level
  MHPCA$eta <- expand.grid(
    d = groups,
    j = obs,
    r = regional,
    stringsAsFactors = FALSE
  ) %>% purrr::pmap_dfr(
    function(d, j, r) {
      dat <- dplyr::filter(data, Group == d, reg == r, Observation == j)
      spline_fit <- stats::smooth.spline(dat$func, dat$y_centered, nknots = nknots)
      data.frame(
        Group = d,
        reg = r,
        Observation = j,
        func = functional,
        eta = stats::predict(spline_fit, functional)$y,
        stringsAsFactors = FALSE
      )
    })

  # center data by group-region-level-specific shifts
  data <- dplyr::full_join(
    data,
    MHPCA$eta,
    by = c("Observation", "Group", "reg", "func")
  ) %>%
    dplyr::mutate(y_centered2 = y_centered - eta)

  tictoc::toc(quiet = quiet)


  # 2. Estimation of Covariances -------------------------------------------

  if (quiet == FALSE) {
    cat("2. Estimation of Covariances\n")
  }

  tictoc::tic("    a. Raw Covariances")
  MHPCA$covar <- matrix(list(), 3, 1)
  names(MHPCA$covar) <- c("total", "between", "within")

  # _a. raw covariances ----
  MHPCA$covar$total <- purrr::map(
    groups,
    function(d) {
      # subset data by group
      data = data.table::data.table(data[which(data$Group == d), ])

      # obtain unique group IDs
      group_subj = as.matrix(unique(data$Subject))

      # number of subjects in group
      lid = length(group_subj)

      # assign an index to the columns of the marginal design matrix
      data[, row_ind := ((match(Subject, group_subj) - 1) +
                           lid * (match(Observation, obs) - 1) + 1)]
      data[, rt_ind := ((func_ind - 1) + f_tot * (reg_ind - 1) + 1)]

      # form sparse marginal design matrix
      sparse_mat = Matrix::sparseMatrix(
        i = data[, row_ind],
        j = data[, rt_ind],
        x = data[, y_centered2],
        dims = c(lid * j_tot, r_tot * f_tot)
      )

      # form sparse marginal indicator matrix
      sparse_mat_ind = Matrix::sparseMatrix(
        i = data[, row_ind],
        j = data[, rt_ind],
        x = 1,
        dims = c(lid * j_tot, r_tot * f_tot)
      )

      # calculate covariance matrix
      as.matrix(Matrix::crossprod(sparse_mat, sparse_mat)
                / Matrix::crossprod(sparse_mat_ind, sparse_mat_ind))
    }
  )

  MHPCA$covar$between <- purrr::map(
    groups,
    function(d) {
      # subset data by group
      data = data.table::data.table(data[which(data$Group == d), ])

      # obtain unique group IDs
      group_subj = as.matrix(unique(data$Subject))

      # number of subjects in group
      lid = length(group_subj)

      # assign an index to the columns of the marginal design matrix
      data[, row_ind := (match(Subject, group_subj))]
      data[, rt_ind := ((func_ind - 1) + f_tot * (reg_ind - 1) + 1)]

      sparse_mats <- vector("list", j_tot)
      sparse_mat_inds <- vector("list", j_tot)

      for (j1 in 1:j_tot) {
        sparse_mats[[j1]] = Matrix::sparseMatrix(
          i = data[which(data$Observation == obs[j1]), row_ind],
          j = data[which(data$Observation == obs[j1]), rt_ind],
          x = data[which(data$Observation == obs[j1]), y_centered2],
          dims = c(lid, r_tot * f_tot)
        )

        sparse_mat_inds[[j1]] = Matrix::sparseMatrix(
          i = data[which(data$Observation == obs[j1]), row_ind],
          j = data[which(data$Observation == obs[j1]), rt_ind],
          x = 1,
          dims = c(lid, r_tot * f_tot)
        )
      }

      big_sum = 0
      for (j1 in 1:j_tot) {
        for (j2 in 1:j_tot) {
          if (j1 < j2) {
            big_sum = big_sum +
              as.matrix(
                Matrix::crossprod(sparse_mats[[j1]], sparse_mats[[j2]]) /
                  Matrix::crossprod(sparse_mat_inds[[j1]], sparse_mat_inds[[j2]])
              )
          }}}
      big_sum
    }
  )

  MHPCA$covar$within <- purrr::map(
    groups,
    function(d) MHPCA$covar$total[[d]] - MHPCA$covar$between[[d]]
  )

  tictoc::toc(quiet = quiet)

  # _b,c. marginal covariances and smoothing ----
  marginal_covar <- function(
    which.cov, # which level of covariance (total, within, between)
    smooth,    # 1 to turn on covariance smoothing, 0 ow
    d          # group indicator
  ) {

    cov_mat <- expand.grid(
      t = functional,
      r = regional,
      t_prime = functional,
      r_prime = regional
    ) %>%
      dplyr::select(r, t, r_prime, t_prime) %>%
      dplyr::mutate(covar = as.vector(MHPCA$covar[[which.cov]][[d]]))

    if (smooth == TRUE) {

      cov_mat <- cov_mat %>%
        dplyr::filter(r == r_prime) %>%
        dplyr::group_by(t, t_prime) %>%
        dplyr::summarise(covar = mean(covar), .groups = "drop_last") %>%
        dplyr::ungroup()

      # vectorize covariance for smoothing
      cov_vec_d <- cov_mat %>%
        dplyr::filter(stats::complete.cases(.))

      if (which.cov == "within") {
        cov_vec_d <- dplyr::filter(cov_vec_d, t != t_prime)
      }

      # smooth the pooled sample covariances
      cov_vec_s <- mgcv::gam(
        covar ~ mgcv::te(t, t_prime, k = nknots, bs = "ps"),
        data = cov_vec_d
      )

      # form 2D grid to predict covariance function
      newdata <- dplyr::select(cov_mat, t, t_prime)

      # estimated marginal covariance function
      cov_mat_s <- matrix(stats::predict(cov_vec_s, newdata = newdata), nrow = f_tot)

      # symmetrize covariance function
      cov_mat_s <- (cov_mat_s + t(cov_mat_s)) / 2

      if (which.cov == "within") {
        # calculate measurement error variance
        # extract diagonal entries from pooled sample covariance
        marg_diag <- cov_mat[which(cov_mat$t == cov_mat$t_prime), ] %>%
          as.matrix()

        loess_diag <- suppressWarnings(fANCOVA::loess.as(
          marg_diag[, "t"],
          marg_diag[, "covar"],
          degree = 1,
          criterion = "gcv",
          user.span = NULL,
          plot = FALSE
        ))

        sigma_2 <- abs(mean(stats::predict(loess_diag, functional) - diag(cov_mat_s)))

        return(list(
          "sigma_hat" = cov_mat,
          "sigma_tilde" = cov_mat_s,
          "sigma_2d" = sigma_2
        ))
      } else {
        return(list(
          "sigma_hat" = cov_mat,
          "sigma_tilde" = cov_mat_s
        ))
      }

    } else {

      cov_mat <- cov_mat %>%
        dplyr::filter(t == t_prime) %>%
        dplyr::group_by(r, r_prime) %>%
        dplyr::summarise(covar = mean(covar), .groups = "drop_last") %>%
        dplyr::ungroup() %>%
        tidyr::spread(r_prime, covar) %>%
        tibble::column_to_rownames("r") %>%
        as.matrix()

      cov_mat_s <- cov_mat

      if (which.cov == "within") {
        cov_mat_s <- cov_mat_s -
          diag(MHPCA$marg$within$functional[[d]]$sigma_2d, r_tot)
      }

      cov_mat_s <- (cov_mat_s + t(cov_mat_s)) / 2

      return(list(
        "sigma_hat" = cov_mat,
        "sigma_tilde" = cov_mat_s
      ))
    }
  }

  tictoc::tic("    b,c. Estimation of Marginal Covariances and Smoothing")
  for (d in groups) {
    MHPCA$marg$between$functional[[d]] <- matrix(list(), 2, 1)
    names(MHPCA$marg$between$functional[[d]]) <-
      c("covariance", "eigendecomp")
    MHPCA$marg$between$regional[[d]] <- matrix(list(), 2, 1)
    names(MHPCA$marg$between$regional[[d]]) <-
      c("covariance", "eigendecomp")
    MHPCA$marg$within$functional[[d]] <- matrix(list(), 2, 1)
    names(MHPCA$marg$within$functional[[d]]) <-
      c("covariance", "eigendecomp")
    MHPCA$marg$within$regional[[d]] <- matrix(list(), 2, 1)
    names(MHPCA$marg$within$regional[[d]]) <-
      c("covariance", "eigendecomp")

    MHPCA$marg$between$functional[[d]] <- marginal_covar("between", TRUE, d)
    MHPCA$marg$between$regional[[d]]   <- marginal_covar("between", FALSE, d)
    MHPCA$marg$within$functional[[d]]  <- marginal_covar("within", TRUE, d)
    MHPCA$marg$within$regional[[d]]    <- marginal_covar("within", FALSE, d)
  }
  tictoc::toc(quiet = quiet)


  # 3. Estimation of Marginal Eigencomponents ------------------------------

  tictoc::tic("3. Estimation of Marginal Eigencomponents")

  eigenfun <- function(covariance, marginal_domain, covariance_function) {
    # compute eigendecomp
    eigen_temp <- eigen(covariance, symmetric = TRUE)

    # obtain positive eigenvalues
    eigen_temp$values <- eigen_temp$values[which(eigen_temp$values > 0)]

    # obtain eigenvectors associated with positive eigenvalues
    eigen_temp$vectors <- as.matrix(
      eigen_temp$vectors[, 1:length(eigen_temp$values)]
    )

    if (covariance_function == TRUE) {
      for (r in 1:length(eigen_temp$values)) {
        # normalize the eigenfunctions over the domain
        eigen_temp$vectors[, r] <- eigen_temp$vectors[, r] /
          sqrt(pracma::trapz(marginal_domain, eigen_temp$vectors[, r] ^ 2))
      }
    } else {
      rownames(eigen_temp$vectors) <- colnames(covariance)
    }

    # calculate number of components to use
    K <- which.max(cumsum(eigen_temp$values)/sum(eigen_temp$values) > 0.99)

    list(
      "values"  = eigen_temp$values,
      "vectors" = eigen_temp$vectors,
      "FVE"     = K
    )
  }

  for (d in groups) {
    # employ FPCA on functional between marginal covariance
    MHPCA$marg$between$functional[[d]]$eigendecomp <-
      eigenfun(MHPCA$marg$between$functional[[d]]$sigma_tilde, functional, TRUE)
    # employ PCA on regional between marginal covariance
    MHPCA$marg$between$regional[[d]]$eigendecomp <-
      eigenfun(MHPCA$marg$between$regional[[d]]$sigma_tilde, regional, FALSE)

    # employ FPCA on functional within marginal covariance
    MHPCA$marg$within$functional[[d]]$eigendecomp <-
      eigenfun(MHPCA$marg$within$functional[[d]]$sigma_tilde, functional, TRUE)
    # employ PCA on regional within marginal covariance
    MHPCA$marg$within$regional[[d]]$eigendecomp <-
      eigenfun(MHPCA$marg$within$regional[[d]]$sigma_tilde, regional, FALSE)
  }

  tictoc::toc(quiet = quiet)


  # 4. Estimation of Variance Components -----------------------------------

  if (quiet == FALSE) {
    cat("4. Estimation of Variance Components\n")
  }

  tictoc::tic("    a. Fit big model")
  # form vectorized versions of the multidimensional orthonormal basis
  prod_surf <- function(level, d) {
    list_surf <- data.frame(
      reg = rep(regional, each = f_tot),
      func = rep(functional, times = r_tot),
      stringsAsFactors = FALSE
    )

    reg_FVE <- MHPCA$marg[[level]]$regional[[d]]$eigendecomp$FVE
    fun_FVE <- MHPCA$marg[[level]]$functional[[d]]$eigendecomp$FVE

    for (x in 1:reg_FVE) {
      for (y in 1:fun_FVE) {
        # xth region marginal eigenvector
        v_x <- MHPCA$marg[[level]]$regional[[d]]$eigendecomp$vectors[, x]

        # yth functional marginal eigenfunction
        phi_y <- MHPCA$marg[[level]]$functional[[d]]$eigendecomp$vectors[, y]

        varname <- paste0(level, "_", x, "_", y)

        # multidimensional orthonormal basis
        list_surf <- dplyr::mutate(
          list_surf,
          !!varname := v_x[match(reg, regional)] *
            phi_y[match(func, functional)]
        ) %>%
          dplyr::select(reg, func, sort(tidyselect::peek_vars()))
      }
    }
    return(list_surf)
  }

  between_prod_surf <- purrr::map(groups, ~ prod_surf("between", .x))
  within_prod_surf  <- purrr::map(groups, ~ prod_surf("within", .x))

  # create group-specific design matrix
  MHPCA$data <- purrr::map(
    groups,
    function(d) {
      dplyr::filter(data, Group == d) %>%
        dplyr::inner_join(between_prod_surf[[d]], by = c("reg", "func")) %>%
        dplyr::inner_join(within_prod_surf[[d]], by = c("reg", "func"))
    }
  )

  # _a. calculate variance components and blups ----
  MHPCA$model <- matrix(list(), 2, 1)
  names(MHPCA$model) <- c("full", "final")
  # calculate variance components and measurement error variance by mm alg
  MHPCA$model$full <- purrr::map(groups, ~ mm(MHPCA, .x, maxiter, epsilon, j_tot))
  # MHPCA$model$full <- map(groups, ~ mm_complete(MHPCA, .x, maxiter, epsilon))

  tictoc::toc(quiet = quiet)

  if (reduce == TRUE) {
    # _b. choosing number of components ----
    # using the FVE_dB and FVE_dW described in the paper, choose the appropriate
    # number of components to keep at each level

    tictoc::tic("    b. Choose number of components")
    MHPCA[["FVE"]] <- vector("list", 6)
    names(MHPCA[["FVE"]]) <- c(
      "D_between","D_within", "FVE_dG", "G_prime", "FVE_dH", "H_prime"
    )

    # keeping track of the number of iterations from the mm alg
    # iters <- map_dbl(
    #   groups,
    #   ~ length(MHPCA$model$full[[.x]]$Lambda1)
    # )

    MHPCA$FVE$D_between <- purrr::map_dfr(
      groups,
      ~ data.frame(
        # D_between = sum(diag(MHPCA$model$full[[.x]]$Lambda1[[iters[.x]]]))
        D_between = sum(diag(MHPCA$model$full[[.x]]$Lambda1))
      ),
      .id = "Group"
    )
    MHPCA$FVE$D_within <- purrr::map_dfr(
      groups,
      # ~ data.frame(D_within = sum(diag(MHPCA$model$full[[.x]]$Lambda2[[iters[.x]]]))),
      ~ data.frame(D_within = sum(diag(MHPCA$model$full[[.x]]$Lambda2))),
      .id = "Group"
    )
    MHPCA$FVE$FVE_dG <- purrr::map(
      groups,
      ~ data.frame(
        # FVE = cumsum(sort(diag(MHPCA$model$full[[.x]]$Lambda1[[iters[.x]]]),
        FVE = cumsum(sort(diag(MHPCA$model$full[[.x]]$Lambda1),
                          decreasing = TRUE)) /
          dplyr::pull(dplyr::filter(MHPCA$FVE$D_between, Group == .x), D_between)
      ) %>%
        tibble::rownames_to_column("term") %>%
        dplyr::mutate(term = as.numeric(term)),
      .id = "Group"
    )
    MHPCA$FVE$FVE_dH <- purrr::map(
      groups,
      ~ data.frame(
        # FVE = cumsum(sort(diag(MHPCA$model$full[[.x]]$Lambda2[[iters[.x]]]),
        FVE = cumsum(sort(diag(MHPCA$model$full[[.x]]$Lambda2),
                          decreasing = TRUE)) /
          dplyr::pull(dplyr::filter(MHPCA$FVE$D_within, Group == .x), D_within)
      ) %>%
        tibble::rownames_to_column("term") %>%
        dplyr::mutate(term = as.numeric(term)),
      .id = "Group"
    )

    MHPCA$FVE$G_prime <- purrr::map(
      groups,
      ~  ifelse(
        nrow(MHPCA$FVE$FVE_dG[[.x]][MHPCA$FVE$FVE_dG[[.x]]$FVE > fve_cutoff, ]) == 0,
        max(MHPCA$FVE$FVE_dG[[.x]][, "term"]),
        MHPCA$FVE$FVE_dG[[.x]][which.max(MHPCA$FVE$FVE_dG[[.x]]$FVE > fve_cutoff), "term"]
      )
    )

    MHPCA$FVE$H_prime <- purrr::map(
      groups,
      ~ ifelse(
        nrow(MHPCA$FVE$FVE_dH[[.x]][MHPCA$FVE$FVE_dH[[.x]]$FVE > fve_cutoff, ]) == 0,
        max(MHPCA$FVE$FVE_dH[[.x]][, "term"]),
        MHPCA$FVE$FVE_dH[[.x]][which.max(MHPCA$FVE$FVE_dH[[.x]]$FVE > fve_cutoff), "term"]
      )
    )

    keep_between <- purrr::map(
      groups,
      function(d) {
        # lambdas <- diag(MHPCA$model$full[[d]]$Lambda1[[iters[d]]])
        lambdas <- diag(MHPCA$model$full[[d]]$Lambda1)
        names(lambdas) <- MHPCA$data[[d]] %>%
          dplyr::select(dplyr::starts_with("between")) %>%
          names()
        lambdas <- sort(lambdas, decreasing = TRUE)
        names(lambdas[1:MHPCA$FVE$G_prime[[d]]])
      }
    )

    keep_within  <- purrr::map(
      groups,
      function(d) {
        # lambdas <- diag(MHPCA$model$full[[d]]$Lambda2[[iters[d]]])
        lambdas <- diag(MHPCA$model$full[[d]]$Lambda2)
        names(lambdas) <- MHPCA$data[[d]] %>%
          dplyr::select(dplyr::starts_with("within")) %>%
          names()
        lambdas <- sort(lambdas, decreasing = TRUE)
        names(lambdas[1:MHPCA$FVE$H_prime[[d]]])
      }
    )

    MHPCA$data <- purrr::map(
      groups,
      function(d) {
        dplyr::select(
          MHPCA$data[[d]], Observation, Group, Subject, reg, func,
          y, func_ind, reg_ind, index, Overall_Mean, eta, y_centered2,
          dplyr::all_of(keep_between[[d]]), dplyr::all_of(keep_within[[d]])
        )
      }
    )
    tictoc::toc(quiet = quiet)

    tictoc::tic("    c. Final Model")
    # _c. refit model using only G' and H' components ----
    MHPCA$model$final <- purrr::map(groups, ~ mm(MHPCA, .x, maxiter, epsilon, j_tot))
    tictoc::toc(quiet = quiet)

    # _d. prediction for each record ----
    tictoc::tic("    d. Prediction")
    MHPCA$data <- purrr::map(groups, ~ predictor(MHPCA, .x, functional, regional))
    tictoc::toc(quiet = quiet)
  }

  # Return Results ---------------------------------------------------------

  tictoc::toc(quiet = quiet)

  return(MHPCA)
}
