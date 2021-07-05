
# MHPCA Bootstrap Within Test ---------------------------------------------

#' @title Multilevel-Hybrid Principal Component Analysis Bootstrap Test
#' @description Function for performing the parametric bootstrap test for
#'  comparing groups as laid out in 'Multilevel-Hybrid Principal Component
#'  Analysis' by Campos et. al (20??).
#'
#' @param MHPCA MHPCA output from MHPCA function
#' @param B number of bootstrap samples (scalar)
#' @param region region to consider for hypothesis test (scalar)
#' @param group groups to consider for the hypothesis test (character)
#' @param nknots number of knots for smoothing splines
#' @param quiet print timing messages (logical)
#'
#' @importFrom dplyr group_by slice summarize n mutate filter summarise rename select pull starts_with mutate_all all_of ungroup
#' @importFrom furrr future_map
#' @importFrom future plan multiprocess
#' @importFrom Matrix rowMeans
#' @importFrom pracma trapz
#' @importFrom purrr map_dfr pmap_dfr map_dfc map_dbl
#' @importFrom stats setNames rnorm smooth.spline predict
#' @importFrom stringr str_remove str_split str_pad
#' @importFrom tictoc tic toc
#'
#' @export
MHPCA_bootstrap_within <- function(MHPCA, B, region, group, nknots, quiet = FALSE) {

  # to avoid issues with non-standard evaluation in tidyeval, set "global
  # variables" to NULL and remove them. this won't cause an issue with the rest
  # of the code.
  functional   <- NULL
  regional     <- NULL
  Subject      <- NULL
  Group        <- NULL
  reg          <- NULL
  func         <- NULL
  Repetition  <- NULL
  overall      <- NULL
  Overall_Mean <- NULL

  rm(list = c(
    "functional", "regional", "Subject", "Group", "reg", "func", "Repetition",
    "overall", "Overall_Mean"
  ))


  tictoc::tic("Bootstrap Procedure")

  # 0. Setup Parametric Bootstrap Components ------------------------------

  tictoc::tic("0. Setup Parametric Bootstrap Components")

  # Define global variables
  functional = unique(MHPCA$data[[group]]$func)
  f_tot <- length(functional)  # total grid points in functional domain
  regional = unique(MHPCA$data[[group]]$reg)
  r_tot <- length(regional)    # total grid points in regional domain
  obs <- unique(MHPCA$data[[group]]$Repetition)
  names(obs) <- obs
  groups <- factor(names(MHPCA$data))
  names(groups) <- groups

  samp_size <- purrr::map_dfr(MHPCA$data, ~.x) %>%
    dplyr::group_by(Subject) %>%
    dplyr::slice(1) %>%
    dplyr::group_by(Group) %>%
    dplyr::summarize(count = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(Group = factor(Group))

  N <- samp_size$count
  names(N) <- samp_size$Group

  # load fixed effects from MHPCA decomposition results
  mu  <- MHPCA[["mu"]]
  eta <- MHPCA[["eta"]]

  # obtain region shift under null hypothesis by taking the point-wise average
  # over obs in region r in group d and in other regions just taking the
  # estimated eta
  eta_null <- purrr::map_dfr(
    groups,
    function(g) {
      purrr::map_dfr(
        regional,
        function(r) {
          if (r %in% region & g %in% group) {
            purrr::map_dfr(
              obs,
              # take point-wise average over obs
              ~ MHPCA[["eta"]] %>%
                dplyr::filter(Group == g, reg == r) %>%
                dplyr::group_by(Group, reg, func) %>%
                dplyr::summarise(eta_null = mean(eta), .groups = "drop") %>%
                dplyr::mutate(Repetition = !!.x)
            )
          } else {
            # just filtering out the appropriate eta
            MHPCA[["eta"]] %>%
              dplyr::filter(Group == g, reg == r) %>%
              dplyr::rename(eta_null = eta)
          }
        }) %>%
        dplyr::select(Repetition, Group, reg, func, eta_null)
    })

  phi_l    <- matrix(list(), nrow = length(groups))
  phi_m    <- matrix(list(), nrow = length(groups))
  v_k      <- matrix(list(), nrow = length(groups))
  v_p      <- matrix(list(), nrow = length(groups))
  lambda_1 <- matrix(list(), nrow = length(groups))
  lambda_2 <- matrix(list(), nrow = length(groups))
  sigma_2  <- matrix(list(), nrow = length(groups))

  for (d in groups) {
    # eigencomponents from MHPCA decomposition results
    phi_l[[d]] <- MHPCA[[4]][[1]][["functional"]][[d]][["eigendecomp"]]$vectors
    phi_m[[d]] <- MHPCA[[4]][[2]][["functional"]][[d]][["eigendecomp"]]$vectors
    v_k[[d]]   <- MHPCA[[4]][[1]][["regional"]][[d]][["eigendecomp"]]$vectors
    v_p[[d]]   <- MHPCA[[4]][[2]][["regional"]][[d]][["eigendecomp"]]$vectors

    lambda_1[[d]] <- MHPCA$model$final[[d]]$Lambda1
    lambda_2[[d]] <- MHPCA$model$final[[d]]$Lambda2
    sigma_2[[d]] <- MHPCA$model$final[[d]]$sigma2
  }

  tictoc::toc(quiet = quiet)


  # 1. Bootstrap Procedure --------------------------------------------------

  tictoc::tic("1. Bootstrap Procedure")

  # parametric bootstrap procedure
  future::plan(future::multiprocess)
  boots <- furrr::future_map(
    1:B,
    function(b) {
      # generate data
      y <- purrr::map_dfr(
        as.character(groups),
        ~ expand.grid(
          g = as.character(.x),
          c = obs,
          r = regional,
          i = 1:N[.x]
        )
      ) %>%
        purrr::pmap_dfr(function(g, c, r, i) {
          eta_null <- eta_null %>%
            dplyr::filter(Repetition == c, Group == g, reg == r) %>%
            dplyr::pull("eta_null")

          # set up level 1 and level 2 vectors
          level_1 <- vector("numeric", f_tot)
          level_2 <- vector("numeric", f_tot)

          kl <- names(dplyr::select(MHPCA$data[[g]], dplyr::starts_with("between"))) %>%
            stringr::str_remove("between_") %>%
            stringr::str_split("_", simplify = TRUE) %>%
            data.frame() %>%
            dplyr::mutate_all(as.numeric) %>%
            stats::setNames(c("k", "l"))

          pm <- names(dplyr::select(MHPCA$data[[g]], dplyr::starts_with("within"))) %>%
            stringr::str_remove("within_") %>%
            stringr::str_split("_", simplify = TRUE) %>%
            data.frame() %>%
            dplyr::mutate_all(as.numeric) %>%
            stats::setNames(c("p", "m"))

          # add level 1 components
          for (j in 1:nrow(kl)) {
            sim_score <- stats::rnorm(1, 0, sqrt(lambda_1[[g]][j, j]))
            level_1 <- level_1 + sim_score *
              v_k[[g]][match(r, regional), kl[j, "k"]] *
              phi_l[[g]][, kl[j, "l"]]
          }

          # add level 2 components
          for (j in 1:nrow(pm)) {
            sim_score <- stats::rnorm(1, 0, sqrt(lambda_2[[g]][j, j]))
            level_2 <- level_2 + sim_score *
              v_p[[g]][match(r, regional), pm[j, "p"]] *
              phi_m[[g]][, pm[j, "m"]]
          }

          data.frame(
            "Repetition" = c,
            "Group"     = g,
            "reg" = r,
            "Subject"   = paste("Subject", stringr::str_pad(i, width = 2, pad = "0")),
            "func"      = functional,
            "y"       = mu + eta_null + level_1 + level_2 +
              stats::rnorm(f_tot, 0, sqrt(sigma_2[[g]])),
            stringsAsFactors = FALSE
          )
        })

      # estimate mean function
      mu <- purrr::map_dfc(
        groups,
        function(g) {
          dat <- dplyr::filter(y, Group == g)
          spline_fit <- stats::smooth.spline(x = dat$func, y = dat$y, nknots = nknots)
          group_mean <- stats::predict(spline_fit, functional)$y
        }
      ) %>%
        dplyr::mutate(overall = Matrix::rowMeans(dplyr::select(., dplyr::all_of(groups)))) %>%
        dplyr::pull(overall)

      # center data by overall mean function
      y <- y %>%
        dplyr::group_by(Subject, Group, reg, Repetition) %>%
        dplyr::mutate(
          Overall_Mean = MHPCA[["mu"]],
          y_centered = y - Overall_Mean
        ) %>%
        dplyr::ungroup()

      # estimate region shifts

      # obtain group-condition-region-specific shifts by fitting a spline on
      # the centered data within each group-condition-region
      eta_b <- expand.grid(
        d = groups,
        c = obs,
        e = regional
      ) %>% purrr::pmap_dfr(function(d, c, e) {
        dat <- dplyr::filter(y, Group == d, reg == e, Repetition == c)
        spline_fit <- stats::smooth.spline(x = dat$func, y = dat$y_centered, nknots = nknots)
        data.frame(
          Group = d,
          reg = e,
          Repetition = c,
          func = functional,
          eta = stats::predict(spline_fit, functional)$y,
          stringsAsFactors = FALSE
        )
      })

      eta_b_null <- eta_b %>%
        dplyr::filter(Group == group, reg == region) %>%
        dplyr::group_by(Group, reg, func) %>%
        dplyr::summarize(eta_bar = mean(eta), .groups = "drop")

      return(list(
        eta_b = eta_b,
        eta_b_null = eta_b_null
      ))
    }, .progress = !quiet)


  tictoc::toc(quiet = quiet)

  # 2. Calculate p-values ---------------------------------------------------

  tictoc::tic("2. Calculate p-values")

  # calculate distribution of test statistic under null hypothesis
  ts_boots <- purrr::map_dbl(
    boots,
    function(b) {
      eta_data <- b$eta_b %>%
        dplyr::filter(Group == group, reg == region) %>%
        split(.$Repetition)
      ts <- sqrt(sum(purrr::map_dbl(
        eta_data,
        ~ pracma::trapz(functional, (.x$eta - b$eta_b_null$eta_bar)^2)
      )))
    })

  # calculate test statistic for observed data

  ts_data <- sqrt(sum(purrr::map_dbl(
    dplyr::filter(MHPCA$eta, Group == group, reg == region) %>%
      split(.$Repetition),
    ~ pracma::trapz(functional, (
      .x$eta - eta_null %>%
        dplyr::filter(
          Group == group,
          reg == region,
          Repetition == obs[1]
        ) %>%
        dplyr::pull(eta_null)
    ) ^ 2)
  )))

  # obtain p-value
  p_val <- mean(ts_boots > ts_data)

  tictoc::toc(quiet = quiet)

  # 3. Output Results -------------------------------------------------------

  # return a list with the p-value and common region shift under the null from
  # each run of the bootstrap procedure
  tictoc::toc(quiet = quiet)

  return(list(
    p_val = p_val,
    boots = boots
  ))
}


# MHPCA Bootstrap Between Test --------------------------------------------

#' @title Multilevel-Hybrid Principal Component Analysis Bootstrap Test
#' @description Function for performing the parametric bootstrap test for
#'  comparing groups as laid out in 'Multilevel-Hybrid Principal Component
#'  Analysis' by Campos et. al (20??).
#'
#' @param MHPCA MHPCA output from MHPCA function
#' @param B number of bootstrap samples (scalar)
#' @param region region to consider for hypothesis test (scalar)
#' @param group groups to consider for the hypothesis test (character)
#' @param nknots number of knots for smoothing splines
#' @param quiet print timing messages (logical)
#'
#' @importFrom dplyr group_by slice summarize n mutate filter pull summarise select starts_with mutate_all all_of ungroup
#' @importFrom furrr future_map
#' @importFrom future plan multiprocess
#' @importFrom Matrix rowMeans
#' @importFrom pracma trapz
#' @importFrom purrr map_dfr pmap_dfr map_dfc map_dbl
#' @importFrom stats setNames rnorm smooth.spline predict
#' @importFrom stringr str_remove str_split str_pad
#' @importFrom tictoc tic toc
#' @importFrom tidyr pivot_wider
#'
#' @export
MHPCA_bootstrap_between <- function(MHPCA, B, region, group, nknots, quiet = FALSE) {

  # to avoid issues with non-standard evaluation in tidyeval, set "global
  # variables" to NULL and remove them. this won't cause an issue with the rest
  # of the code.
  functional   <- NULL
  regional     <- NULL
  Subject      <- NULL
  Group        <- NULL
  reg          <- NULL
  func         <- NULL
  Repetition  <- NULL
  overall      <- NULL
  Overall_Mean <- NULL
  diff_df      <- NULL
  `.data`      <- NULL

  rm(list = c(
    "functional", "regional", "Subject", "Group", "reg", "func", "Repetition",
    "overall", "Overall_Mean", "diff_df", "`.data`"
  ))

  tictoc::tic("Bootstrap Procedure")

  # 0. Setup Parametric Bootstrap Components ------------------------------

  tictoc::tic("0. Setup Parametric Bootstrap Components")

  # Define global variables
  f_tot <- length(functional)  # total grid points in functional domain
  r_tot <- length(regional)    # total grid points in regional domain
  obs <- unique(MHPCA$data[[group[[1]]]]$Repetition)
  names(obs) <- obs
  groups <- factor(names(MHPCA$data))
  names(groups) <- groups

  samp_size <- purrr::map_dfr(MHPCA$data, ~.x) %>%
    dplyr::group_by(Subject) %>%
    dplyr::slice(1) %>%
    dplyr::group_by(Group) %>%
    dplyr::summarize(count = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(Group = factor(Group))

  N <- samp_size$count
  names(N) <- samp_size$Group

  # load fixed effects from MHPCA decomposition results
  mu  <- MHPCA[["mu"]]
  eta <- MHPCA[["eta"]]

  # calculate the average difference in etas
  eta_bar <- MHPCA$eta %>%
    dplyr::filter(Group %in% group, reg %in% region) %>%
    tidyr::pivot_wider(names_from = Repetition, values_from = eta) %>%
    dplyr::mutate(diff = .data[[obs[1]]] - .data[[obs[2]]]) %>%
    {. ->> diff_df} %>% # save the differences to calculate the test stat later
    dplyr::group_by(func) %>%
    dplyr::summarize(eta_bar = mean(diff), .groups = "drop") %>%
    dplyr::pull(eta_bar)

  # obtain region shift under null hypothesis by taking the point-wise average
  # over obs in region r in group d and in other regions just taking the
  # estimated eta
  eta_null <- purrr::map_dfr(groups, function(g) {
    purrr::map_dfr(regional, function(r) {
      if (r %in% region & g %in% group) {
        purrr::map_dfr(obs, function(c) {
          MHPCA[["eta"]] %>%
            dplyr::filter(Group == g, reg == r) %>%
            dplyr::group_by(Group, reg, func) %>%
            # calculating the average eta
            dplyr::summarise(eta_null = mean(eta), .groups = "drop")  %>%
            mutate(
              Repetition = c,
              # adding/subtracting half of the average difference
              eta_null = eta_null + (-1) ^ match(c, obs) * (1/2) * eta_bar
            )
        })
      } else {
        # just filtering out the appropriate eta
        MHPCA[["eta"]] %>%
          filter(Group == g, reg == r) %>%
          dplyr::filter(eta_null = eta)
      }
    }) %>%
      dplyr::select(Repetition, Group, reg, func, eta_null)
  })

  phi_l    <- matrix(list(), nrow = length(groups))
  phi_m    <- matrix(list(), nrow = length(groups))
  v_k      <- matrix(list(), nrow = length(groups))
  v_p      <- matrix(list(), nrow = length(groups))
  lambda_1 <- matrix(list(), nrow = length(groups))
  lambda_2 <- matrix(list(), nrow = length(groups))
  sigma_2  <- matrix(list(), nrow = length(groups))

  for (d in groups) {
    # eigencomponents from MHPCA decomposition results
    phi_l[[d]] <- MHPCA[[4]][[1]][["functional"]][[d]][["eigendecomp"]]$vectors
    phi_m[[d]] <- MHPCA[[4]][[2]][["functional"]][[d]][["eigendecomp"]]$vectors
    v_k[[d]]   <- MHPCA[[4]][[1]][["regional"]][[d]][["eigendecomp"]]$vectors
    v_p[[d]]   <- MHPCA[[4]][[2]][["regional"]][[d]][["eigendecomp"]]$vectors

    # random effects error variance
    lambda_1[[d]] <- MHPCA$model$final[[d]]$Lambda1[[length(MHPCA$model$final[[d]]$Lambda1)]]
    lambda_2[[d]] <- MHPCA$model$final[[d]]$Lambda2[[length(MHPCA$model$final[[d]]$Lambda2)]]

    # measurement error variance
    sigma_2[[d]] <- MHPCA$model$final[[d]]$sigma2[length(MHPCA$model$final[[d]]$sigma2)]
  }

  tictoc::toc(quiet = quiet)


  # 1. Bootstrap Procedure --------------------------------------------------

  tictoc::tic("1. Bootstrap Procedure")

  # parametric bootstrap procedure
  future::plan(future::multiprocess)
  boots <- furrr::future_map(
    1:B,
    function(b) {
      # generate data
      y <- purrr::map_dfr(
        as.character(groups),
        ~ expand.grid(
          g = as.character(.x),
          c = obs,
          r = regional,
          i = 1:N[.x]
        )
      ) %>%
        purrr::pmap_dfr(function(g, c, r, i) {
          eta_null <- eta_null %>%
            dplyr::filter(Repetition == c, Group == g, reg == r) %>%
            dplyr::pull("eta_null")

          # set up level 1 and level 2 vectors
          level_1 <- vector("numeric", length(functional))
          level_2 <- vector("numeric", length(functional))

          kl <- names(dplyr::select(MHPCA$data[[g]], dplyr::starts_with("between"))) %>%
            stringr::str_remove("between_") %>%
            stringr::str_split("_", simplify = TRUE) %>%
            data.frame() %>%
            dplyr::mutate_all(as.numeric) %>%
            stats::setNames(c("k", "l"))

          pm <- names(dplyr::select(MHPCA$data[[g]], dplyr::starts_with("within"))) %>%
            stringr::str_remove("within_") %>%
            stringr::str_split("_", simplify = TRUE) %>%
            data.frame() %>%
            dplyr::mutate_all(as.numeric) %>%
            stats::setNames(c("p", "m"))

          # add level 1 components
          for (j in 1:nrow(kl)) {
            sim_score <- stats::rnorm(1, 0, sqrt(lambda_1[[g]][j, j]))
            level_1 <- level_1 + sim_score *
              v_k[[g]][match(r, regional), kl[j, "k"]] *
              phi_l[[g]][, kl[j, "l"]]
          }

          # add level 2 components
          for (j in 1:nrow(pm)) {
            sim_score <- stats::rnorm(1, 0, sqrt(lambda_2[[g]][j, j]))
            level_2 <- level_2 + sim_score *
              v_p[[g]][match(r, regional), pm[j, "p"]] *
              phi_m[[g]][, pm[j, "m"]]
          }

          data.frame(
            "Repetition" = c,
            "Group"     = g,
            "reg" = r,
            "Subject"   = paste("Subject", stringr::str_pad(i, width = 2, pad = "0")),
            "func"      = functional,
            "y"       = mu + eta_null + level_1 + level_2 +
              stats::rnorm(f_tot, 0, sqrt(sigma_2[[g]])),
            stringsAsFactors = FALSE
          )
        })

      # estimate mean function
      mu <- purrr::map_dfc(
        groups,
        function(g) {
          dat <- dplyr::filter(y, Group == g)
          spline_fit <- stats::smooth.spline(x = dat$func, y = dat$y, nknots = nknots)
          group_mean <- stats::predict(spline_fit, functional)$y
        }
      ) %>%
        dplyr::mutate(overall = Matrix::rowMeans(dplyr::select(., dplyr::all_of(groups)))) %>%
        dplyr::pull(overall)

      # center data by overall mean function
      y <- y %>%
        dplyr::group_by(Subject, Group, reg, Repetition) %>%
        dplyr::mutate(
          Overall_Mean = MHPCA[["mu"]],
          y_centered = y - Overall_Mean
        ) %>%
        dplyr::ungroup()

      # estimate region shifts

      # obtain group-condition-region-specific shifts by fitting a spline on the
      # centered data within each group-condition-region
      eta_diff_b <- expand.grid(
        d = group,
        c = obs,
        e = region,
        stringsAsFactors = FALSE
      ) %>% purrr::pmap_dfr(function(d, c, e) {
        dat <- dplyr::filter(y, Group == d, reg == e, Repetition == c)
        spline_fit <- stats::smooth.spline(x = dat$func, y = dat$y_centered, nknots = nknots)
        data.frame(
          Group = d,
          reg = e,
          Repetition = c,
          func = functional,
          eta = stats::predict(spline_fit, functional)$y,
          stringsAsFactors = FALSE
        )
      }) %>%
        tidyr::pivot_wider(names_from = Repetition, values_from = eta) %>%
        dplyr::mutate(diff = .data[[obs[1]]] - .data[[obs[2]]])

      eta_bar_b <- eta_diff_b %>%
        dplyr::group_by(func) %>%
        dplyr::summarize(eta_bar = mean(diff), .groups = "drop")

      return(list(
        eta_diff_b = eta_diff_b,
        eta_bar_b  = eta_bar_b
      ))
    })


  tictoc::toc(quiet = quiet)

  # 2. Calculate p-values ---------------------------------------------------

  tictoc::tic("2. Calculate p-values")

  # calculate distribution of test statistic under null hypothesis
  ts_boots <- purrr::map_dbl(
    boots,
    function(b) {
      eta_data <- b$eta_diff_b %>%
        split(.$Group, drop = TRUE)
      ts <- sqrt(sum(purrr::map_dbl(
        eta_data,
        ~ pracma::trapz(functional, (.x$diff - b$eta_bar_b$eta_bar)^2)
      )))
    })

  # calculate test statistic for observed data
  ts_data <- sqrt(sum(
    diff_df %>%
      split(.$Group, drop = TRUE) %>%
      purrr::map_dfr(~ pracma::trapz(functional, (.x$diff - eta_bar)^2))
  ))

  # obtain p-value
  p_val <- mean(ts_boots > ts_data)

  tictoc::toc(quiet = quiet)

  # 3. Output Results -------------------------------------------------------

  # return a list with the p-value and common region shift under the null from
  # each run of the bootstrap procedure
  tictoc::toc(quiet = quiet)

  return(list(
    p_val = p_val,
    boots = boots
  ))
}

