#' Functional Intraclass Correlation
#'
#' @param MHPCA MHPCA object
#'
#' @return fICC (data.frame)
#' @export
#' @importFrom purrr map_dfr
functional_ICC <- function(MHPCA) {


  purrr::map_dfr(
    MHPCA$model$full,
    ~ data.frame(
      rho_dW = sum(diag(.x$Lambda1)) /
        (sum(diag(.x$Lambda1)) + sum(diag(.x$Lambda2))
        )
    ), .id = "Group"
  )
}

#' Functional Intraclass Correlation Bootstrap
#'
#' @param MHPCA MHPCA object
#' @param nboots number of bootstraps to perform
#' @param alpha alpha level
#' @param maxiter maximum iterations
#' @param epsilon epsilon for MM algorithm
#' @param fve_cutoff fraction of variance cutoff for M-HPCA
#' @param nknots number of knots for splines
#' @param quiet logical for printing timing messages
#' @param file_name name of file to save individual bootstrapped results in if desired
#'
#' @return fICC percentiles (data.frame)
#' @export
#' @importFrom dplyr tibble filter select rename group_by summarize
#' @importFrom future plan multisession sequential
#' @importFrom furrr future_map_dfr
#' @importFrom parallelly availableCores
#' @importFrom purrr map_dfr
#' @importFrom readr write_csv
#' @importFrom stats quantile
#' @importFrom tictoc tic toc
#' @importFrom tidyr pivot_wider
functional_ICC_bootstrap <- function(
  MHPCA, nboots, alpha = 0.05, maxiter, epsilon, fve_cutoff,
  nknots, quiet = FALSE, file_name = NULL
) {

  # to avoid issues with non-standard evaluation in tidyeval, set "global
  # variables" to NULL and remove them. this won't cause an issue with the rest
  # of the code.
  Repetition        <- NULL
  Group             <- NULL
  subject_resampled <- NULL
  reg               <- NULL
  func              <- NULL
  y                 <- NULL
  rho_dW            <- NULL
  probs             <- NULL
  x                 <- NULL


  rm(list = c(
    "Repetition", "Group", "subject_resampled", "reg", "func", "y", "rho_dW",
    "probs", "x"
  ))

  tictoc::tic("Functional ICC Bootstrap")
  my_quantile <- function(x, probs) {
    dplyr::tibble(x = stats::quantile(x, probs), probs = probs)
  }

  percentiles <- c(alpha / 2, 1 - alpha / 2)

  if(!is.null(file_name)) {
    readr::write_csv(
      data.frame(
        Group  = character(),
        rho_dW = double()
      ),
      file = file_name
    )
  }

  future::plan(future::multisession, workers = parallelly::availableCores() - 2, gc = TRUE)
  fICC <- furrr::future_map_dfr(1:nboots, function(b) {
    # 1. resample subjects
    resampled_data <- purrr::map_dfr(MHPCA$data, function(data) {
      subs <- unique(data$Subject)
      # sample the subjects with replacement
      sampled_subs <- sample(subs, replace = TRUE)
      # M-HPCA has issues when the subjects are repeated so add a character
      # to make them distinct
      names(sampled_subs) <- paste(sampled_subs, 1:length(sampled_subs), sep = "_")

      purrr::map_dfr(
        sampled_subs,
        ~ dplyr::filter(data, Subject == .x),
        .id = "subject_resampled"
      ) %>%
        dplyr::select(Repetition, Group, subject_resampled, reg, func, y) %>%
        dplyr::rename(Subject = subject_resampled)
    })

    # 2. estimate model components
    resampled_MHPCA <- MHPCA_decomp(
      data       = resampled_data,
      maxiter    = maxiter,
      epsilon    = epsilon,
      fve_cutoff = fve_cutoff,
      nknots     = nknots,
      reduce     = FALSE,
      quiet      = TRUE
    )

    rm(resampled_data)

    # 3. calculate functional ICC
    fICC <- functional_ICC(resampled_MHPCA)

    rm(resampled_MHPCA)

    if(!is.null(file_name)) {
      readr::write_csv(
        data.frame(fICC),
        file = file_name,
        append = TRUE
      )
    }

    fICC
  }, .id = "iteration", .progress = !quiet) %>%

    # 4. output confidence interval
    dplyr::group_by(Group) %>%
    dplyr::summarize(my_quantile(rho_dW, percentiles), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = probs, values_from = x)
  future::plan(future::sequential)
  tictoc::toc(quiet = quiet)
  return(fICC)
}
