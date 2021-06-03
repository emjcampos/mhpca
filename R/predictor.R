#' Subject trajectory estimation
#'
#' @param MHPCA MHPCA object
#' @param d group
#' @param functional functional dimension
#' @param regional regional dimension
#'
#' @importFrom dplyr select starts_with slice mutate pull filter
#' @importFrom Matrix as.matrix kronecker diag
#' @importFrom purrr map_dfr
predictor <- function(MHPCA, d, functional, regional) {

  # to avoid issues with non-standard evaluation in tidyeval, set "global
  # variables" to NULL and remove them. this won't cause an issue with the rest
  # of the code.
  Subject      <- NULL
  Repetition  <- NULL

  rm(list = c(
    "Subject", "Repetition"
  ))


  subjects        <- unique(MHPCA$data[[d]]$Subject)      # subject IDs
  names(subjects) <- subjects                             # name subjects
  n               <- length(subjects)                     # number of subjects
  f_tot           <- length(functional)                   # grid pts functional
  r_tot           <- length(regional)                     # grid pts regional
  obs             <- unique(MHPCA$data[[d]]$Repetition)  # conditions
  names(obs)      <- obs                                  # name conditions
  c_tot           <- length(obs)                          # number of conditions

  Phi1 = MHPCA$data[[d]] %>%
    dplyr::select(dplyr::starts_with("between")) %>%
    dplyr::slice(1:(f_tot * r_tot)) %>%
    Matrix::as.matrix()

  Phi2 = MHPCA$data[[d]] %>%
    dplyr::select(dplyr::starts_with("within")) %>%
    dplyr::slice(1:(f_tot * r_tot)) %>%
    Matrix::as.matrix()

  z = Matrix::kronecker(rep(1, c_tot), Phi1) %>%
    data.frame() %>%
    dplyr::mutate(Repetition = rep(obs, each = r_tot * f_tot))

  w = Matrix::kronecker(Matrix::diag(c_tot), Phi2) %>%
    data.frame() %>%
    dplyr::mutate(Repetition = rep(obs, each = r_tot * f_tot))

  z_i = vector("list", n)
  names(z_i) = subjects
  w_i = vector("list", n)
  names(w_i) = subjects

  c_i = map(
    subjects,
    ~ MHPCA$data[[d]][MHPCA$data[[d]]$Subject == subjects[.x], ] %>%
      dplyr::pull(Repetition) %>% unique()
  )

  for (i in subjects) {
    z_i[[i]] = dplyr::filter(z, Repetition %in% c_i[[i]]) %>% dplyr::select(-Repetition) %>% Matrix::as.matrix()
    w_i[[i]] = dplyr::filter(w, Repetition %in% c_i[[i]]) %>% dplyr::select(-Repetition) %>% Matrix::as.matrix()
  }

  purrr::map_dfr(
    subjects,
    function(i) {
      df <- dplyr::filter(MHPCA$data[[d]], Subject == i)
      df$predicted <- df$Overall_Mean + df$eta +
        as.vector(z_i[[i]] %*% MHPCA$model$final[[d]]$zeta_i[[i]]) +
        as.vector(w_i[[i]] %*% MHPCA$model$final[[d]]$xi_i[[i]])
      return(df)
    })
}
