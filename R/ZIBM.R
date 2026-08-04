#' The main function estimates direct and indirect effects and standard errors, and calculates unadjusted and adjusted p-values (Benjamini-Hochberg adjustment).
#'
#' @param MicrobData Matrix/data frame of microbiome data (subjects x taxa).
#' @param CovData Data frame with outcome, exposure, library size, and
#'   confounder columns, in the same subject order as MicrobData.
#' @param lib_name Column Name of library size variable within CovData.
#' @param y_name Column name in CovData for outcome Y.
#' @param x_name Column name in CovData for covariate of interest X.
#' @param conf_name Name of confounders within CovData. Default is NULL, meaning no confounder.
#' @param k_range Integer vector of candidate mixture-component counts to try per taxon; best k chosen by AIC.
#' @param num_cores Number of cores for parallel fitting.
#' @param zero_prop_NIE2 Proportion-of-zeros threshold for using the zero-inflated variant. Default 0.1.
#' @param zero_count_NIE2 Absolute zero-count threshold. Default 4 * (length(conf_name) + 2)/
#' @param x4_inter Whether to include the interaction term \mjeqn{\beta_4}{} Default is TRUE.
#' @param x5_inter Whether to include the interaction term \mjeqn{\beta_5}{} Default is TRUE.
#' @return List with:
#' \describe{
#'   \item{list_save}{A list per-taxon of res_fin_med (AIC-selected fit, or NA if all k failed) and
#'   res_list_med (fits for every k).}
#'   \item{nTaxa}{Number of taxa.}
#'   \item{nSub}{Number of subjects.}
#'   \item{taxon_ori_name}{Original taxon names.}
#' }
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom NlcOptim solnl
#' @importFrom numDeriv jacobian hessian grad
#' @examples
#' \donttest{
#' n_sub  <- 40
#' n_taxa <- 10
#' lib_size   <- sample(8000:12000, n_sub, replace = TRUE)
#' taxon_prop <- matrix(rbeta(n_sub * n_taxa, 0.5, 20), nrow = n_sub)
#' MicrobData <- matrix(
#'   rbinom(n_sub * n_taxa, size = lib_size, prob = taxon_prop),
#'   nrow = n_sub, ncol = n_taxa,
#'   dimnames = list(NULL, paste0("taxon", seq_len(n_taxa)))
#' )
#'
#' CovData <- data.frame(
#'   Y   = rnorm(n_sub),
#'   X   = rbinom(n_sub, 1, 0.5),
#'   Z   = rnorm(n_sub),
#'   lib = lib_size
#' )
#'
#' res <- ZIBM(
#'   MicrobData = MicrobData,
#'   CovData    = CovData,
#'   lib_name   = "lib",
#'   y_name     = "Y",
#'   x_name     = "X",
#'   conf_name  = "Z",
#'   k_range    = 1,
#'   num_cores  = 1
#' )
#' res_table <- summarize_results(ZIBM_res = res)
#' }
#' @export


ZIBM <- function(MicrobData,
           CovData,
           lib_name,
           y_name,
           x_name,
           conf_name,
           k_range,
           num_cores,
           zero_prop_NIE2 = 0.1,
           zero_count_NIE2 = 4 * (length(conf_name) + 2),
           x4_inter = TRUE,
           x5_inter = TRUE) {
    num_taxon <- ncol(MicrobData)
    num_sub <- nrow(MicrobData)
    taxon_ori_name <- colnames(MicrobData)

    yi_vec <- CovData[, y_name]
    xi_vec <- CovData[, x_name]
    li_vec <- CovData[, lib_name]
    conf_mat <- as.matrix(CovData[, conf_name, drop = FALSE])

    trial <- seq_len(num_taxon)

    registerDoParallel(num_cores)

    i <- numeric(0)

    list_save <- foreach(i = trial) %dopar% {
      obs_m_vec <- MicrobData[, i]
      temp_name <- colnames(MicrobData)[i]
      res_list_med <- list()
      AIC_select <- c()
      for (k in k_range) {
        if (k == 1) {
          if (sum(obs_m_vec == 0) > min(zero_prop_NIE2 * length(obs_m_vec), zero_count_NIE2)) {
            res_temp <-
              try(real_data_run_func_nomix(yi_vec, obs_m_vec, xi_vec, li_vec, conf_mat, x4_inter, x5_inter, k),
                  TRUE)
          } else {
            res_temp <-
              try(real_data_run_func_nz_nomix(yi_vec, obs_m_vec, xi_vec, li_vec, conf_mat, x4_inter, x5_inter, k),
                  TRUE)
          }
        } else {
          if (sum(obs_m_vec == 0) > min((zero_prop_NIE2 * length(obs_m_vec)), zero_count_NIE2)) {
            res_temp <-
              try(real_data_run_func(yi_vec, obs_m_vec, xi_vec, li_vec, conf_mat, x4_inter, x5_inter, k),
                  TRUE)
          } else {
            res_temp <-
              try(real_data_run_func_nz(yi_vec, obs_m_vec, xi_vec, li_vec, conf_mat, x4_inter, x5_inter, k),
                  TRUE)
          }
        }
        if (inherits(res_temp, "try-error")) {
          AIC_select[k] <- NA
        } else if (any(is.na(res_temp$mediation_effect))) {
          AIC_select[k] <- NA
        } else {
          AIC_select[k] <- res_temp$AIC_est
        }
        res_temp$taxon_name <- temp_name
        res_list_med[[k]] <- res_temp
      }
      if (all(is.na(AIC_select))) {
        res_fin_med <- NA
      } else {
        res_fin_med <- res_list_med[[which.min(AIC_select)]]
      }

      return(list(res_fin_med = res_fin_med, res_list_med = res_list_med))
    }

    return(
      list(
        list_save = list_save,
        nTaxa = num_taxon,
        nSub = num_sub,
        taxon_ori_name = taxon_ori_name
      )
    )
  }
