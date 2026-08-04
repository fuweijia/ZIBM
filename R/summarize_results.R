#' Summarize mediation effects across taxa
#'
#' @param ZIBM_res Output list from \code{ZIBM()}.
#' @param p_adjust_method Method passed to \code{p.adjust}. Default "BH".
#' @return A Data frame with one row per taxon: taxa_name, mediation_effect_NIE1/NIE2/NDE/NIE,
#' NIE_sd_NIE1/NIE2/NDE/NIE, p_value_NIE1/NIE2/NDE/NIE, adjusted_p_value_NIE1/NIE2/NDE/NIE.
#' Taxa where all k failed to fit are dropped.
#' @importFrom stats setNames
#' @export
#'
summarize_results <- function(ZIBM_res, p_adjust_method = "BH") {
  eff_names <- c("NIE1", "NIE2", "NDE", "NIE")
  rows <- lapply(ZIBM_res$list_save, function(taxon_res) {
    res <- taxon_res$res_fin_med
    if (!is.list(res)) {
       return(NULL)
     }
    sd_vec <- res$NIE_sd
    if (inherits(sd_vec, "try-error") || length(sd_vec) != 4) {
      sd_vec <- setNames(rep(NA_real_, 4), eff_names)
    }

    row <- c(
      list(taxa_name = res$taxon_name),
      setNames(as.list(res$mediation_effect[eff_names]), paste0("mediation_effect_", eff_names)),
      setNames(as.list(sd_vec[eff_names]), paste0("NIE_sd_", eff_names)),
      setNames(as.list(res$p_value[eff_names]), paste0("p_value_", eff_names))
    )
    as.data.frame(row, stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, rows)

  for (eff in eff_names) {
    out[[paste0("adjusted_p_value_", eff)]] <-
      p.adjust(out[[paste0("p_value_", eff)]], method = "BH")
  }
  return(out)
}
