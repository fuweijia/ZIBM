# ZIBM
An R package for mediation analysis of microbiome relative abundance data using a zero-inflated beta mixture model. The package provides functions for model fitting, parameter estimation, inference, and model selection.

## Installation
```
install.packages("remotes")
remotes::install_github("fuweijia/ZIBM")
```
## Usage

The main function estimates direct and indirect effects and standard errors, and calculates unadjusted and adjusted p-values using Benjamini-Hochberg adjustment. Use microbiome data, CovData (data frame with outcome, exposure, library size, and confounders) as inputs.
```
n_sub  <- 40
n_taxa <- 10
lib_size   <- sample(8000:12000, n_sub, replace = TRUE)
taxon_prop <- matrix(rbeta(n_sub * n_taxa, 0.5, 20), nrow = n_sub)
MicrobData <- matrix(
   rbinom(n_sub * n_taxa, size = lib_size, prob = taxon_prop),
   nrow = n_sub, ncol = n_taxa,
   dimnames = list(NULL, paste0("taxon", seq_len(n_taxa)))
 )
CovData <- data.frame(
   Y = rnorm(n_sub),
   X = rbinom(n_sub, 1, 0.5),
   Z = rnorm(n_sub),
  lib = lib_size
)
res <- ZIBM(
  MicrobData = MicrobData,
  CovData    = CovData,
  lib_name   = "lib",
  y_name     = "Y",
  x_name     = "X",
  conf_name  = "Z",
  k_range    = c(1:3),
  num_cores  = 4
)

# Summary table of mediation effects, SEs, and p-values across taxa
summary_tab <- summarize_mediation(res)
head(summary_tab)
```
## Output structure
 
`ZIBM()` returns:
 
```
list(
  list_save        # one entry per taxon
  nTaxa            # number of taxa
  nSub             # number of subjects
  taxon_ori_name   # original taxon names
)
```

The `list_save[[i]]` contains:
 
- `res_fin_med` — the AIC-selected model fit for that taxon 
- `res_list_med` — model fits for every candidate `k` in `k_range`
Each fit object includes a decomposed mediation effect with four components:
 
| Component | Interpretation |
|---|---|
| `NIE1` | Indirect effect through mediator presence (zero vs. non-zero) |
| `NIE2` | Indirect effect through mediator abundance, given presence |
| `NDE`  | Natural direct effect |
| `NIE`  | Total natural indirect effect (`NIE1 + NIE2`) |
 
along with matching standard errors (`NIE_sd`), p-values (`p_value`) and adjusted p-values (`adjusted_p_value`).
 

