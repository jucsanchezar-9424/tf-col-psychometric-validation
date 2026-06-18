# ============================================================
# TF-Col - Master script: runs the full analysis in order.
# Run from the repository ROOT:  source("R/00_run_all.R")
# Requires data/tfcol_dataset.xlsx (not shared; see data/README.md)
# ============================================================
source("R/01_efa_cfa_irt_convergent.R")  # EFA, CFA, IRT-GRM, convergent validity
source("R/02_bifactor.R")                 # Bifactor (CFA subsample, n = 209)
source("R/03_table1.R")                   # Table 1 (sample characteristics)

# Record software versions (Methods / AE #4)
writeLines(capture.output(sessionInfo()), file.path("outputs", "sessionInfo.txt"))
cat("\nDone. Outputs written to outputs/. Versions saved to outputs/sessionInfo.txt\n")
