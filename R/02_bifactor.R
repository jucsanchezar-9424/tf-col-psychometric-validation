# ============================================================
# TF-Col - Bifactor sensitivity analysis (Schmid-Leiman, psych::omega)
# Reproduces: Supplementary Tables S9 (indices) and S10 (SL loadings)
# IMPORTANT: the bifactor model is estimated on the SAME CFA subsample
# (n = 209) used in the main analysis, reconstructed deterministically
# from outputs/split_indices_efa.csv. Running omega() on the full sample
# yields different indices and does NOT reproduce the published values.
# Run 01_efa_cfa_irt_convergent.R first (it creates split_indices_efa.csv),
# or ensure outputs/split_indices_efa.csv is present.
# ============================================================

pkgs <- c("readxl", "dplyr", "janitor", "psych", "readr")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

DATA_FILE <- "data/tfcol_dataset.xlsx"   # Dataset NO incluido (ver data/README.md)
out_dir   <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir)

# 1. Datos y 20 ítems ----------------------------------------
datos <- janitor::clean_names(readxl::read_excel(DATA_FILE))

items_IF  <- paste0("imf",  sprintf("%02d", 1:10), "_if")
items_IPS <- paste0("imps", sprintf("%02d", 1:10), "_ips")
items_all <- c(items_IF, items_IPS)

faltantes <- items_all[!items_all %in% names(datos)]
if (length(faltantes) > 0) stop("Ítems no encontrados: ", paste(faltantes, collapse = ", "))

df <- datos %>%
  dplyr::select(dplyr::all_of(items_all)) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.numeric(as.character(.x))))

# 2. Reconstruir la submuestra CFA (n = 209) -----------------
# El split principal usó set.seed(123) y guardó los índices EFA.
# La submuestra CFA = filas que NO están en split_indices_efa.csv.
split_path <- file.path(out_dir, "split_indices_efa.csv")
if (!file.exists(split_path))
  stop("Falta ", split_path, ". Corre primero 01_efa_cfa_irt_convergent.R.")
idx_efa <- readr::read_csv(split_path, show_col_types = FALSE)$id
df_cfa  <- df[-idx_efa, ]
cat("Submuestra CFA: n =", nrow(df_cfa), "(esperado 209)\n")

# 3. Bifactor (Schmid-Leiman) sobre la submuestra CFA --------
R_poly <- psych::polychoric(df_cfa)$rho
omega_res <- suppressWarnings(
  psych::omega(m = R_poly, nfactors = 2, rotate = "oblimin",
               digits = 3, title = "TF-Col Bifactor", plot = FALSE)
)

# 4. Índices (S9) --------------------------------------------
omega_t  <- omega_res$omega.tot           # ωt total
omega_h  <- omega_res$omega_h             # ωh jerárquico (guion bajo)
ws_IF    <- omega_res$omega.group[2, "group"]   # ωs específico IF
ws_IPS   <- omega_res$omega.group[3, "group"]   # ωs específico IPS
ecv      <- as.numeric(omega_res$ECV)     # ECV general
# PUC: proporción de correlaciones no contaminadas (estructural)
k <- length(items_all); g1 <- length(items_IF); g2 <- length(items_IPS)
within  <- choose(g1, 2) + choose(g2, 2)
total   <- choose(k, 2)
puc     <- (total - within) / total

s9 <- data.frame(
  Index = c("omega total (wt)", "omega hierarchical (wh)",
            "ws Financial Impact (ws-IF)", "ws Psychosocial Impact (ws-IPS)",
            "ECV (general factor)", "PUC"),
  Value = round(c(omega_t, omega_h, ws_IF, ws_IPS, ecv, puc), 3),
  N_cfa_subsample = nrow(df_cfa)
)
readr::write_csv(s9, file.path(out_dir, "supp_table_s9_bifactor_indices.csv"))

# 5. Schmid-Leiman loadings (S10) ----------------------------
sl <- as.data.frame(unclass(omega_res$schmid$sl))
sl$item   <- rownames(sl)
sl$Domain <- ifelse(sl$item %in% items_IF, "Financial Impact", "Psychosocial Impact")
g_col  <- grep("^g$",  names(sl), value = TRUE)
f_cols <- grep("^F",   names(sl), value = TRUE)
h2_col <- grep("^h2$", names(sl), value = TRUE)
s10 <- data.frame(
  item   = sl$item,
  Domain = sl$Domain,
  g      = round(sl[[g_col]], 3),
  `F1*`  = round(sl[[f_cols[1]]], 3),
  `F2*`  = round(sl[[f_cols[2]]], 3),
  h2     = round(sl[[h2_col]], 3),
  check.names = FALSE
)
readr::write_csv(s10, file.path(out_dir, "supp_table_s10_schmid_leiman_loadings.csv"))

# 6. Resultado + verificación contra los valores del manuscrito
cat("\n--- Bifactor indices (CFA subsample) ---\n"); print(s9, row.names = FALSE)
esperado <- c(0.972, 0.700, 0.284, 0.297, 0.694, 0.526)
obtenido <- s9$Value
cat("\nValores del manuscrito: wt=0.972 wh=0.700 ws-IF=0.284 ws-IPS=0.297 ECV=0.694 PUC=0.526\n")
if (all(abs(obtenido - esperado) <= 0.005)) {
  cat("OK: reproduce los valores publicados.\n")
} else {
  cat("ATENCIÓN: difiere del manuscrito. Diferencias:\n")
  print(data.frame(index = s9$Index, obtenido, esperado, dif = round(obtenido - esperado, 3)),
        row.names = FALSE)
}
