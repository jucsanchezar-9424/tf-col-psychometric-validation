# ============================================================
# TF-Col - Main analysis (EFA / CFA / IRT-GRM / convergent validity)
# Reproduces: Tables 3, 4, 6 and Figures 1-3, plus Supp. S2-S4.
# Run from the repository root. Requires data/tfcol_dataset.xlsx
# (not shared; see data/README.md).
# ============================================================


#########################################################
# TF-Col: Cross-validated EFA/CFA (clean, reproducible)
# Dataset: datos_raw.xlsx (348x85)
#########################################################

rm(list = ls())

pkgs <- c("readxl","dplyr","stringr","psych","lavaan","semTools","GPArotation","semPlot","readr","janitor")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(readxl)
library(dplyr)
library(stringr)
library(psych)
library(lavaan)
library(semTools)
library(GPArotation)
library(semPlot)
library(readr)
library(janitor)

set.seed(123)

# ---------------------------
# 1) Load data
# ---------------------------

library(readxl)
DATA_FILE <- "data/tfcol_dataset.xlsx"   # Dataset NO incluido por confidencialidad (ver data/README.md)
datos_raw <- read_excel(DATA_FILE)
# View(datos_raw)

datos <- datos_raw




# (Opcional) limpiar nombres para trabajar más fácil
# Esto crea nombres tipo: tipo_de_poblacion, 5_sexo_paciente, etc.
datos_clean <- janitor::clean_names(datos)

# ---------------------------
# 2) Identify TF-Col item columns
# ---------------------------
# Opción A (recomendada si tus ítems tienen un prefijo consistente):
# Ej: tf01, tf02..., o imf01_if..., imps01_ips...
# Ajusta el patrón según cómo estén nombrados en tu archivo.
item_candidates <- names(datos_clean)[
  str_detect(names(datos_clean), "^imf\\d+_if$|^imps\\d+_ips$|^tf\\d+$|^item\\d+$")
]

cat("Items detectados automáticamente:\n")
print(item_candidates)

# >>>> Si detectó exactamente 20 items, usamos esa lista:
if (length(item_candidates) == 20) {
  items <- item_candidates
} else {
  # Opción B: define manualmente los 20 ítems (más seguro)
  # IMPORTANTE: reemplaza estos nombres por los reales en datos_clean
  items <- c(
    "imf01_if","imf02_if","imf03_if","imf04_if","imf05_if",
    "imf06_if","imf07_if","imf08_if","imf09_if","imf10_if",
    "imps01_ips","imps02_ips","imps03_ips","imps04_ips","imps05_ips",
    "imps06_ips","imps07_ips","imps08_ips","imps09_ips","imps10_ips"
  )
}

# View(datos_clean)

stopifnot(all(items %in% names(datos_clean)))

# ---------------------------
# 3) Prepare ordinal item data
# ---------------------------
df <- datos_clean %>%
  select(all_of(items))

# Convertir a ordered factors (ordinal).
# Si vienen como "1","2","3" en texto, esto lo maneja.
df <- df %>%
  mutate(across(everything(), ~{
    x <- .x
    # si está como carácter, intenta convertir a numeric conservando NA
    if (is.character(x)) {
      x <- str_trim(x)
      x <- na_if(x, "")
      suppressWarnings(xn <- as.numeric(x))
      if (sum(!is.na(xn)) > 0) x <- xn
    }
    # factor ordinal
    return(as.ordered(x))
  }))

# Quitar filas con muchos NA (si aplica)
# (No lo hagas si ya tienes buena completitud; esto es opcional)
# df <- df %>% filter(rowSums(is.na(.)) <= 2)

n_total <- nrow(df)

# ---------------------------
# 4) Split sample (40% EFA, 60% CFA)
# ---------------------------
idx_efa <- sample(seq_len(n_total), size = round(0.40 * n_total))
df_efa <- df[idx_efa, ]
df_cfa <- df[-idx_efa, ]

out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir)
write_csv(tibble(id = idx_efa), file.path(out_dir, "split_indices_efa.csv"))

# ---------------------------
# 5) EFA (polychoric + oblimin)
# ---------------------------
R_poly <- psych::polychoric(df_efa)$rho
saveRDS(R_poly, file.path(out_dir, "efa_polychoric_matrix.rds"))

kmo <- psych::KMO(R_poly)
bart <- psych::cortest.bartlett(R_poly, n = nrow(df_efa))
kmo

writeLines(c(
  paste0("EFA n = ", nrow(df_efa)),
  paste0("KMO overall MSA: ", round(kmo$MSA, 3)),
  paste0("Bartlett chi2: ", round(bart$chisq, 2), " | df: ", bart$df, " | p: ", format.pval(bart$p.value, digits = 3))
), con = file.path(out_dir, "efa_factorability.txt"))

png(file.path(out_dir, "efa_parallel_analysis.png"), width = 1400, height = 900, res = 150)
psych::fa.parallel(R_poly, n.obs = nrow(df_efa), fm = "ml", fa = "fa",
                   main = "Parallel analysis (EFA sample)")
dev.off()

efa <- psych::fa(R_poly, nfactors = 1, n.obs = nrow(df_efa), fm = "ml", rotate = "oblimin")
saveRDS(efa, file.path(out_dir, "efa_model_object.rds"))

# Export loadings (Table 3 base)
efa_load <- as.data.frame(unclass(efa$loadings))
efa_load$item <- rownames(efa_load)
efa_load <- efa_load %>% relocate(item)
write_csv(efa_load, file.path(out_dir, "table3_efa_loadings_raw.csv"))

# Communalities
efa_comm <- tibble(item = names(efa$communality), communality = as.numeric(efa$communality))
write_csv(efa_comm, file.path(out_dir, "efa_communalities.csv"))

# EFA fit indices (TLI / RMSEA from psych::fa)
efa_fit <- tibble(
  n_efa = nrow(df_efa),
  TLI = as.numeric(efa$TLI),
  RMSEA = as.numeric(efa$RMSEA[1]),
  RMSEA_low90 = as.numeric(efa$RMSEA[2]),
  RMSEA_high90 = as.numeric(efa$RMSEA[3])
)
write_csv(efa_fit, file.path(out_dir, "efa_fit_indices.csv"))

png(file.path(out_dir, "efa_diagram.png"), width = 1600, height = 1000, res = 150)
psych::fa.diagram(efa, simple = TRUE)
dev.off()


efa$Vaccounted
write.csv(efa$Vaccounted, file.path(out_dir, "efa_variance_accounted.csv"))

print(efa$loadings, cutoff = 0.30)
capture.output(print(efa$loadings, cutoff = 0.30, digits = 3),
               file = file.path(out_dir, "efa_loadings_print.txt"))



# ---------------------------
# 6) CFA (WLSMV, ordinal)
# ---------------------------
# Define your 2-factor model:
# IMPORTANTE: aquí necesitas asignación ítem->factor.
# Si tus nombres ya distinguen IF vs IPS (como imf vs imps), esto funciona.
# Si no, te doy abajo cómo generarlo automáticamente.

if_all <- items[str_detect(items, "^imf")]
ips_all <- items[str_detect(items, "^imps")]

stopifnot(length(if_all) > 0, length(ips_all) > 0)

model_bi <- paste0(
  "IF  =~ ", paste(if_all, collapse = " + "), "\n",
  "Ips =~ ", paste(ips_all, collapse = " + ")
)

model_uni <- paste0(
  "TF =~ ", paste(items, collapse = " + ")
)

cfa_bi <- lavaan::cfa(model_bi, data = df_cfa, ordered = colnames(df_cfa),
                      estimator = "WLSMV", std.lv = TRUE)
cfa_uni <- lavaan::cfa(model_uni, data = df_cfa, ordered = colnames(df_cfa),
                       estimator = "WLSMV", std.lv = TRUE)

saveRDS(cfa_bi, file.path(out_dir, "cfa_bi_object.rds"))
saveRDS(cfa_uni, file.path(out_dir, "cfa_uni_object.rds"))

fit_bi <- lavaan::fitMeasures(cfa_bi, c("cfi","tli","rmsea","rmsea.ci.lower","rmsea.ci.upper","srmr"))
fit_uni <- lavaan::fitMeasures(cfa_uni, c("cfi","tli","rmsea","rmsea.ci.lower","rmsea.ci.upper","srmr"))

fit_tbl <- bind_rows(
  tibble(model = "Unidimensional", n_cfa = nrow(df_cfa), !!!as.list(fit_uni)),
  tibble(model = "Correlated 2-factor", n_cfa = nrow(df_cfa), !!!as.list(fit_bi))
)
write_csv(fit_tbl, file.path(out_dir, "table4_cfa_fit_indices.csv"))

# WLSMV nested comparison
diff_test <- lavaan::lavTestLRT(cfa_uni, cfa_bi)
capture.output(diff_test, file = file.path(out_dir, "cfa_model_comparison_lrt.txt"))

# Loadings + SE (for Table 3 CFA columns or supplement)
pe_bi <- lavaan::parameterEstimates(cfa_bi, standardized = TRUE)
load_bi <- pe_bi %>%
  dplyr::filter(op == "=~") %>%
  dplyr::transmute(
    factor = lhs,
    item   = rhs,
    loading_std = std.all,
    SE = se,
    CI_low = ci.lower,
    CI_high = ci.upper
  )

readr::write_csv(load_bi, file.path(out_dir, "cfa_standardized_loadings_se_ci.csv"))



# Factor correlation (and CI if available)
psi_std <- lavaan::lavInspect(cfa_bi, "std")$psi
write_csv(as.data.frame(psi_std), file.path(out_dir, "cfa_factor_correlations_std.csv"))

ci_bi <- tryCatch(confint(cfa_bi), error = function(e) NULL)
if (!is.null(ci_bi)) {
  ci_rows <- ci_bi[grep("IF~~Ips|Ips~~IF", rownames(ci_bi)), , drop = FALSE]
  write.csv(ci_rows, file.path(out_dir, "cfa_factor_correlation_ci.csv"))
}

png(file.path(out_dir, "cfa_bi_diagram.png"), width = 1800, height = 1100, res = 150)
semPlot::semPaths(cfa_bi, "std", whatLabels = "std",
                  layout = "tree", edge.label.cex = 0.8,
                  residuals = FALSE, intercepts = FALSE, nCharNodes = 0)
dev.off()




png(file.path(out_dir, "cfa_bi_diagram_loadings.png"),
    width = 2600, height = 1500, res = 220)

semPaths(
  cfa_bi,
  what = "path",
  whatLabels = "std",          # cargas estandarizadas
  style = "lisrel",
  layout = "tree2",
  rotation = 2,
  residuals = FALSE,
  intercepts = FALSE,
  thresholds = FALSE,
  nCharNodes = 0,
  sizeLat = 13,
  sizeMan = 6,
  edge.label.cex = 0.75,       # <- baja el tamaño de números
  label.cex = 0.85,
  edge.color = "gray35",
  fade = FALSE,
  curvePivot = TRUE,
  mar = c(6, 6, 6, 6),
  digits = 2                   # <- 2 decimales (clave)
)

dev.off()




pdf(file.path(out_dir, "cfa_bi_diagram_clean.pdf"),
    width = 11, height = 6.5)

semPaths(
  cfa_bi,
  what = "path",
  whatLabels = "none",
  style = "lisrel",
  layout = "tree2",
  rotation = 2,
  residuals = FALSE,
  intercepts = FALSE,
  thresholds = FALSE,
  nCharNodes = 0,
  sizeLat = 12,
  sizeMan = 6,
  edge.color = "gray30",
  color = list(lat = "white", man = "white"),
  label.cex = 0.9,
  mar = c(6, 6, 6, 6)
)

dev.off()

# ---------------------------
# 7) Reproducibility
# ---------------------------
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("Done. Outputs saved to: ", normalizePath(out_dir))


labs <- c(paste0("IF", 1:10), paste0("IPS", 1:10))
names(labs) <- c(paste0("imf", sprintf("%02d",1:10), "_if"),
                 paste0("imps", sprintf("%02d",1:10), "_ips"))


# View(df_cfa)





pkgs <- c("readxl","dplyr","stringr","psych","lavaan","semTools","GPArotation",
          "semPlot","readr","janitor","mirt","ggplot2","tidyr","purrr")




######################################################
######cosntruir un path mas decente para publicar.
######################################################

library(semPlot)

# Etiquetas cortas para ítems
labs_items <- c(paste0("IF", 1:10), paste0("IPS", 1:10))
names(labs_items) <- c(paste0("imf", sprintf("%02d", 1:10), "_if"),
                       paste0("imps", sprintf("%02d", 1:10), "_ips"))

# Construir labels en el orden interno que usa semPlot
spm <- semPlot::semPlotModel(cfa_bi)
node_names <- spm@Vars$name

node_labels <- node_names
node_labels[node_labels %in% names(labs_items)] <- labs_items[node_labels[node_labels %in% names(labs_items)]]
node_labels[node_labels == "IF"]  <- "Financial\nImpact"
node_labels[node_labels == "Ips"] <- "Psychosocial\nImpact"

# ---------- PDF (recomendado para publicación) ----------
pdf(file.path(out_dir, "Figure_1_CFA_path_diagram_PAPER.pdf"),
    width = 12, height = 6.5)

semPlot::semPaths(
  cfa_bi,
  what = "path",
  whatLabels = "std",
  nodeLabels = node_labels,
  style = "lisrel",
  layout = "tree2",
  rotation = 2,
  reorder = FALSE,
  
  residuals = FALSE,
  intercepts = FALSE,
  thresholds = FALSE,
  
  nCharNodes = 0,
  fade = FALSE,
  
  sizeLat = 12,
  sizeMan = 5,
  edge.label.cex = 0.6,
  label.cex = 1.05,
  edge.color = "gray25",
  color = list(lat = "white", man = "white"),
  mar = c(6, 7, 5, 2),
  digits = 2
)

dev.off()

# ---------- PNG (vista rápida) ----------
png(file.path(out_dir, "Figure_1_CFA_path_diagram_PAPER.png"),
    width = 3600, height = 1950, res = 300)

semPlot::semPaths(
  cfa_bi,
  what = "path",
  whatLabels = "std",
  nodeLabels = node_labels,
  style = "lisrel",
  layout = "tree2",
  rotation = 2,
  reorder = FALSE,
  residuals = FALSE,
  intercepts = FALSE,
  thresholds = FALSE,
  nCharNodes = 0,
  fade = FALSE,
  sizeLat = 14,
  sizeMan = 7,
  edge.label.cex = 0.85,
  label.cex = 1.05,
  edge.color = "gray25",
  color = list(lat = "white", man = "white"),
  mar = c(6, 7, 5, 2),
  digits = 2
)

dev.off()





#### vamos mejorando

pdf("Figure_1_CFA_path_diagram_PAPER_loadings.pdf",
    width = 11, height = 6.5)

semPaths(
  cfa_bi,
  what = "path",
  whatLabels = "std",           # cargas estandarizadas
  style = "lisrel",
  layout = "tree2",
  rotation = 2,
  residuals = FALSE,
  intercepts = FALSE,
  thresholds = FALSE,
  nCharNodes = 0,
  shapeLat = "circle",
  shapeMan = "ellipse", 
  
  # tamaños ajustados
  sizeLat = 10,
  sizeMan = 6.5,
  edge.label.cex = 0.95,         # 
  label.cex = 0.75,
  digits = 2,
  
  edge.color = "gray35",
  color = list(lat = "white", man = "white"),
  
  mar = c(6, 6, 6, 6)
)

dev.off()

#another




library(semPlot)
library(qgraph)

# ------------------------------------------------------------
# 1) Mapa de labels cortos para ítems + nombres bonitos factores
# ------------------------------------------------------------
item_map <- c(paste0("IF", 1:10), paste0("IPS", 1:10))
names(item_map) <- c(paste0("imf", sprintf("%02d", 1:10), "_if"),
                     paste0("imps", sprintf("%02d", 1:10), "_ips"))

# Modelo semPlot (para saber el orden real de nodos)
spm <- semPlot::semPlotModel(cfa_bi)

nodeLabels <- spm@Vars$name
nodeLabels <- ifelse(nodeLabels == "IF",  "Financial\nImpact",
                     ifelse(nodeLabels == "Ips", "Psychosocial\nImpact",
                            dplyr::recode(nodeLabels, !!!item_map, .default = nodeLabels)))

# ------------------------------------------------------------
# 2) Crear el gráfico SIN dibujar, para editar layout (posiciones)
# ------------------------------------------------------------
qg <- semPlot::semPaths(
  cfa_bi,
  what = "path",
  whatLabels = "std",
  style = "lisrel",
  layout = "tree2",
  rotation = 2,              # izquierda -> derecha
  residuals = FALSE,
  intercepts = FALSE,
  thresholds = FALSE,
  nCharNodes = 0,
  
  # Etiquetas controladas
  nodeLabels = nodeLabels,
  
  # Óvalos para ítems y latentes
  shapeLat = "ellipse",
  shapeMan = "ellipse",
  
  # Tamaños (AJUSTA aquí si quieres más/menos)
  sizeLat = 12,              # tamaño de los factores
  sizeMan = 5.2,             # tamaño de los ítems (sube/baja)
  label.cex = 0.95,          # tamaño texto dentro de nodos (ítems y factores)
  
  # Cargas más pequeñas + redondeadas
  edge.label.cex = 0.65,     # <- baja esto para que NO se monten las cargas
  digits = 2,
  
  edge.color = "gray30",
  fade = FALSE,
  curvePivot = TRUE,
  mar = c(6, 6, 6, 6),
  
  DoNotPlot = TRUE
)

# ------------------------------------------------------------
# 3) Ajustar el layout: más espacio entre óvalos y líneas más cortas
# ------------------------------------------------------------
lay <- qg$layout

# (A) Acortar líneas: acercar columnas (reduce distancia en X)
lay[, 1] <- lay[, 1] * 0.6   # prueba 0.70–0.85 según tu gusto

# (B) Separar óvalos de ítems en Y (evita que se monten)
man_idx <- which(spm@Vars$manifest)   # índices de variables observadas
ord <- man_idx[order(lay[man_idx, 2], decreasing = TRUE)]

# Re-espaciar de forma uniforme (más separación)
lay[ord, 2] <- seq(from =  3, to = -3, length.out = length(ord))

qg$layout <- lay

# ------------------------------------------------------------
# 4) Exportar en alta resolución (recomendado para paper)
# ------------------------------------------------------------
png(file.path(out_dir, "figure1_cfa_path_final.png"),
    width = 3200, height = 2000, res = 350)
plot(qg)
dev.off()

print(out_dir)

pdf(file.path(out_dir, "figure1_cfa_path_final.pdf"),
    width = 12, height = 7.5)
plot(qg)
dev.off()



crear_path_diagram <- function(modelo, item_map, outfile, sin_cargas = FALSE) {
  require(semPlot)
  require(dplyr)
  
  # 1. Preparar labels
  spm <- semPlot::semPlotModel(modelo)
  nodeLabels <- spm@Vars$name
  nodeLabels <- ifelse(nodeLabels == "IF",  "Financial\nImpact",
                       ifelse(nodeLabels == "Ips", "Psychosocial\nImpact",
                              recode(nodeLabels, !!!item_map, .default = nodeLabels)))
  
  # 2. Crear objeto semPlot sin graficar
  qg <- semPaths(
    modelo,
    what = "path",
    whatLabels = ifelse(sin_cargas, "none", "std"),
    style = "lisrel",
    layout = "tree2",
    rotation = 2,
    residuals = FALSE,
    intercepts = FALSE,
    thresholds = FALSE,
    nCharNodes = 0,
    nodeLabels = nodeLabels,
    shapeLat = "ellipse",
    shapeMan = "rectangle",
    sizeLat = 12,
    sizeMan = 6,
    label.cex = 0.9,
    edge.label.cex = ifelse(sin_cargas, 0, 0.7),
    edge.color = "gray30",
    fade = FALSE,
    mar = c(7, 7, 7, 7),
    curvePivot = TRUE,
    digits = 2,
    DoNotPlot = TRUE
  )
  
  # 3. Ajuste de layout
  lay <- qg$layout
  lay[, 1] <- lay[, 1] * 0.75
  man_idx <- which(spm@Vars$manifest)
  ord <- man_idx[order(lay[man_idx, 2], decreasing = TRUE)]
  lay[ord, 2] <- seq(from = 4.5, to = -4.5, length.out = length(ord))
  qg$layout <- lay
  
  # 4. Exportar PNG + PDF
  png(paste0(outfile, ".png"), width = 3600, height = 2200, res = 350, type = "cairo")
  par(mar = c(7, 7, 7, 7), xpd = NA)
  plot(qg)
  dev.off()
  
  pdf(paste0(outfile, ".pdf"), width = 12, height = 7.5)
  par(mar = c(7, 7, 7, 7), xpd = NA)
  plot(qg)
  dev.off()
  
  message("Exportado a: ", outfile, ".png/pdf")
}

# Mapeo de ítems
item_map <- setNames(c(paste0("IF", 1:10), paste0("IPS", 1:10)),
                     c(paste0("imf", sprintf("%02d",1:10), "_if"),
                       paste0("imps", sprintf("%02d",1:10), "_ips")))

# Crear path diagram con cargas
crear_path_diagram(cfa_bi, item_map, "outputs/Figure1_CFA_with_loadings", sin_cargas = FALSE)

# Versión sin números
crear_path_diagram(cfa_bi, item_map, "outputs/Figure1_CFA_clean", sin_cargas = TRUE)




message("Figura guardada en: ", normalizePath("outputs/Figure_1_CFA_clean.pdf"))



# ------------------------------------------------------------
# 0) Labels bonitos
# ------------------------------------------------------------
item_map <- c(paste0("IF", 1:10), paste0("IPS", 1:10))
names(item_map) <- c(paste0("imf", sprintf("%02d", 1:10), "_if"),
                     paste0("imps", sprintf("%02d", 1:10), "_ips"))

spm <- semPlot::semPlotModel(cfa_bi)

nodeLabels <- spm@Vars$name
nodeLabels <- ifelse(nodeLabels == "IF",  "Financial\nImpact",
                     ifelse(nodeLabels == "Ips", "Psychosocial\nImpact",
                            dplyr::recode(nodeLabels, !!!item_map, .default = nodeLabels)))

# ------------------------------------------------------------
# 1) Crear objeto sin plot
# ------------------------------------------------------------
qg <- semPlot::semPaths(
  cfa_bi,
  what = "path",
  whatLabels = "std",
  style = "lisrel",
  layout = "tree2",
  rotation = 2,
  residuals = FALSE,
  intercepts = FALSE,
  thresholds = FALSE,
  nCharNodes = 0,
  
  nodeLabels = nodeLabels,
  
  shapeLat = "ellipse",
  shapeMan = "ellipse",
  
  sizeLat = 13,        # <- un poco más grande para que los 2 factores se vean parejos
  sizeMan = 5.6,       # <- sube un poco el tamaño de ítems
  label.cex = 0.95,
  
  edge.label.cex = 0.60,  # <- cargas más pequeñas
  digits = 2,
  
  edge.color = "gray30",
  fade = FALSE,
  curvePivot = TRUE,
  mar = c(7, 7, 7, 7),
  
  DoNotPlot = TRUE
)

# ------------------------------------------------------------
# 2) Ajustar layout
# ------------------------------------------------------------
lay <- qg$layout

# (A) Acortar líneas: acercar columnas (reduce distancia en X)
lay[, 1] <- lay[, 1] * 0.85    # <- prueba 0.65–0.85

# (B) Separar óvalos de ítems en Y
man_idx <- which(spm@Vars$manifest)
ord <- man_idx[order(lay[man_idx, 2], decreasing = TRUE)]

# Rango más seguro (evita que se salga del panel al exportar)
lay[ord, 2] <- seq(from = 4, to = -4, length.out = length(ord))

qg$layout <- lay

# ------------------------------------------------------------
# 3) Exportar (CLAVE: par(xpd=NA) + márgenes)
# ------------------------------------------------------------
png(file.path(out_dir, "figure1_cfa_path_final.png"),
    width = 4000, height = 1330, res = 350)
par(mar = c(7, 7, 7, 7), xpd = NA)
plot(qg)
dev.off()

pdf(file.path(out_dir, "figure1_cfa_path_final.pdf"),
    width = 12.5, height = 8.5)
par(mar = c(7, 7, 7, 7), xpd = NA)
plot(qg)
dev.off()










# ==============================
# Figure 1 (CFA path diagram) - clean, journal-style (NO loadings)
# ==============================
# Requiere: semPlot, dplyr (ya los tienes)

library(semPlot)
library(dplyr)

# 1) Mapeo de etiquetas: items cortos + factores con nombre
item_map <- c(paste0("IF", 1:10), paste0("IPS", 1:10))
names(item_map) <- c(paste0("imf", sprintf("%02d", 1:10), "_if"),
                     paste0("imps", sprintf("%02d", 1:10), "_ips"))

spm <- semPlot::semPlotModel(cfa_bi)

nodeLabels <- spm@Vars$name
nodeLabels <- ifelse(nodeLabels == "IF",  "Financial\nImpact",
                     ifelse(nodeLabels == "Ips", "Psychosocial\nImpact",
                            dplyr::recode(nodeLabels, !!!item_map, .default = nodeLabels)))

# 2) Generar el objeto del gráfico SIN cargas y SIN dibujar (para export limpio)
qg <- semPlot::semPaths(
  cfa_bi,
  what        = "path",
  whatLabels  = "none",       # <-- CLAVE: sin números
  style       = "lisrel",
  layout      = "tree2",
  rotation    = 2,            # izquierda -> derecha (limpio)
  residuals   = FALSE,
  intercepts  = FALSE,
  thresholds  = FALSE,
  nCharNodes  = 0,
  nodeLabels  = nodeLabels,
  
  # formas (paper-like)
  shapeLat    = "ellipse",
  shapeMan    = "rectangle",  # <-- recomendado: evita que se monten como óvalos
  
  # tamaños (ajusta suave)
  sizeLat     = 12,
  sizeMan     = 6,
  label.cex   = 0.95,
  
  edge.color  = "gray25",
  fade        = FALSE,
  mar         = c(6, 6, 6, 6),
  
  DoNotPlot   = TRUE
)

# 3) Ajuste opcional del layout: acercar columnas y ordenar items
lay <- qg$layout

# Acorta líneas (acerca items a factores)
lay[, 1] <- lay[, 1] * 0.80   # prueba 0.70–0.90

# Aumenta separación vertical entre ítems para evitar solapes
man_idx <- which(spm@Vars$manifest)
ord <- man_idx[order(lay[man_idx, 2], decreasing = TRUE)]
lay[ord, 2] <- seq(from = 4.5, to = -4.5, length.out = length(ord))

qg$layout <- lay

# 4) EXPORTAR (Cairo evita “se exporta pero no se ve / se corta”)
# PNG alta resolución
png(file.path(out_dir, "Figure_1_CFA_clean.png"),
    width = 3600, height = 2200, res = 350, type = "cairo")
plot(qg)
dev.off()

# PDF vectorial nítido
pdf(file.path(out_dir, "Figure_1_CFA_clean.pdf"),
    width = 12, height = 7.5)
plot(qg)
dev.off()


# =========================================================
# 8) IRT (GRM) – TF-Col (ordinal, 5 categories)
# Outputs:
#  - Table 4: IRT summary ranges
#  - Figure 1: Test information + standard error curves
#  - Supplementary Table S2: full item parameters
# =========================================================
install.packages("mirt")
library(mirt)
install.packages("ggplot2")
library(ggplot2)
install.packages("tidyr")
library(tidyr)
install.packages("purrr")
library(purrr)

# ---------------------------
# 8.1 Prepare data for IRT
# ---------------------------
# Convert ordered factors to numeric 1..K
df_num <- df %>%
  mutate(across(everything(), ~ as.numeric(.x)))  # assumes levels 1..5

# Keep complete cases for IRT calibration (recommended)
df_irt <- df_num[complete.cases(df_num), ]
n_irt <- nrow(df_irt)
writeLines(paste0("IRT complete-case N = ", n_irt),
           con = file.path(out_dir, "irt_sample_size.txt"))

# Sanity checks
stopifnot(all(sapply(df_irt, is.numeric)))
stopifnot(all(sapply(df_irt, function(x) all(x %in% 1:5))))

# ---------------------------
# 8.2 Fit GRM (1-factor)
# ---------------------------
# Technical settings for stable estimation
grm_1f <- mirt(df_irt, 1, itemtype = "graded",
               technical = list(NCYCLES = 2000),
               verbose = FALSE)

saveRDS(grm_1f, file.path(out_dir, "irt_grm_1factor_object.rds"))

# Convergence note
writeLines(c(
  paste0("Converged: ", extract.mirt(grm_1f, "converged")),
  paste0("LogLik: ", extract.mirt(grm_1f, "logLik"))
), con = file.path(out_dir, "irt_model_notes.txt"))

# ---------------------------
# 8.3 Extract full parameters (Supplementary Table S2)
# ---------------------------
# IRTpars=TRUE returns a (a, b1..b4) parameterization for GRM
pars <- coef(grm_1f, IRTpars = TRUE, simplify = TRUE)$items
pars_df <- as.data.frame(pars) %>%
  tibble::rownames_to_column("item")

# Keep only a and b thresholds (names may be a1, b1..b4 depending on version)
# mirt usually names discrimination as 'a1' and thresholds as 'b1','b2',...
names(pars_df) <- gsub("\\.", "_", names(pars_df))

# Export full parameters
readr::write_csv(pars_df, file.path(out_dir, "supp_table_s2_irt_full_parameters.csv"))

# ---------------------------
# 8.4 Summary ranges (Table 4 in main manuscript)
# ---------------------------
# Identify columns
a_col <- intersect(names(pars_df), c("a1","a"))
b_cols <- names(pars_df)[stringr::str_detect(names(pars_df), "^b\\d+")]

stopifnot(length(a_col) == 1, length(b_cols) >= 2)

a_vals <- pars_df[[a_col]]
b_vals <- pars_df[, b_cols, drop = FALSE]

irt_summary <- tibble::tibble(
  n_irt = n_irt,
  discrimination_min = min(a_vals, na.rm = TRUE),
  discrimination_max = max(a_vals, na.rm = TRUE),
  discrimination_median = median(a_vals, na.rm = TRUE),
  threshold_min = min(as.matrix(b_vals), na.rm = TRUE),
  threshold_max = max(as.matrix(b_vals), na.rm = TRUE)
)

readr::write_csv(irt_summary, file.path(out_dir, "table4_irt_parameter_summary.csv"))

# (Optional) ranges by threshold (nice for reviewers)
irt_b_by_k <- purrr::map_dfr(b_cols, \(cc) {
  tibble::tibble(
    threshold = cc,
    min = min(b_vals[[cc]], na.rm = TRUE),
    max = max(b_vals[[cc]], na.rm = TRUE)
  )
})
readr::write_csv(irt_b_by_k, file.path(out_dir, "irt_threshold_ranges_by_category.csv"))

# ---------------------------
# 8.5 Figure 1: Test Information + Standard Error curves
# ---------------------------
theta <- seq(-4, 4, by = 0.05)

# testinfo() gives information across theta
info <- testinfo(grm_1f, Theta = theta)
se <- 1 / sqrt(info)

fig_df <- tibble::tibble(theta = theta, information = info, se = se)

# Long format for plotting
fig_long <- fig_df %>%
  tidyr::pivot_longer(cols = c("information","se"),
                      names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         information = "Test information",
                         se = "Standard error"))

p <- ggplot(fig_long, aes(x = theta, y = value)) +
  geom_line(linewidth = 1.0, color = "gray20") +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  labs(x = expression(theta), y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "figure1_test_information_standard_error.png"),
       p, width = 6.5, height = 7.0, dpi = 300)

ggsave(file.path(out_dir, "figure1_test_information_standard_error.pdf"),
       p, width = 6.5, height = 7.0)

# =========================================================
# 8.6 Item fit (S-X² / G²) – "¿el ítem se comporta como espera el modelo?"
# Exporta:
#  - outputs/irt_itemfit_sx2_g2.csv  (tabla completa por ítem)
#  - outputs/irt_itemfit_summary.csv (resumen: media/rango)
# =========================================================

library(mirt)
library(dplyr)
library(readr)
library(tibble)

# Asegurar que el objeto exista
stopifnot(exists("grm_1f"))

fit_it <- itemfit(grm_1f)  # S-X2 / G2 según corresponda a tu versión de mirt

fit_it

fit_df <- as.data.frame(fit_it) %>%
  rename(pregunta = item)

# Export tabla completa
write_csv(fit_df, file.path(out_dir, "irt_itemfit_sx2_g2.csv"))

# Resumen útil para texto (solo columnas numéricas)
fit_summary <- fit_df %>%
  summarise(across(where(is.numeric),
                   list(mean = ~mean(.x, na.rm=TRUE),
                        min  = ~min(.x, na.rm=TRUE),
                        max  = ~max(.x, na.rm=TRUE)),
                   .names = "{.col}_{.fn}"))

write_csv(fit_summary, file.path(out_dir, "irt_itemfit_summary.csv"))

cat("\nItem fit listo:\n",
    "- outputs/irt_itemfit_sx2_g2.csv\n",
    "- outputs/irt_itemfit_summary.csv\n")

#TAbla para articulo (ordenada)

library(dplyr)
library(readr)
library(tibble)

# Convertir a data.frame (ya trae columna 'item')
fit_df <- as.data.frame(fit_it)

# Ordenar y redondear (opcional, recomendado para suplemento)
fit_df_clean <- fit_df %>%
  mutate(
    S_X2 = round(S_X2, 3),
    RMSEA_S_X2 = round(RMSEA.S_X2, 3),
    p_S_X2 = round(p.S_X2, 3)
  ) %>%
  select(
    item,
    S_X2,
    df.S_X2,
    RMSEA_S_X2,
    p_S_X2
  ) %>%
  arrange(item)

# Exportar tabla suplementaria
write_csv(
  fit_df_clean,
  file.path(out_dir, "supp_table_s4_item_fit_sx2.csv")
)

# Vista rápida
print(fit_df_clean)

# ============================
# Wright map / Person-Item Map
# ============================

install.packages("WrightMap")
library(WrightMap)
library(mirt)
library(dplyr)
library(tibble)

# theta (habilidad/persona)
theta_hat <- fscores(grm_1f, method = "EAP")[,1]

# parámetros IRT (a, b1-b4)
pars <- coef(grm_1f, IRTpars = TRUE, simplify = TRUE)$items %>%
  as.data.frame() %>%
  rownames_to_column("item")

# quedarnos con los thresholds b1..b4
b_cols <- grep("^b", names(pars), value = TRUE)
thr_mat <- as.matrix(pars[, b_cols])
rownames(thr_mat) <- pars$item


theta_df <- data.frame(
  theta = fscores(grm_1f, method = "EAP")[,1]
)


pars <- coef(grm_1f, IRTpars = TRUE, simplify = TRUE)$items %>%
  as.data.frame() %>%
  tibble::rownames_to_column("item")

b_cols <- grep("^b", names(pars), value = TRUE)

item_df <- pars %>%
  select(item, all_of(b_cols)) %>%
  tidyr::pivot_longer(
    cols = starts_with("b"),
    names_to = "threshold",
    values_to = "location"
  )



# ======================================
# Publication-ready Wright map (GRM) - mirt
# ======================================

library(mirt)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(readr)

stopifnot(exists("grm_1f"))
if (!exists("out_dir")) out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir)

# ---------------------------
# 1) Person locations (theta)
# ---------------------------
theta_hat <- as.numeric(fscores(grm_1f, method = "EAP")[,1])

# ---------------------------
# 2) Item thresholds (b1..b4)
# ---------------------------
pars_items <- coef(grm_1f, IRTpars = TRUE, simplify = TRUE)$items
thr_cols <- grep("^b", colnames(pars_items), value = TRUE)
thr_mat <- as.data.frame(pars_items[, thr_cols, drop = FALSE]) %>%
  tibble::rownames_to_column("item_raw")

# 3) Short labels
thr_mat <- thr_mat %>%
  mutate(
    item = item_raw,
    item = gsub("^imf0?([0-9]+)_if$", "IF\\1", item),
    item = gsub("^imps0?([0-9]+)_ips$", "IPS\\1", item)
  )

# Order items: IF1..IF10 then IPS1..IPS10
ord <- c(paste0("IF", 1:10), paste0("IPS", 1:10))
thr_mat$item <- factor(thr_mat$item, levels = rev(ord))  # rev -> top-to-bottom

# Long format: one row per threshold
thr_long <- thr_mat %>%
  pivot_longer(cols = all_of(thr_cols), names_to = "threshold", values_to = "b") %>%
  mutate(threshold = factor(threshold, levels = thr_cols))

# ---------------------------
# 4) Plot settings
# ---------------------------
x_limits <- range(c(theta_hat, thr_long$b), na.rm = TRUE)
pad <- 0.3
x_limits <- c(x_limits[1] - pad, x_limits[2] + pad)

# ---------------------------
# 5) Left panel: persons (hist)
# ---------------------------
p_person <- ggplot(data.frame(theta = theta_hat), aes(x = theta)) +
  geom_histogram(binwidth = 0.25, fill = "grey70", color = "grey35") +
  coord_flip() +
  scale_x_continuous(limits = x_limits) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(x = expression(theta), y = "Respondents")

# ---------------------------
# 6) Right panel: item thresholds
# ---------------------------
p_items <- ggplot(thr_long, aes(x = b, y = item)) +
  geom_point(aes(shape = threshold), size = 2.2, color = "grey15") +
  scale_shape_manual(values = c(16, 17, 15, 18)) +  # b1..b4 different shapes
  scale_x_continuous(limits = x_limits) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  labs(x = expression(theta), y = NULL) +
  guides(shape = guide_legend(nrow = 1))

# ---------------------------
# 7) Combine panels (no extra packages)
#    We'll save separately OR use patchwork if available
# ---------------------------
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
library(patchwork)

p_wright <- p_person + p_items +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(
    title = "Person–Item Map (TF-Col, Graded Response Model)",
    subtitle = "Points represent category thresholds (b1–b4) for each item; left panel shows respondent trait distribution"
  )

# ---------------------------
# 8) Export (publication-ready)
# ---------------------------
ggsave(file.path(out_dir, "figure_wright_map_tfcol_grm.png"),
       p_wright, width = 11, height = 6.5, dpi = 300)

ggsave(file.path(out_dir, "figure_wright_map_tfcol_grm.pdf"),
       p_wright, width = 11, height = 6.5)


################################
#Validez de criterio y constructo
################################

library(dplyr)
library(stringr)
library(readr)
library(tibble)

# =========================
# 0) Variables (según tus nombres reales)
# =========================
tf_total <- "puntaje_tf"
tf_if    <- "puntaje_if"
tf_ips   <- "puntaje_ips"

cost_total <- "puntaje_cost"

peds_total <- "score_total_ql"
peds_psych <- "score_psico_social_ql"

group_var  <- "tipo_de_poblacion"   # Adulto / Infantil (en clean_names)

# =========================
# 1) Forzar a numérico (robusto para Excel)
# =========================
num_vars <- c(tf_total, tf_if, tf_ips, cost_total, peds_total, peds_psych)

datos_num <- datos_clean %>%
  mutate(across(all_of(num_vars), ~{
    x <- as.character(.x)
    x <- str_trim(x)
    x <- na_if(x, "")
    x <- str_replace_all(x, ",", ".")
    suppressWarnings(as.numeric(x))
  }))

# =========================
# 2) Función Spearman + IC 95% (Fisher z approx)
# =========================
spearman_ci_fisher <- function(data, x, y, conf = 0.95) {
  d <- data %>%
    select(all_of(c(x, y))) %>%
    filter(complete.cases(.))
  
  n <- nrow(d)
  if (n < 10) {
    return(tibble(x=x, y=y, n=n, rho=NA_real_, ci_low=NA_real_, ci_high=NA_real_))
  }
  
  r <- suppressWarnings(cor(d[[x]], d[[y]], method = "spearman"))
  r <- max(min(r, 0.9999), -0.9999)
  
  z <- atanh(r)
  se <- 1/sqrt(n - 3)
  alpha <- 1 - conf
  zcrit <- qnorm(1 - alpha/2)
  
  tibble(
    x = x, y = y, n = n,
    rho = r,
    ci_low = tanh(z - zcrit*se),
    ci_high = tanh(z + zcrit*se)
  )
}

# =========================
# 3) Definir subconjuntos correctos
# =========================
data_all <- datos_num

# OJO: aquí defines el subconjunto PedsQL.
# Si "Infantil" corresponde a cuidadores/pediatría, esto es lo correcto:
data_peds <- datos_num %>%
  filter(.data[[group_var]] == "Infantil")

# =========================
# 4) Recalcular pares: COST en total / PedsQL en infantil
# =========================
pairs_cost <- tribble(
  ~tf_score, ~comp,        ~dataset,
  tf_total,  cost_total,   "all",
  tf_if,     cost_total,   "all",
  tf_ips,    cost_total,   "all"
)

pairs_peds <- tribble(
  ~tf_score, ~comp,        ~dataset,
  tf_total,  peds_total,   "peds",
  tf_if,     peds_total,   "peds",
  tf_ips,    peds_total,   "peds",
  tf_total,  peds_psych,   "peds",
  tf_if,     peds_psych,   "peds",
  tf_ips,    peds_psych,   "peds"
)

pairs_all <- bind_rows(pairs_cost, pairs_peds)

pretty <- c(
  puntaje_tf = "TF-Col total score",
  puntaje_if = "Financial Impact",
  puntaje_ips = "Psychosocial Impact",
  puntaje_cost = "COST-FACIT total score",
  score_total_ql = "PedsQL total score",
  score_psico_social_ql = "PedsQL psychosocial summary"
)

results_table6 <- pairs_all %>%
  rowwise() %>%
  mutate(res = list(
    spearman_ci_fisher(
      data = ifelse(dataset == "peds", list(data_peds), list(data_all))[[1]],
      x = tf_score,
      y = comp
    )
  )) %>%
  tidyr::unnest(res) %>%
  ungroup() %>%
  mutate(
    TF_Col = recode(x, !!!pretty),
    Comparator = recode(y, !!!pretty),
    rho_ci = sprintf("%.3f (%.3f to %.3f)", rho, ci_low, ci_high)
  ) %>%
  select(TF_Col, Comparator, n, rho, ci_low, ci_high, rho_ci)

print(results_table6)

# Exporta para comparar con tu tabla actual
write_csv(results_table6, file.path(out_dir, "table6_correlations_main_verified.csv"))




############################################################
# Supplementary Table S3 (UPDATED):
# TF-Col (total + domains) vs PedsQL domains (pediatric only)
# Spearman rho + 95% CI (Fisher z)
############################################################

library(dplyr)
library(stringr)
library(purrr)
library(readr)
library(tibble)

# ---------------------------
# 0) REAL variable names in datos_clean
# ---------------------------
tf_total <- "puntaje_tf"
tf_if    <- "puntaje_if"
tf_ips   <- "puntaje_ips"

# PedsQL domains (caregiver-reported pediatric subsample)
peds_domains <- c("df_ql","de_ql","ds_ql","des_ql")

# Pretty labels (for table readability)
pretty <- c(
  puntaje_tf = "TF-Col total score",
  puntaje_if = "Financial Impact",
  puntaje_ips = "Psychosocial Impact",
  df_ql = "PedsQL physical functioning",
  de_ql = "PedsQL emotional functioning",
  ds_ql = "PedsQL social functioning",
  des_ql = "PedsQL school functioning"
)

# ---------------------------
# 1) Force numeric (robust for Excel imports)
# ---------------------------
num_vars <- c(tf_total, tf_if, tf_ips, peds_domains)

datos_num <- datos_clean %>%
  mutate(across(all_of(num_vars), ~ {
    x <- as.character(.x)
    x <- str_trim(x)
    x <- na_if(x, "")
    x <- str_replace_all(x, ",", ".")   # decimal comma -> dot
    suppressWarnings(as.numeric(x))
  }))

# ---------------------------
# 2) Optional: restrict to pediatric/caregiver subsample
#    (Recommended, so n matches pediatric PedsQL availability)
#    Adjust this filter if your coding differs.
# ---------------------------
# If your variable is "tipo_de_poblacion" and pediatric rows are "Infantil":
if ("tipo_de_poblacion" %in% names(datos_num)) {
  datos_peds <- datos_num %>%
    filter(str_to_lower(tipo_de_poblacion) == "infantil")
} else {
  # fallback: don't filter if variable not found
  datos_peds <- datos_num
}

# ---------------------------
# 3) Helper: Spearman rho + Fisher z CI (approx)
# ---------------------------
spearman_ci_fisher <- function(data, x, y, conf = 0.95) {
  d <- data %>%
    select(all_of(c(x, y))) %>%
    filter(complete.cases(.))
  
  n <- nrow(d)
  if (n < 10) {
    return(tibble(x=x, y=y, n=n, rho=NA_real_, ci_low=NA_real_, ci_high=NA_real_))
  }
  
  if (!is.numeric(d[[x]]) || !is.numeric(d[[y]])) {
    stop(sprintf("Non-numeric detected in x=%s or y=%s after conversion.", x, y))
  }
  
  r <- suppressWarnings(cor(d[[x]], d[[y]], method = "spearman"))
  r <- max(min(r, 0.9999), -0.9999)  # avoid atanh(±1)
  
  z <- atanh(r)
  se <- 1/sqrt(n - 3)
  alpha <- 1 - conf
  zcrit <- qnorm(1 - alpha/2)
  
  z_low <- z - zcrit * se
  z_high <- z + zcrit * se
  
  tibble(
    x = x, y = y, n = n,
    rho = r,
    ci_low = tanh(z_low),
    ci_high = tanh(z_high)
  )
}

# ---------------------------
# 4) Build UPDATED S3 pairs:
#    INCLUDE TF-Col total + IF + IPS vs each PedsQL domain
# ---------------------------
pairs_s3 <- expand.grid(
  tf_score = c(tf_total, tf_if, tf_ips),
  comparator = peds_domains,
  stringsAsFactors = FALSE
) %>% as_tibble()

supp_s3_updated <- purrr::pmap_dfr(
  pairs_s3,
  ~ spearman_ci_fisher(datos_peds, ..1, ..2)
) %>%
  mutate(
    TF_Col = recode(x, !!!pretty),
    Comparator = recode(y, !!!pretty),
    rho_ci = sprintf("%.3f (%.3f to %.3f)", rho, ci_low, ci_high)
  ) %>%
  select(TF_Col, Comparator, n, rho, ci_low, ci_high, rho_ci) %>%
  # Nice ordering: Total first, then IF, then IPS; within each: Physical, Emotional, Social, School
  mutate(
    TF_Col = factor(TF_Col, levels = c("TF-Col total score","Financial Impact","Psychosocial Impact")),
    Comparator = factor(Comparator, levels = c(
      "PedsQL physical functioning",
      "PedsQL emotional functioning",
      "PedsQL social functioning",
      "PedsQL school functioning"
    ))
  ) %>%
  arrange(TF_Col, Comparator) %>%
  mutate(TF_Col = as.character(TF_Col),
         Comparator = as.character(Comparator))

# ---------------------------
# 5) Export
# ---------------------------
out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir)

write_csv(supp_s3_updated, file.path(out_dir, "supp_table_s3_pedsql_domains_UPDATED.csv"))

cat("\n=== Supplementary Table S3 (UPDATED) ===\n")
print(supp_s3_updated)


