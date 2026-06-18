# ============================================================
# TABLE 1 — Sample characteristics of the study population
# Manuscrito: VIHRI-LA-2026-0036 (TF-Col)
# Genera la Tabla 1 de forma reproducible desde el dataset.
# ============================================================
# INSTRUCCIONES:
# 1. Coloca este script en la misma carpeta que el archivo Excel.
# 2. Ajusta ARCHIVO si el nombre difiere.
# 3. Corre el script completo. Produce Table1_sample_characteristics.csv
#    e imprime la tabla en consola.
# ============================================================

# 0. PAQUETES -------------------------------------------------
pkgs <- c("readxl", "dplyr")
nuevos <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(nuevos) > 0) install.packages(nuevos)
invisible(lapply(pkgs, library, character.only = TRUE))

DATA_FILE <- "data/tfcol_dataset.xlsx"   # primera hoja = datos (no incluido; ver data/README.md)

# 1. CARGA ----------------------------------------------------
df <- readxl::read_excel(DATA_FILE)

# 2. VARIABLES DE TRABAJO -------------------------------------
pop    <- df[["Tipo de población"]]
adult  <- pop == "Adulto"      # paciente adulto = respondiente
infant <- pop == "Infantil"    # cuidador respondiente; paciente = niño

edad_pac  <- suppressWarnings(as.numeric(df[["6_edad_paciente"]]))  # edad del paciente (adulto o niño)
edad_cuid <- suppressWarnings(as.numeric(df[["8_edad_cuidador"]]))  # edad del cuidador
sexo_pac  <- df[["5_sexo_paciente"]]
sexo_cuid <- df[["7_sexo_cuidador"]]
educ      <- df[["28_niveleducativo"]]
dx        <- df[["11_dx_cat_paciente"]]

N  <- nrow(df)             # 348
nA <- sum(adult)           # 200
nC <- sum(infant)          # 148 (cuidadores / pacientes pediátricos)

# 3. HELPERS DE FORMATO ---------------------------------------
fmt_ms <- function(x) {                 # media ± DE (DE muestral), 1 decimal
  x <- x[!is.na(x)]
  sprintf("%.1f \u00b1 %.1f", mean(x), stats::sd(x))
}
fmt_np <- function(n, d) {              # n (%) con 1 decimal
  if (n == 0) return("\u2014")
  sprintf("%d (%.1f)", n, 100 * n / d)
}
cnt <- function(x, key) sum(x == key, na.rm = TRUE)
DASH <- "\u2014"                         # em dash (—)

# 4. CONSTRUCCIÓN DE LA TABLA ---------------------------------
tab <- data.frame(Characteristic = character(), `Total (N = 348)` = character(),
                  `Adult patients (n = 200)` = character(),
                  `Caregivers (n = 148)` = character(),
                  check.names = FALSE, stringsAsFactors = FALSE)
add <- function(tab, a, b, c, d) rbind(tab, setNames(
  data.frame(a, b, c, d, stringsAsFactors = FALSE), names(tab)))

# Tipo de respondiente
tab <- add(tab, "Type of respondent, n (%)", "", "", "")
tab <- add(tab, "  Adult patient",                  fmt_np(nA, N), DASH, DASH)
tab <- add(tab, "  Caregiver of pediatric patient", fmt_np(nC, N), DASH, DASH)

# Edades: tres sujetos distintos, etiquetados explícitamente
tab <- add(tab, "Adult patient age, years, mean \u00b1 SD",     DASH, fmt_ms(edad_pac[adult]),  DASH)
tab <- add(tab, "Caregiver age, years, mean \u00b1 SD",         DASH, DASH, fmt_ms(edad_cuid))
tab <- add(tab, "Pediatric patient age, years, mean \u00b1 SD", DASH, DASH, fmt_ms(edad_pac[infant]))

# Sexo del respondiente
tab <- add(tab, "Sex of respondent, n (%)", "", "", "")
tab <- add(tab, "  Female", DASH, fmt_np(cnt(sexo_pac[adult], "mujer"),  nA), fmt_np(cnt(sexo_cuid, "mujer"),  nC))
tab <- add(tab, "  Male",   DASH, fmt_np(cnt(sexo_pac[adult], "hombre"), nA), fmt_np(cnt(sexo_cuid, "hombre"), nC))

# Nivel educativo (a nivel total)
tab <- add(tab, "Educational level, n (%)", "", "", "")
edu_map <- list(c("No formal education","Sin educación"), c("Primary education","Primaria"),
                c("Secondary education","Secundaria"), c("Technical/technologist","Técnico o tecnólogo"),
                c("Undergraduate degree","Pregrado"), c("Postgraduate degree","Posgrado"))
for (m in edu_map) tab <- add(tab, paste0("  ", m[1]), fmt_np(cnt(educ, m[2]), N), DASH, DASH)

# Diagnóstico de cáncer
tab <- add(tab, "Cancer diagnosis, n (%)", "", "", "")
dx_map <- list(c("Leukemia","Leucemia"), c("Thyroid cancer","Cáncer de tiroides"),
               c("Breast cancer","Cáncer de mama"), c("Lymphoma","Linfoma"))
principales <- vapply(dx_map, function(m) m[2], character(1))
for (m in dx_map) {
  tab <- add(tab, paste0("  ", m[1]),
             fmt_np(cnt(dx, m[2]), N),
             fmt_np(cnt(dx[adult], m[2]), nA),
             fmt_np(cnt(dx[infant], m[2]), nC))
}
otros   <- !(dx %in% principales)
tab <- add(tab, "  Other cancers*",
           fmt_np(sum(otros, na.rm = TRUE), N),
           fmt_np(sum(otros & adult, na.rm = TRUE), nA),
           fmt_np(sum(otros & infant, na.rm = TRUE), nC))

# 5. EXPORTAR E IMPRIMIR --------------------------------------
if(!dir.exists("outputs")) dir.create("outputs")
write.csv(tab, "outputs/Table1_sample_characteristics.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
print(tab, row.names = FALSE, right = FALSE)

# 6. AUTO-VERIFICACIÓN (valores reportados en el manuscrito) --
stopifnot(
  N == 348, nA == 200, nC == 148,
  fmt_ms(edad_pac[adult])  == "58.3 \u00b1 13.5",
  fmt_ms(edad_cuid)        == "36.1 \u00b1 9.5",
  fmt_ms(edad_pac[infant]) == "8.2 \u00b1 4.4"
)
cat("\nOK: cifras consistentes con el manuscrito.\n")
cat("Edad disponible: pacientes adultos", sum(!is.na(edad_pac[adult])), "/", nA,
    "| cuidadores", sum(!is.na(edad_cuid)), "/", nC,
    "| pacientes pediátricos", sum(!is.na(edad_pac[infant])), "/", nC, "\n")

# Versiones (para reportar en Methods, AE #4):
cat("\nR:", R.version.string, "\n")
cat("readxl:", as.character(packageVersion("readxl")),
    "| dplyr:", as.character(packageVersion("dplyr")), "\n")
