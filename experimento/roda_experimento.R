# =============================================================================
# EXECUCAO DOS EXPERIMENTOS - IC M4 HOURLY
# Um unico script para todos os modelos (originais e novos) e variantes da revisao.
#
# Terminal:
#   Rscript roda_experimento.R <MODELO> <H> [VARIANTES] [PARTE]
#     MODELO     ARIMA | ARIMA_FOURIER | SARIMA | MLP | ELM   (LSTM desativado: ver abaixo)
#     H          horizonte de previsao (3, 5, 10)
#     VARIANTES  "base" ou combinacao com "+":
#                  S24   redes com a entrada sazonal x[t-24]            -> prefixo <MODELO>_S24
#                  S24W  redes com x[t-24] e x[t-168] (dia e semana)    -> prefixo <MODELO>_S24W
#                  CT    alvo contaminado no treino/validacao das redes -> prefixo <MODELO>_CT
#                  FACT  + cenarios so-gaussiano e so-outlier (10 cenarios no total)
#                  REP3 / REP5   3 ou 5 realizacoes da contaminacao
#                  N100 / N200 / N414   100, 200 ou todas as series do M4 Hourly
#                (FACT, REPx e Nx nao mudam o prefixo: o arquivo inclui os resultados originais)
#     PARTE      "k/n" roda so a k-esima de n partes das series (rodada paralela)
#                "junta" so junta os checkpoints e grava o Excel
#   Exemplos:
#     Rscript roda_experimento.R ARIMA_FOURIER 5
#     Rscript roda_experimento.R MLP 3 S24
#     Rscript roda_experimento.R ELM 10 CT+FACT 2/4
#     Rscript roda_experimento.R ELM 10 CT+FACT junta
#   Em paralelo (Linux/macOS): bash roda_paralelo.sh MLP 3 S24 4
#
# RStudio: edite a secao "SEM ARGUMENTOS" abaixo e rode o arquivo inteiro (Source).
#
# Teste rapido (confere se tudo roda em ~1-2 min antes de uma rodada longa):
#   IC_TESTE_RAPIDO=1 Rscript roda_experimento.R MLP 5 S24
# (2 series, 2 janelas, 3 filtros; grava com prefixo TESTE_<...>, pode apagar depois)
# =============================================================================

# ---------- SEM ARGUMENTOS (RStudio) ----------
MODELO <- "ARIMA_FOURIER"
H <- 5L
VARIANTES <- "base"
PARTE <- NA_character_          # ex.: "1/4", ou "junta"

# ---------- Grades de hiperparametros ----------
# CONFIRA com os seus scripts originais (MLP.R, ELM.R) e copie a grade usada.
# ELM: valores selecionados nos seus resultados (nhid = 5, 10, 20, 50).
GRIDS <- list(
  ARIMA = data.frame(dummy = 1), ARIMA_FOURIER = data.frame(dummy = 1), SARIMA = data.frame(dummy = 1),
  ELM = data.frame(nhid = c(5L, 10L, 20L, 50L)),
  MLP = expand.grid(size = c(5L, 10L), decay = c(0.001, 0.01), maxit = 300L)       # <- COPIE A SUA GRADE
)
# LSTM desativado (nao cabe no prazo). Para reativar, acrescente a GRIDS:
#   LSTM = data.frame(units = c(16L, 32L), dropout = c(0.1, 0.1), learning_rate = c(0.001, 0.001),
#                     epochs = c(60L, 120L), batch_size = c(32L, 32L))

# ---------- Argumentos ----------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) MODELO <- args[1]
if (length(args) >= 2L) H <- as.integer(args[2])
if (length(args) >= 3L) VARIANTES <- args[3]
if (length(args) >= 4L) PARTE <- args[4]
if (toupper(MODELO) == "LSTM") stop("LSTM desativado neste script (nao cabe no prazo). Veja o comentario em GRIDS para reativar.")
if (!MODELO %in% names(GRIDS)) stop("MODELO deve ser um de: ", paste(names(GRIDS), collapse = ", "))
if (!is.finite(H) || H < 1L) stop("H deve ser um inteiro positivo.")
vars <- toupper(strsplit(VARIANTES, "+", fixed = TRUE)[[1]])
vars <- setdiff(vars, "BASE")
validas <- c("S24", "S24W", "CT", "FACT", "REP3", "REP5", "N100", "N200", "N414")
if (length(setdiff(vars, validas))) stop("Variantes desconhecidas: ", paste(setdiff(vars, validas), collapse = ", "))
redes <- MODELO %in% c("MLP", "ELM")
if (!redes && any(vars %in% c("S24", "S24W", "CT")))
  stop("S24, S24W e CT so se aplicam a MLP e ELM (os modelos ARIMA ja treinam so com a serie contaminada).")
if (all(c("S24", "S24W") %in% vars)) stop("Use S24 ou S24W, nao os dois.")

source("IC_Common.R")

# ---------- Configuracao da rodada ----------
PRED_HORIZON <- H
prefixo <- MODELO
if ("S24" %in% vars)  { NN_SEASONAL_LAGS <- 24L;           prefixo <- paste0(prefixo, "_S24") }
if ("S24W" %in% vars) { NN_SEASONAL_LAGS <- c(24L, 168L);  prefixo <- paste0(prefixo, "_S24W") }
if ("CT" %in% vars)   { TARGET_MODE <- "corrupted";        prefixo <- paste0(prefixo, "_CT") }
if ("FACT" %in% vars) SCENARIO_SET <- "factorial"
if ("REP3" %in% vars) N_REPLICATIONS <- 3L
if ("REP5" %in% vars) N_REPLICATIONS <- 5L
if ("N100" %in% vars) N_SERIES <- 100L
if ("N200" %in% vars) N_SERIES <- 200L
if ("N414" %in% vars) N_SERIES <- 414L

if (nzchar(Sys.getenv("IC_TESTE_RAPIDO"))) {
  message(">>> TESTE RAPIDO: 2 series, 2 janelas, 3 filtros")
  N_WINDOWS <- 2L; N_VALIDATION_ORIGINS <- 1L; filter_names <- filter_names[1:3]
  N_SERIES <- 2L; prefixo <- paste0("TESTE_", prefixo)
  if (identical(SCENARIO_SET, "factorial")) factorial_scenarios <- factorial_scenarios[c("Gauss_High", "Out_High")]
}

grid <- GRIDS[[MODELO]]
message(sprintf("Modelo %s | h = %d | variantes: %s | prefixo: %s | checkpoint: %s",
                MODELO, H, if (length(vars)) paste(vars, collapse = "+") else "base", prefixo,
                checkpoint_dir_for(prefixo)))

# ---------- Execucao ----------
if (!is.na(PARTE) && tolower(PARTE) == "junta") {
  collect_experiment(prefixo)
} else if (!is.na(PARTE)) {
  kn <- as.integer(strsplit(PARTE, "/", fixed = TRUE)[[1]])
  if (length(kn) != 2L || any(is.na(kn)) || kn[1] < 1L || kn[1] > kn[2]) stop("PARTE deve ser 'k/n', ex.: 2/4")
  idx <- split(seq_len(N_SERIES), cut(seq_len(N_SERIES), kn[2], labels = FALSE))[[kn[1]]]
  message(sprintf("Parte %d de %d: series %d a %d", kn[1], kn[2], min(idx), max(idx)))
  run_experiment(MODELO, grid, prefixo, series_idx = idx, collect = FALSE)
  message("Parte concluida. Quando todas terminarem: Rscript roda_experimento.R ", MODELO, " ", H, " ",
          VARIANTES, " junta")
} else {
  run_experiment(MODELO, grid, prefixo)
}
