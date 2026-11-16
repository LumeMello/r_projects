# =============================================================================
# FIGURAS DO ARTIGO - IC M4 HOURLY / RESILIENCIA A CORRUPCAO
# Todos os textos das figuras estao em ingles (lista TXT abaixo).
# Rodar na mesma pasta dos resultados (<MODELO>_industrial_sensors_results_<N>pts.xlsx)
# e do IC_Common.R. Gera, em FIG_DIR, PNG (300 dpi) e PDF de:
#   0. Perfil de contaminacao (serie limpa x corrompida, 4 cenarios)
#   1. Diagramas de Diferenca Critica (Demsar / Nemenyi)
#   2. Heatmap de Delta RMSE% mediano vs. baseline (asterisco = Wilcoxon-Holm)
#   3. Fronteira de Pareto: acuracia vs. custo computacional
#   4. Painel diagnostico de degradacao: (a) Error Ratio por passo, (b) 48 h sobrepostas
#   5. Barras empilhadas do teste Diebold-Mariano serie a serie
# e as tabelas LaTeX de ranking de Friedman (uma por modelo) em TAB_DIR.
# Pacotes: ggplot2, patchwork, ggrepel, readxl (e os do IC_Common.R para o painel 4b
# e para medir o custo dos filtros quando nao houver tabela de overhead).
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

# ------------------------------- Configuracao --------------------------------
RESULTS_DIR <- "."
FIG_DIR <- "figuras"
# Modelos (prefixo do arquivo). Acrescente variantes rodadas com outro output_prefix,
# ex.: "MLP_CT" (alvo contaminado) ou "MLP_S24" (lags sazonais).
MODELS <- c("ARIMA", "ARIMA_FOURIER", "SARIMA", "MLP", "ELM", "LSTM")
N_WINDOWS_RESULTS <- 10L          # horizonte = pontos de teste / N_WINDOWS_RESULTS
ALPHA <- 0.05

# Overhead por filtro (ms). Se este CSV existir (colunas Filter, Overhead_ms, com os
# nomes de filtro do experimento: ma, ema, LOESS, ...), ele e usado (ex.: Tabela 7).
# Senao, o custo e medido aplicando cada filtro e salvo em OVERHEAD_MEASURED_CSV.
OVERHEAD_CSV <- "tabela7_overhead.csv"
OVERHEAD_MEASURED_CSV <- "overhead_medido.csv"
OVERHEAD_BENCH_SERIES <- 10L

# Capacidade de atualizacao incremental (edite se a sua classificacao for outra).
# O(1): o valor novo sai do estado anterior + a observacao nova (janela/recursao).
# Batch: precisa reprocessar a serie (ou um trecho longo) a cada observacao.
INCREMENTAL <- c(ma = "O(1)", ema = "O(1)", qes = "O(1)", kalman = "O(1)",
                 LOESS = "Batch", LOWESS = "Batch", spline = "Batch", winsor = "Batch",
                 fft = "Batch", wavelet = "Batch", emd = "Batch", remd = "Batch",
                 hp = "Batch", seasonal = "Batch", gabor = "Batch")

# Painel diagnostico
DIAG_WORST <- c("LOWESS", "remd", "gabor", "spline")
DIAG_CONSERVATIVE <- c("LOESS")   # o baseline entra sempre como referencia
DIAG_SCENARIO_OVERLAY <- "High_Noise"
DIAG_SERIES <- NULL               # NULL = primeira serie com ciclo de 24 h detectado
PARETO_SCENARIO <- "High_Noise"
PARETO_DELTA_CAP <- 100           # eixo y do Pareto (Delta RMSE%) e cortado aqui; acima disso o
                                  # ponto fica na borda e o rotulo mostra o valor real
SHOW_TITLES <- TRUE               # FALSE remove titulo/subtitulo de dentro das figuras (legenda no LaTeX)
CONTAM_SERIES <- NULL             # serie da figura 0 (NULL = a mesma do painel diagnostico)
CONTAM_LAST_N <- 200L
CONTAM_SCENARIOS <- c("Clean", "Mod_Noise", "High_Noise", "Sev_Noise")  # paineis da figura 0
TAB_DIR <- "tabelas"
TABLE_SIG <- "wilcoxon"           # negrito nas tabelas: "wilcoxon" (Wilcoxon-Holm) ou "nemenyi" (rank < baseline - CD)
HEAT_LIMIT <- 50                  # Delta RMSE% no extremo da escala de cor (valores acima saturam)

FILTER_ORDER <- c("non_filter", "ma", "ema", "LOESS", "LOWESS", "spline", "winsor", "qes",
                  "fft", "wavelet", "emd", "remd", "hp", "kalman", "seasonal", "gabor")
FILTER_LABEL <- c(non_filter = "Baseline", ma = "MA", ema = "EMA", LOESS = "Prog. Smooth",
                  LOWESS = "LOWESS", spline = "Spline", winsor = "Winsor", qes = "QES",
                  fft = "FFT", wavelet = "Wavelet", emd = "EMD", remd = "REMD", hp = "HP",
                  kalman = "Kalman", seasonal = "Seasonal", gabor = "Gabor")
SCENARIO_ORDER <- c("Clean", "Mod_Noise", "High_Noise", "Sev_Noise",
                    "Gauss_Mod", "Gauss_High", "Gauss_Sev", "Out_Mod", "Out_High", "Out_Sev")
PLOT_SCENARIOS <- NULL            # NULL = todos os cenarios presentes; ou um subconjunto

# Correcao para multiplos testes (Wilcoxon filtro vs. baseline):
#   "family"      Holm entre os 15 filtros de cada modelo x cenario x horizonte (original)
#   "global_holm" Holm sobre TODOS os testes do artigo (todos os modelos, cenarios e horizontes)
#   "global_bh"   Benjamini-Hochberg (FDR) sobre todos os testes
# As tres colunas saem em figuras/tabela_delta_rmse_wilcoxon.csv; esta opcao escolhe a usada
# nos asteriscos do heatmap e no negrito das tabelas.
MULTIPLICITY <- "family"
DM_MULTIPLICITY <- "none"         # "bh": BH entre as series de cada filtro (barras do DM)

# Filtros causais medidos com verifica_causalidade.R (os demais revisam o passado).
CAUSAL_FILTERS <- c("ma", "ema", "kalman")
SPLIT_CAUSALITY <- TRUE           # separa causais e nao causais nas figuras e tabelas

# Textos das figuras (ingles)
TXT <- list(
  rank = "Mean rank (Friedman)",
  delta = "Median \u0394RMSE% vs. baseline",
  delta_short = "Median \u0394RMSE%",
  heat_caption = "* Paired Wilcoxon signed-rank test vs. baseline, %s p \u2264 0.05.\nBlue = error reduction; red = degradation. Causal = output at t uses only data up to t.",
  mult_family = "Holm-adjusted (within each model, scenario and horizon)",
  mult_global_holm = "Holm-adjusted across all tests",
  mult_global_bh = "Benjamini-Hochberg-adjusted (FDR) across all tests",
  causal = "Causal", noncausal = "Non-causal",
  cd_group = "%s  |  h = %d  |  %s  |  %s",
  cd_sub = "N = %d series, k = %d methods; bars connect methods that are not significantly different (Nemenyi, \u03b1 = %.2f)",
  overhead = "Filter overhead (ms, log scale)",
  incremental = "Update mode",
  frontier = "Pareto frontier",
  pareto_title = "Accuracy vs. computational cost",
  pareto_sub = "%s scenario, averaged over h = %s. Dashed line = baseline (no filter).\nSolid line = %s (lower cost and lower error).",
  pareto_cap = "Points above %s%% are drawn at the top edge; the label gives the actual value.",
  step = "Step within the forecast window",
  ratio = "Error ratio (filter RMSE / baseline RMSE)",
  er_title = "(a) %s: per-step error ratio (h = %d)",
  er_sub = "Median across series of the filter RMSE / baseline RMSE at each step; > 1 = worse than the baseline (log scale).",
  ov_title = "(b) Last 48 h of the training window: series %s, %s scenario",
  ov_sub = "Coloured line = filtered series. '24 h cycle retained' = amplitude of the 24 h sinusoid in the filtered series / in the clean series.",
  ov_strip = "%s (%s)\n24 h cycle retained: %.0f%%",
  hours = "Hour (last 48 h of the training window)",
  value = "Value",
  dm_title = "Series-level Diebold-Mariano test (h = %d)",
  dm_sub = "Share of the %d series in each class; filters ordered by mean share of significant improvement.",
  dm_better = "Significantly better",
  dm_equal = "Not significantly different",
  dm_worse = "Significantly worse",
  dm_share = "Share of series",
  dm_caption = "Diebold-Mariano test (squared-error loss) vs. baseline; one-sided p < 0.05 in each direction%s.",
  clean = "Clean series", corrupted = "Corrupted series", filtered = "Filtered",
  contam_title = "Contamination profile: series %s",
  contam_sub = "Last %d observations; dual-component contamination applied across the full series prior to windowing.",
  contam_x = "Hourly observation", contam_y = "Value"
)

# Cores (paleta validada para daltonismo; texto sempre em tinta neutra)
COL <- list(ink = "#0b0b0b", ink2 = "#52514e", grid = "#e4e3de", surface = "#fcfcfb",
            neg = "#1c5cab", mid = "#f0efec", pos = "#b8322f",
            good = "#0ca30c", neutral = "#c9c8c2", bad = "#d03b3b",
            o1 = "#2a78d6", batch = "#eb6834",
            cat = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"))

theme_article <- function(base = 9) {
  theme_minimal(base_size = base) %+replace%
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = COL$grid, linewidth = 0.3),
          axis.text = element_text(colour = COL$ink2),
          axis.title = element_text(colour = COL$ink2),
          plot.title = element_text(colour = COL$ink, face = "bold", hjust = 0, size = base + 2,
                                    margin = margin(b = 4)),
          plot.subtitle = element_text(colour = COL$ink2, hjust = 0, margin = margin(b = 6)),
          plot.caption = element_text(colour = COL$ink2, hjust = 0, size = base - 1),
          strip.text = element_text(colour = COL$ink, face = "bold", margin = margin(3, 3, 3, 3)),
          legend.position = "top", legend.title = element_text(colour = COL$ink2),
          legend.text = element_text(colour = COL$ink2),
          plot.background = element_rect(fill = "white", colour = NA))
}

save_fig <- function(p, name, width, height) {
  dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = width, height = height, dpi = 300, bg = "white")
  pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = width, height = height, bg = "white", device = pdf_dev)
  message("  figura: ", file.path(FIG_DIR, paste0(name, ".png")))
}

titles <- function(title = NULL, subtitle = NULL) {
  if (SHOW_TITLES) labs(title = title, subtitle = subtitle) else labs(title = NULL, subtitle = NULL)
}

scen_levels <- function(x) { x <- unique(as.character(x)); c(SCENARIO_ORDER[SCENARIO_ORDER %in% x], sort(setdiff(x, SCENARIO_ORDER))) }
causal_group <- function(f) factor(ifelse(f %in% c("non_filter", CAUSAL_FILTERS), TXT$causal, TXT$noncausal),
                                   levels = c(TXT$causal, TXT$noncausal))
mult_label <- function() switch(MULTIPLICITY, family = TXT$mult_family, global_holm = TXT$mult_global_holm,
                                global_bh = TXT$mult_global_bh, stop("MULTIPLICITY invalido"))

flab <- function(f) unname(ifelse(f %in% names(FILTER_LABEL), FILTER_LABEL[f], f))
scen_lab <- function(s) gsub("_", " ", s)

# ------------------------------ Leitura dos dados ----------------------------
load_results <- function() {
  # Nome esperado: <MODELO>_industrial_sensors_results_<N>pts.xlsx (aceita prefixo, ex.: "c07115eb-ARIMA_...").
  pat <- sprintf("^(?:.*[-_ ])?(%s)_industrial_sensors_results_([0-9]+)pts\\.xlsx$",
                 paste(MODELS[order(-nchar(MODELS))], collapse = "|"))
  files <- list.files(RESULTS_DIR, pattern = "\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
  files <- files[!startsWith(basename(files), "~$") & grepl(pat, basename(files), perl = TRUE)]
  if (!length(files)) stop("Nenhum arquivo <MODELO>_industrial_sensors_results_<N>pts.xlsx em ", normalizePath(RESULTS_DIR))
  info <- t(vapply(basename(files), function(b) {
    m <- regmatches(b, regexec(pat, b, perl = TRUE))[[1]]; c(m[2], m[3])
  }, character(2)))
  key <- paste(info[, 1], info[, 2])
  if (anyDuplicated(key)) {
    # Dois arquivos para o mesmo modelo e horizonte: fica o mais recente.
    mt <- file.mtime(files); o <- order(mt, decreasing = TRUE)
    keep <- o[!duplicated(key[o])]
    message("  Aviso: arquivos repetidos para o mesmo modelo/horizonte; usando o mais recente: ",
            paste(basename(files[setdiff(seq_along(files), keep)]), collapse = ", "), " ignorado(s)")
    files <- files[sort(keep)]; info <- info[sort(keep), , drop = FALSE]
  }
  res <- list(); err <- list()
  for (k in seq_along(files)) {
    f <- files[k]
    model <- info[k, 1]; h <- as.integer(info[k, 2]) %/% N_WINDOWS_RESULTS
    r <- as.data.frame(readxl::read_excel(f, sheet = "Resultados_Gerais"))
    r$Model <- model; r$Horizon <- h
    if (!"Replication" %in% names(r)) r$Replication <- 1L
    res[[f]] <- r
    if ("Erros_Teste" %in% readxl::excel_sheets(f)) {
      e <- as.data.frame(readxl::read_excel(f, sheet = "Erros_Teste"))
      if (nrow(e)) {
        names(e)[names(e) == "Horizon"] <- "Step"
        e$Model <- model; e$Horizon <- h
        if (!"Replication" %in% names(e)) e$Replication <- 1L
        err[[f]] <- e[, c("Model", "Horizon", "Series_ID", "Scenario", "Replication", "Filter", "Window", "Step", "Error")]
      }
    } else {
      message("  ", basename(f), ": sem aba Erros_Teste (Error Ratio por passo indisponivel para esse arquivo)")
    }
    message("  lido: ", basename(f), " (", model, ", h = ", h, ", ", nrow(r), " linhas)")
  }
  cols <- c("Model", "Horizon", "Series_ID", "Scenario", "Replication", "Filter", "Filter_Param",
            "Status", "RMSE", "DM_Stat", "DM_Pvalue")
  raw <- do.call(rbind, lapply(res, function(r) r[, cols]))
  raw <- raw[raw$Status == "CONVERGED" & is.finite(raw$RMSE), ]
  err <- if (length(err)) do.call(rbind, err) else NULL
  if (!is.null(PLOT_SCENARIOS)) {
    raw <- raw[raw$Scenario %in% PLOT_SCENARIOS, ]
    if (!is.null(err)) err <- err[err$Scenario %in% PLOT_SCENARIOS, ]
  }
  # Com varias realizacoes da contaminacao, as analises por serie (ranks, Wilcoxon,
  # heatmap, Pareto) usam o RMSE medio das realizacoes; o DM usa cada realizacao.
  res <- raw
  if (max(raw$Replication) > 1L) {
    key <- c("Model", "Horizon", "Series_ID", "Scenario", "Filter")
    res <- stats::aggregate(RMSE ~ Model + Horizon + Series_ID + Scenario + Filter, raw, mean)
    k1 <- do.call(paste, raw[key]); k2 <- do.call(paste, res[key])
    md <- tapply(raw$Filter_Param, k1, function(v) names(sort(table(v), decreasing = TRUE))[1])
    res$Filter_Param <- unname(md[k2])
    res$Status <- "CONVERGED"
    message("  ", max(raw$Replication), " realizacoes da contaminacao: ranks e Wilcoxon usam o RMSE medio por serie")
  }
  list(res = res, raw = raw, err = err)
}

# Matriz series x filtros (so series completas, como no IC_Common.R)
rmse_matrix <- function(res, model, h, scen, filters = FILTER_ORDER) {
  d <- res[res$Model == model & res$Horizon == h & res$Scenario == scen, c("Series_ID", "Filter", "RMSE")]
  if (!nrow(d)) return(NULL)
  w <- reshape(d, idvar = "Series_ID", timevar = "Filter", direction = "wide")
  names(w) <- sub("^RMSE\\.", "", names(w))
  fs <- intersect(filters, names(w))
  m <- as.matrix(w[, fs, drop = FALSE]); rownames(m) <- w$Series_ID
  m[stats::complete.cases(m), , drop = FALSE]
}

combos <- function(res) {
  u <- unique(res[, c("Model", "Horizon", "Scenario")])
  u <- u[order(match(u$Model, MODELS), u$Horizon, match(u$Scenario, scen_levels(u$Scenario))), ]
  rownames(u) <- NULL; u
}

# ------------------------- 1. Diferenca Critica ------------------------------
nemenyi_cd <- function(k, n, alpha = ALPHA) stats::qtukey(1 - alpha, k, Inf) / sqrt(2) * sqrt(k * (k + 1) / (6 * n))

cd_cliques <- function(r, cd) {
  r <- sort(r); k <- length(r); out <- list(); last_j <- 0L
  for (i in seq_len(k)) {
    j <- max(which(r - r[i] < cd))
    if (j > i && j > last_j) { out[[length(out) + 1L]] <- c(i, j); last_j <- j }
  }
  out
}

plot_cd <- function(m, title) {
  k <- ncol(m); n <- nrow(m)
  ranks <- colMeans(t(apply(m, 1, rank, ties.method = "average")))
  cd <- nemenyi_cd(k, n)
  r <- sort(ranks); nm <- names(r)
  nl <- ceiling(k / 2); step <- 0.55
  left <- seq_len(nl); right <- (nl + 1L):k
  cl <- cd_cliques(ranks, cd)
  bar_gap <- 0.16
  top <- 0.45 + bar_gap * max(1L, length(cl))
  lab_y <- c(-(seq_along(left)) * step - top, -(rev(seq_along(right))) * step - top)
  lab_x <- c(rep(0.6, length(left)), rep(k + 0.4, length(right)))
  is_left <- seq_len(k) <= nl
  el <- data.frame(x = r, y0 = 0, y1 = lab_y, xe = lab_x, hj = ifelse(is_left, 1, 0),
                   base = nm == "non_filter")
  el$lab <- ifelse(is_left, sprintf("%s  %.2f", flab(nm), r), sprintf("%.2f  %s", r, flab(nm)))
  bars <- if (length(cl)) data.frame(x0 = r[vapply(cl, `[`, 1, 1)] - 0.06, x1 = r[vapply(cl, `[`, 1, 2)] + 0.06,
                                      y = -0.25 - (seq_along(cl) - 1) * bar_gap) else NULL
  ticks <- data.frame(x = seq_len(k))
  p <- ggplot() +
    geom_segment(aes(x = 1, xend = k, y = 0, yend = 0), colour = COL$ink, linewidth = 0.4) +
    geom_segment(data = ticks, aes(x = x, xend = x, y = 0, yend = 0.12), colour = COL$ink, linewidth = 0.3) +
    geom_text(data = ticks, aes(x = x, y = 0.3, label = x), colour = COL$ink2, size = 2.6) +
    geom_segment(aes(x = 1, xend = 1 + cd, y = 0.75, yend = 0.75), colour = COL$ink, linewidth = 0.6) +
    geom_segment(aes(x = c(1, 1 + cd), xend = c(1, 1 + cd), y = 0.68, yend = 0.82), colour = COL$ink, linewidth = 0.4) +
    annotate("text", x = 1 + cd / 2, y = 0.98, label = sprintf("CD = %.2f", cd), colour = COL$ink, size = 2.7) +
    geom_segment(data = el, aes(x = x, xend = x, y = y0, yend = y1), colour = COL$ink2, linewidth = 0.3) +
    geom_segment(data = el, aes(x = x, xend = xe, y = y1, yend = y1), colour = COL$ink2, linewidth = 0.3) +
    geom_text(data = el, aes(x = xe + ifelse(hj == 1, -0.1, 0.1), y = y1, label = lab, hjust = hj,
                             fontface = ifelse(base, "bold", "plain")), colour = COL$ink, size = 2.7)
  if (!is.null(bars))
    p <- p + geom_segment(data = bars, aes(x = x0, xend = x1, y = y, yend = y), colour = COL$ink,
                          linewidth = 0.9, lineend = "round")
  p + scale_x_continuous(limits = c(-3.2, k + 4.2)) +
    scale_y_continuous(limits = c(min(lab_y) - 0.4, 1.15)) +
    labs(title = title, subtitle = if (SHOW_TITLES) sprintf(TXT$cd_sub, n, k, ALPHA) else NULL) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold", colour = COL$ink, size = 10),
          plot.subtitle = element_text(colour = COL$ink2, size = 7.5, margin = margin(b = 2)),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(6, 6, 6, 6))
}

fig_cd <- function(res) {
  message("1. Diagramas de Diferenca Critica")
  cb <- combos(res)
  for (key in unique(paste(cb$Model, cb$Horizon))) {
    sub <- cb[paste(cb$Model, cb$Horizon) == key, ]
    plots <- lapply(seq_len(nrow(sub)), function(i) {
      m <- rmse_matrix(res, sub$Model[i], sub$Horizon[i], sub$Scenario[i])
      if (is.null(m) || nrow(m) < 2L || ncol(m) < 3L) return(NULL)
      plot_cd(m, sprintf("%s  |  h = %d  |  %s", sub$Model[i], sub$Horizon[i], scen_lab(sub$Scenario[i])))
    })
    plots <- Filter(Negate(is.null), plots)
    if (!length(plots)) next
    p <- wrap_plots(plots, ncol = 1)
    save_fig(p, sprintf("1_cd_%s_h%02d", sub$Model[1], sub$Horizon[1]), 7.2, 3.2 * length(plots))
    if (SPLIT_CAUSALITY) {
      # Rankings separados (cada um com o baseline): k e CD mudam com o tamanho do grupo.
      groups <- list(causal = c("non_filter", CAUSAL_FILTERS),
                     noncausal = c("non_filter", setdiff(FILTER_ORDER, c("non_filter", CAUSAL_FILTERS))))
      for (g in names(groups)) {
        gl <- if (g == "causal") TXT$causal else TXT$noncausal
        pg <- lapply(seq_len(nrow(sub)), function(i) {
          m <- rmse_matrix(res, sub$Model[i], sub$Horizon[i], sub$Scenario[i], filters = groups[[g]])
          if (is.null(m) || nrow(m) < 2L || ncol(m) < 3L) return(NULL)
          plot_cd(m, sprintf(TXT$cd_group, sub$Model[i], sub$Horizon[i], scen_lab(sub$Scenario[i]), gl))
        })
        pg <- Filter(Negate(is.null), pg)
        if (!length(pg)) next
        hh <- if (g == "causal") 2.0 else 3.0
        save_fig(wrap_plots(pg, ncol = 1), sprintf("1_cd_%s_h%02d_%s", sub$Model[1], sub$Horizon[1], g), 7.2, hh * length(pg))
      }
    }
  }
}

# ---------------------- 2. Heatmap de Delta RMSE% ----------------------------
delta_table <- function(res) {
  out <- list()
  cb <- combos(res)
  for (i in seq_len(nrow(cb))) {
    m <- rmse_matrix(res, cb$Model[i], cb$Horizon[i], cb$Scenario[i])
    if (is.null(m) || !"non_filter" %in% colnames(m) || nrow(m) < 2L) next
    fs <- setdiff(colnames(m), "non_filter")
    base <- m[, "non_filter"]
    raw_p <- vapply(fs, function(f) tryCatch(stats::wilcox.test(m[, f], base, paired = TRUE, exact = FALSE)$p.value,
                                             error = function(e) NA_real_), numeric(1))
    out[[i]] <- data.frame(Model = cb$Model[i], Horizon = cb$Horizon[i], Scenario = cb$Scenario[i],
                           Filter = fs,
                           Delta = vapply(fs, function(f) stats::median((m[, f] / base - 1) * 100), numeric(1)),
                           P_Raw = raw_p, P_Holm_Family = stats::p.adjust(raw_p, "holm"), N = nrow(m))
  }
  d <- do.call(rbind, out); rownames(d) <- NULL
  d$P_Holm_Global <- stats::p.adjust(d$P_Raw, "holm")
  d$Q_BH_Global <- stats::p.adjust(d$P_Raw, "BH")
  d$P_Holm <- switch(MULTIPLICITY, family = d$P_Holm_Family, global_holm = d$P_Holm_Global,
                     global_bh = d$Q_BH_Global)          # coluna usada nas figuras/tabelas
  d$Causality <- ifelse(d$Filter %in% CAUSAL_FILTERS, "causal", "non_causal")
  d
}

# Quantos filtros sao significativos (melhor ou pior) em cada correcao: apoio ao texto.
multiplicity_summary <- function(dt) {
  f <- function(p) is.finite(p) & p <= ALPHA
  agg <- function(col) stats::aggregate(cbind(Better = f(dt[[col]]) & dt$Delta < 0, Worse = f(dt[[col]]) & dt$Delta > 0) ~
                                          Model + Horizon + Scenario, dt, sum)
  a <- agg("P_Holm_Family"); b <- agg("P_Holm_Global"); c <- agg("Q_BH_Global")
  names(a)[4:5] <- c("Better_HolmFamily", "Worse_HolmFamily")
  names(b)[4:5] <- c("Better_HolmGlobal", "Worse_HolmGlobal")
  names(c)[4:5] <- c("Better_BHGlobal", "Worse_BHGlobal")
  out <- Reduce(function(x, y) merge(x, y, by = c("Model", "Horizon", "Scenario")), list(a, b, c))
  out$Tests_Total <- nrow(dt)
  out[order(match(out$Model, MODELS), out$Horizon, match(out$Scenario, scen_levels(out$Scenario))), ]
}

fig_heatmap <- function(dt) {
  message("2. Heatmap de Delta RMSE%")
  d <- dt
  d$X <- factor(paste0("h=", d$Horizon), levels = paste0("h=", sort(unique(d$Horizon))))
  d$Scenario <- factor(scen_lab(d$Scenario), levels = scen_lab(scen_levels(d$Scenario)))
  d$Y <- factor(flab(d$Filter), levels = rev(flab(setdiff(FILTER_ORDER, "non_filter"))))
  d$Model <- factor(d$Model, levels = intersect(MODELS, d$Model))
  d$Group <- causal_group(d$Filter)
  lim <- HEAT_LIMIT
  d$sig <- ifelse(is.finite(d$P_Holm) & d$P_Holm <= ALPHA, "*", "")
  d$txtcol <- ifelse(abs(pmin(pmax(d$Delta, -lim), lim)) / lim > 0.55, "white", COL$ink)
  p <- ggplot(d, aes(X, Y, fill = Delta)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sig, colour = txtcol), size = 4, vjust = 0.75, fontface = "bold") +
    scale_colour_identity() +
    scale_fill_gradient2(low = COL$neg, mid = COL$mid, high = COL$pos, midpoint = 0,
                         limits = c(-lim, lim), oob = scales::squish,
                         breaks = c(-lim, -lim / 2, 0, lim / 2, lim),
                         labels = function(x) ifelse(abs(x) == lim, paste0(ifelse(x > 0, "\u2265 +", "\u2264 "), x, "%"),
                                                     paste0(ifelse(x > 0, "+", ""), x, "%")),
                         name = TXT$delta_short,
                         guide = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
    (if (SPLIT_CAUSALITY) facet_grid(Model + Group ~ Scenario, scales = "free", space = "free", switch = "x")
     else facet_grid(Model ~ Scenario, scales = "free_x", space = "free_x", switch = "x")) +
    titles(TXT$delta) + labs(x = NULL, y = NULL, caption = sprintf(TXT$heat_caption, mult_label())) +
    theme_article() +
    theme(panel.grid = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_text(size = 7.5),
          strip.placement = "outside", panel.spacing.x = unit(0.35, "lines"),
          legend.key.width = unit(1.6, "cm"), legend.key.height = unit(0.3, "cm"))
  save_fig(p, "2_heatmap_delta_rmse", max(7.2, 1.6 + 0.62 * nlevels(d$Scenario) * nlevels(d$X)),
           1.2 + 3.3 * nlevels(d$Model))
}

# -------------------------- 3. Fronteira de Pareto ---------------------------
mode_param <- function(res, f) {
  v <- res$Filter_Param[res$Filter == f]
  if (!length(v)) return(NA_character_)
  names(sort(table(v), decreasing = TRUE))[1]
}

ic_common_loaded <- function() exists("get_filter_candidates", mode = "function") && exists("apply_filter_candidate", mode = "function")

load_ic_common <- function() {
  if (ic_common_loaded()) return(TRUE)
  if (!file.exists("IC_Common.R")) { message("  IC_Common.R nao encontrado nesta pasta."); return(FALSE) }
  ok <- tryCatch({ suppressMessages(source("IC_Common.R")); TRUE },
                 error = function(e) { message("  Falha ao carregar IC_Common.R: ", conditionMessage(e)); FALSE })
  ok && ic_common_loaded()
}

candidate_for <- function(f, desc, freq = 24L) {
  j <- match(f, filter_names)
  cands <- get_filter_candidates(j, freq)
  hit <- Filter(function(c) identical(c$desc, desc), cands)
  if (length(hit)) hit[[1]] else cands[[1]]
}

measure_overhead <- function(res) {
  if (!load_ic_common()) return(NULL)
  message("  medindo o custo de cada filtro em ", OVERHEAD_BENCH_SERIES, " series (uma vez; salvo em ",
          OVERHEAD_MEASURED_CSV, ")")
  sl <- load_m4_hourly(OVERHEAD_BENCH_SERIES)
  n_test <- N_WINDOWS * PRED_HORIZON
  sc <- noise_scenarios[[PARETO_SCENARIO]]
  xs <- lapply(seq_along(sl), function(i) {
    cl <- as.numeric(sl[[i]])
    corrupt_sensor_signal(cl, sc$noise_lvl, sc$outlier_ratio, sc$spike_mag, BASE_SEED + i)[seq_len(length(cl) - n_test)]
  })
  fs <- setdiff(FILTER_ORDER, "non_filter")
  ms <- vapply(fs, function(f) {
    cand <- candidate_for(f, mode_param(res, f))
    t <- vapply(xs, function(x) {
      reps <- 1L
      repeat {
        el <- system.time(for (r in seq_len(reps)) apply_filter_candidate(cand, x))[["elapsed"]]
        if (el >= 0.05 || reps >= 512L) break
        reps <- reps * 4L
      }
      el / reps
    }, numeric(1))
    stats::median(t) * 1000
  }, numeric(1))
  out <- data.frame(Filter = fs, Overhead_ms = pmax(ms, 0.01), Param = vapply(fs, function(f) mode_param(res, f), ""))
  utils::write.csv(out, OVERHEAD_MEASURED_CSV, row.names = FALSE)
  out
}

get_overhead <- function(res) {
  if (file.exists(OVERHEAD_CSV)) { message("  overhead lido de ", OVERHEAD_CSV); return(utils::read.csv(OVERHEAD_CSV)) }
  if (file.exists(OVERHEAD_MEASURED_CSV)) { message("  overhead lido de ", OVERHEAD_MEASURED_CSV); return(utils::read.csv(OVERHEAD_MEASURED_CSV)) }
  measure_overhead(res)
}

pareto_front <- function(x, y) {
  o <- order(x, y); best <- Inf; keep <- logical(length(x))
  for (i in o) if (y[i] < best) { keep[i] <- TRUE; best <- y[i] }
  keep
}

fig_pareto <- function(res, dt) {
  message("3. Fronteira de Pareto")
  oh <- get_overhead(res)
  if (is.null(oh)) { message("  sem overhead: figura 3 pulada"); return(invisible()) }
  rk <- list()
  cb <- combos(res); cb <- cb[cb$Scenario == PARETO_SCENARIO, ]
  for (i in seq_len(nrow(cb))) {
    m <- rmse_matrix(res, cb$Model[i], cb$Horizon[i], cb$Scenario[i])
    if (is.null(m) || nrow(m) < 2L) next
    r <- colMeans(t(apply(m, 1, rank, ties.method = "average")))
    rk[[i]] <- data.frame(Model = cb$Model[i], Horizon = cb$Horizon[i], Filter = names(r), Rank = r)
  }
  rk <- do.call(rbind, rk)
  rk <- aggregate(Rank ~ Model + Filter, rk, mean)
  dl <- aggregate(Delta ~ Model + Filter, dt[dt$Scenario == PARETO_SCENARIO, ], stats::median)
  base_rank <- rk[rk$Filter == "non_filter", c("Model", "Rank")]
  d <- merge(rk[rk$Filter != "non_filter", ], dl, by = c("Model", "Filter"), all.x = TRUE)
  d <- merge(d, oh[, c("Filter", "Overhead_ms")], by = "Filter")
  d$Update <- factor(ifelse(d$Filter %in% names(INCREMENTAL), INCREMENTAL[d$Filter], "Batch"), levels = c("O(1)", "Batch"))
  d$Model <- factor(d$Model, levels = intersect(MODELS, d$Model))
  d$Label <- flab(d$Filter)
  hs <- paste(sort(unique(cb$Horizon)), collapse = ", ")
  one <- function(ycol, ylab, ref, file, cap = Inf) {
    d$Y <- d[[ycol]]
    ref$Model <- factor(ref$Model, levels = levels(d$Model))
    d <- d[is.finite(d$Y), ]
    d$front <- FALSE
    for (mm in levels(d$Model)) { ii <- d$Model == mm; d$front[ii] <- pareto_front(d$Overhead_ms[ii], d$Y[ii]) }
    # Valores extremos (ex.: LOWESS +800%) achatariam o resto: corta o eixo e mostra o valor real no rotulo.
    d$clipped <- d$Y > cap
    d$Yp <- pmin(d$Y, cap)
    d$Label <- ifelse(d$clipped, sprintf("%s (%+.0f%%)", d$Label, d$Y), d$Label)
    fr <- d[d$front, ]; fr <- fr[order(fr$Model, fr$Overhead_ms), ]
    p <- ggplot(d, aes(Overhead_ms, Yp)) +
      geom_hline(data = ref, aes(yintercept = Y), linetype = "dashed", colour = COL$ink2, linewidth = 0.4) +
      geom_step(data = fr, aes(group = Model), direction = "hv", colour = COL$ink, linewidth = 0.5, alpha = 0.6) +
      geom_point(aes(fill = Update, shape = Update), size = 2.6, colour = "white", stroke = 0.6) +
      geom_text_repel(aes(label = Label, fontface = ifelse(front, "bold", "plain")), size = 2.5,
                      colour = COL$ink, min.segment.length = 0, segment.colour = COL$ink2,
                      segment.size = 0.25, box.padding = 0.45, point.padding = 0.25, force = 3,
                      force_pull = 0.5, max.overlaps = Inf, max.iter = 20000, seed = 1) +
      scale_x_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE, drop0trailing = TRUE, trim = TRUE),
                    expand = expansion(mult = 0.08)) +
      scale_y_continuous(expand = expansion(mult = 0.12)) +
      scale_fill_manual(values = c("O(1)" = COL$o1, "Batch" = COL$batch), name = TXT$incremental, drop = FALSE) +
      scale_shape_manual(values = c("O(1)" = 21, "Batch" = 24), name = TXT$incremental, drop = FALSE) +
      facet_wrap(~ Model, ncol = 1, scales = "free_y") +
      titles(TXT$pareto_title, sprintf(TXT$pareto_sub, scen_lab(PARETO_SCENARIO), hs, TXT$frontier)) +
      labs(x = TXT$overhead, y = ylab, caption = if (any(d$clipped)) sprintf(TXT$pareto_cap, cap) else NULL) +
      theme_article()
    save_fig(p, file, 7.2, 1.3 + 3.4 * nlevels(d$Model))
  }
  one("Rank", TXT$rank, setNames(base_rank, c("Model", "Y")), "3_pareto_rank")
  one("Delta", TXT$delta_short, data.frame(Model = levels(d$Model), Y = 0), "3_pareto_delta_rmse", cap = PARETO_DELTA_CAP)
}

# -------------------- 4. Painel diagnostico de degradacao --------------------
error_ratio_table <- function(err, model) {
  e <- err[err$Model == model, ]
  if (!nrow(e)) return(NULL)
  h <- max(e$Horizon); e <- e[e$Horizon == h, ]
  fs <- c(DIAG_WORST, DIAG_CONSERVATIVE)
  r <- aggregate(Error ~ Series_ID + Scenario + Filter + Step, e[e$Filter %in% c(fs, "non_filter"), ],
                 function(v) sqrt(mean(v^2)))
  b <- r[r$Filter == "non_filter", ]; names(b)[names(b) == "Error"] <- "Base"
  r <- merge(r[r$Filter != "non_filter", ], b[, c("Series_ID", "Scenario", "Step", "Base")])
  r$Ratio <- r$Error / r$Base
  r <- r[is.finite(r$Ratio) & r$Ratio > 0, ]
  out <- aggregate(Ratio ~ Scenario + Filter + Step, r, stats::median)
  out$Horizon <- h; out
}

# Serie usada nas figuras 0 e 4b: a indicada ou a primeira com ciclo de 24 h detectado.
pick_series <- function(sl, wanted = DIAG_SERIES) {
  if (!is.null(wanted) && wanted %in% names(sl)) return(wanted)
  n_test <- N_WINDOWS * PRED_HORIZON
  for (id in names(sl)) {
    cl <- as.numeric(sl[[id]]); cl <- cl[seq_len(length(cl) - n_test)]
    if (forecast::findfrequency(cl) == 24) return(id)
  }
  names(sl)[1]
}

overlay_data <- function(res, model) {
  if (!load_ic_common()) return(NULL)
  sl <- load_m4_hourly(N_SERIES)
  n_test <- N_WINDOWS * PRED_HORIZON
  sc <- noise_scenarios[[DIAG_SCENARIO_OVERLAY]]
  pick <- pick_series(sl)
  i <- match(pick, names(sl)); cl <- as.numeric(sl[[i]])
  co <- corrupt_sensor_signal(cl, sc$noise_lvl, sc$outlier_ratio, sc$spike_mag, BASE_SEED + i)
  n <- length(cl) - n_test; x <- co[seq_len(n)]; y <- cl[seq_len(n)]
  amp24 <- function(v) { t <- seq_along(v); unname(sqrt(sum(stats::coef(stats::lm(v ~ sin(2 * pi * t / 24) + cos(2 * pi * t / 24)))[2:3]^2))) }
  a_clean <- amp24(y)
  rs <- res[res$Model == model & res$Scenario == DIAG_SCENARIO_OVERLAY, ]
  fs <- c(DIAG_WORST, DIAG_CONSERVATIVE)
  keep <- (n - 47L):n
  rows <- list()
  for (f in fs) {
    desc <- mode_param(rs, f)
    out <- apply_filter_candidate(candidate_for(f, desc), x)
    if (is.null(out$series)) next
    rows[[f]] <- data.frame(Filter = f, Param = desc, Hour = seq_along(keep),
                            Clean = y[keep], Corrupted = x[keep], Filtered = out$series[keep],
                            Amp = 100 * amp24(out$series) / a_clean)
  }
  d <- do.call(rbind, rows); d$Series_ID <- pick; d
}

# Linha filtrada com a cor de cada filtro, sem entrar na legenda de cores
ggnewscale_free_filtered <- function(ov) {
  geom_line(data = ov, aes(x = Hour, y = Filtered), colour = ov$col, linewidth = 0.8, inherit.aes = FALSE)
}

fig_diagnostic <- function(res, err) {
  message("4. Painel diagnostico de degradacao")
  fs <- c(DIAG_WORST, DIAG_CONSERVATIVE)
  pal <- setNames(COL$cat[seq_along(fs)], flab(fs))
  for (model in intersect(MODELS, unique(res$Model))) {
    pa <- NULL
    er <- if (!is.null(err)) error_ratio_table(err, model) else NULL
    if (!is.null(er) && nrow(er)) {
      er$Scenario <- factor(scen_lab(er$Scenario), levels = scen_lab(SCENARIO_ORDER))
      er$F <- factor(flab(er$Filter), levels = flab(fs))
      pa <- ggplot(er, aes(Step, Ratio, colour = F)) +
        geom_hline(yintercept = 1, linetype = "dashed", colour = COL$ink2, linewidth = 0.4) +
        geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
        scale_y_log10() + scale_x_continuous(breaks = seq_len(max(er$Step))) +
        scale_colour_manual(values = pal, name = NULL) +
        facet_wrap(~ Scenario, nrow = 1) +
        (if (SHOW_TITLES) labs(title = sprintf(TXT$er_title, model, er$Horizon[1]), subtitle = TXT$er_sub) else labs(title = "(a)")) +
        labs(             x = TXT$step, y = TXT$ratio) +
        theme_article()
    } else {
      message("  ", model, ": sem erros por passo (aba Erros_Teste); painel (a) pulado")
    }
    ov <- overlay_data(res, model)
    pb <- NULL
    if (!is.null(ov) && nrow(ov)) {
      strip <- sprintf(TXT$ov_strip, flab(ov$Filter), ov$Param, ov$Amp)
      ov$F <- factor(strip, levels = unique(strip[order(match(ov$Filter, fs))]))
      ov$col <- pal[flab(ov$Filter)]
      ref <- rbind(data.frame(Hour = ov$Hour, F = ov$F, Y = ov$Corrupted, Ref = TXT$corrupted),
                   data.frame(Hour = ov$Hour, F = ov$F, Y = ov$Clean, Ref = TXT$clean))
      ref$Ref <- factor(ref$Ref, levels = c(TXT$clean, TXT$corrupted))
      pb <- ggplot(ov, aes(Hour)) +
        geom_line(data = ref, aes(y = Y, colour = Ref, linewidth = Ref)) +
        scale_colour_manual(values = setNames(c(COL$ink, "#b5b4ae"), c(TXT$clean, TXT$corrupted)), name = NULL) +
        scale_linewidth_manual(values = setNames(c(0.5, 0.35), c(TXT$clean, TXT$corrupted)), name = NULL) +
        ggnewscale_free_filtered(ov) +
        scale_x_continuous(breaks = c(1, 12, 24, 36, 48)) +
        facet_wrap(~ F, nrow = 1) +
        (if (SHOW_TITLES) labs(title = sprintf(TXT$ov_title, ov$Series_ID[1], scen_lab(DIAG_SCENARIO_OVERLAY)),
                               subtitle = TXT$ov_sub) else labs(title = "(b)")) +
        labs(             x = TXT$hours, y = TXT$value) +
        theme_article() + theme(strip.text = element_text(size = 7.5))
    } else {
      message("  painel (b) pulado (precisa do IC_Common.R e dos ts_fil_*.R nesta pasta)")
    }
    parts <- Filter(Negate(is.null), list(pa, pb))
    if (!length(parts)) next
    save_fig(wrap_plots(parts, ncol = 1), sprintf("4_diagnostico_%s", model), 10, 3.6 * length(parts))
  }
}

# ---------------------- 5. Barras empilhadas do DM ---------------------------
fig_dm <- function(res) {
  message("5. Distribuicao Diebold-Mariano serie a serie")
  d <- res[res$Filter != "non_filter", ]
  identical0 <- is.finite(d$DM_Stat) & d$DM_Stat == 0 & is.finite(d$DM_Pvalue) & d$DM_Pvalue == 1
  # DM_Pvalue e unilateral "greater" (filtro mais preciso que o baseline); a outra
  # cauda e 1 - p (mesma estatistica t).
  pb <- ifelse(!is.finite(d$DM_Pvalue) | identical0, 1, d$DM_Pvalue)   # filtro melhor
  pw <- ifelse(!is.finite(d$DM_Pvalue) | identical0, 1, 1 - d$DM_Pvalue) # filtro pior
  if (DM_MULTIPLICITY == "bh") {
    g <- paste(d$Model, d$Horizon, d$Scenario, d$Filter)
    for (k in unique(g)) {
      ii <- which(g == k); adj <- stats::p.adjust(c(pb[ii], pw[ii]), "BH")
      pb[ii] <- adj[seq_along(ii)]; pw[ii] <- adj[length(ii) + seq_along(ii)]
    }
  }
  d$Class <- ifelse(pb < ALPHA, TXT$dm_better, ifelse(pw < ALPHA, TXT$dm_worse, TXT$dm_equal))
  d$Class <- factor(d$Class, levels = c(TXT$dm_worse, TXT$dm_equal, TXT$dm_better))
  for (h in sort(unique(d$Horizon))) {
    x <- d[d$Horizon == h, ]
    tab <- as.data.frame(table(Model = x$Model, Scenario = x$Scenario, Filter = x$Filter, Class = x$Class))
    tot <- aggregate(Freq ~ Model + Scenario + Filter, tab, sum)
    tab <- merge(tab, setNames(tot, c("Model", "Scenario", "Filter", "N")))
    tab <- tab[tab$N > 0, ]
    tab$Share <- tab$Freq / tab$N
    ordf <- aggregate(Share ~ Filter, tab[tab$Class == TXT$dm_better, ], mean)
    ordf <- ordf$Filter[order(ordf$Share)]
    tab$F <- factor(flab(as.character(tab$Filter)), levels = flab(as.character(ordf)))
    tab$Group <- causal_group(as.character(tab$Filter))
    tab$Scenario <- factor(scen_lab(as.character(tab$Scenario)), levels = scen_lab(scen_levels(tab$Scenario)))
    tab$Model <- factor(tab$Model, levels = intersect(MODELS, unique(tab$Model)))
    p <- ggplot(tab, aes(Share, F, fill = Class)) +
      geom_col(width = 0.75, colour = "white", linewidth = 0.4) +
      scale_fill_manual(values = setNames(c(COL$bad, COL$neutral, COL$good), levels(d$Class)), name = NULL,
                        breaks = c(TXT$dm_better, TXT$dm_equal, TXT$dm_worse)) +
      scale_x_continuous(labels = scales::percent, breaks = c(0.25, 0.5, 0.75), expand = c(0, 0)) +
      (if (SPLIT_CAUSALITY) facet_grid(Model + Group ~ Scenario, scales = "free_y", space = "free_y")
       else facet_grid(Model ~ Scenario)) +
      titles(sprintf(TXT$dm_title, h), sprintf(TXT$dm_sub, max(tab$N))) +
      labs(x = TXT$dm_share, y = NULL,
           caption = sprintf(TXT$dm_caption, if (DM_MULTIPLICITY == "bh") ", Benjamini-Hochberg-adjusted across series" else "")) +
      theme_article() + theme(panel.grid.major.y = element_blank(), panel.spacing.x = unit(0.8, "lines"))
    save_fig(p, sprintf("5_dm_series_h%02d", h), max(9, 1.6 + 1.9 * nlevels(tab$Scenario)), 1.6 + 3.1 * nlevels(tab$Model))
  }
}

# ------------------------ 0. Perfil de contaminacao --------------------------
fig_contamination <- function() {
  message("0. Perfil de contaminacao")
  if (!load_ic_common()) { message("  figura 0 pulada (precisa do IC_Common.R)"); return(invisible()) }
  sl <- load_m4_hourly(N_SERIES)
  id <- pick_series(sl, if (is.null(CONTAM_SERIES)) DIAG_SERIES else CONTAM_SERIES)
  i <- match(id, names(sl)); cl <- as.numeric(sl[[i]]); n <- length(cl)
  keep <- max(1L, n - CONTAM_LAST_N + 1L):n
  allsc <- c(noise_scenarios, if (exists("factorial_scenarios")) factorial_scenarios)
  scs <- intersect(CONTAM_SCENARIOS, names(allsc))
  rows <- lapply(scs, function(sn) {
    sc <- allsc[[sn]]
    # mesma semente em todos os cenarios (common random numbers), contaminacao na serie inteira
    co <- corrupt_sensor_signal(cl, sc$noise_lvl, sc$outlier_ratio, sc$spike_mag, BASE_SEED + i)
    rbind(data.frame(Scenario = sn, Obs = keep, Value = co[keep], Kind = TXT$corrupted),
          data.frame(Scenario = sn, Obs = keep, Value = cl[keep], Kind = TXT$clean))
  })
  d <- do.call(rbind, rows)
  d$Scenario <- factor(scen_lab(d$Scenario), levels = scen_lab(scs))
  d$Kind <- factor(d$Kind, levels = c(TXT$corrupted, TXT$clean))
  p <- ggplot(d, aes(Obs, Value, colour = Kind, linewidth = Kind)) +
    geom_line() +
    scale_colour_manual(values = setNames(c(COL$o1, COL$ink), c(TXT$corrupted, TXT$clean)), name = NULL) +
    scale_linewidth_manual(values = setNames(c(0.35, 0.45), c(TXT$corrupted, TXT$clean)), name = NULL) +
    facet_wrap(~ Scenario, ncol = 1) +
    titles(sprintf(TXT$contam_title, id), sprintf(TXT$contam_sub, length(keep))) +
    labs(x = TXT$contam_x, y = TXT$contam_y) +
    theme_article()
  save_fig(p, "0_contamination_profile", 7.2, 1.5 + 1.5 * length(scs))
}

# ------------------- Tabelas LaTeX de ranking de Friedman --------------------
# Uma tabela por modelo: linhas = filtros, colunas = cenario x horizonte, valor = rank
# medio. Negrito = filtro significativamente melhor que o baseline (TABLE_SIG).
# Usa \usepackage{booktabs} e \usepackage{rotating} (sidewaystable) no preambulo.
write_rank_tables <- function(res, dt) {
  message("Tabelas LaTeX de ranking")
  dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
  cb <- combos(res); out_all <- list()
  for (model in intersect(MODELS, unique(cb$Model))) {
    cm <- cb[cb$Model == model, ]
    cm <- cm[order(match(cm$Scenario, scen_levels(cm$Scenario)), cm$Horizon), ]
    cols <- list(); ns <- integer(0)
    for (i in seq_len(nrow(cm))) {
      m <- rmse_matrix(res, model, cm$Horizon[i], cm$Scenario[i])
      if (is.null(m) || nrow(m) < 2L) next
      r <- colMeans(t(apply(m, 1, rank, ties.method = "average")))
      k <- ncol(m); n <- nrow(m); cd <- nemenyi_cd(k, n)
      sig <- setNames(rep(FALSE, length(r)), names(r))
      if (TABLE_SIG == "nemenyi") {
        sig <- r < r[["non_filter"]] - cd
      } else {
        w <- dt[dt$Model == model & dt$Horizon == cm$Horizon[i] & dt$Scenario == cm$Scenario[i], ]
        ok <- w$Filter[is.finite(w$P_Holm) & w$P_Holm <= ALPHA & w$Delta < 0]
        sig[intersect(names(sig), ok)] <- TRUE
      }
      key <- paste(cm$Scenario[i], cm$Horizon[i])
      cols[[key]] <- list(scen = cm$Scenario[i], h = cm$Horizon[i], rank = r, sig = sig)
      ns <- c(ns, n)
    }
    if (!length(cols)) next
    fs <- intersect(FILTER_ORDER, unique(unlist(lapply(cols, function(c) names(c$rank)))))
    scen <- vapply(cols, `[[`, "", "scen"); hs <- vapply(cols, `[[`, 0, "h")
    groups <- rle(scen)
    starts <- cumsum(c(2L, head(groups$lengths, -1L)))
    head1 <- paste(c("", sprintf("\\multicolumn{%d}{c}{%s}", groups$lengths, gsub("_", " ", groups$values))), collapse = " & ")
    rules <- paste(sprintf("\\cmidrule(lr){%d-%d}", starts, starts + groups$lengths - 1L), collapse = " ")
    head2 <- paste(c("Filter", sprintf("$h=%d$", hs)), collapse = " & ")
    body <- vapply(fs, function(f) {
      cells <- vapply(cols, function(c) {
        v <- c$rank[f]
        if (!is.finite(v)) return("--")
        txt <- sprintf("%.2f", v)
        if (isTRUE(c$sig[f])) sprintf("\\textbf{%s}", txt) else txt
      }, "")
      lab <- flab(f)
      paste(c(lab, cells), collapse = " & ")
    }, "")
    body <- paste0(body, " \\\\")
    if (SPLIT_CAUSALITY) {
      nc <- length(cols) + 1L
      sec <- function(t) sprintf("\\midrule\n\\multicolumn{%d}{l}{\\textit{%s}} \\\\", nc, t)
      causal_f <- intersect(fs, CAUSAL_FILTERS); non_f <- setdiff(fs, c("non_filter", CAUSAL_FILTERS))
      body <- c(body[fs == "non_filter"], sec(paste(TXT$causal, "filters")), body[match(causal_f, fs)],
                sec(paste(TXT$noncausal, "filters")), body[match(non_f, fs)])
    }
    crit <- if (TABLE_SIG == "nemenyi") "its mean rank is lower than the baseline's by more than the Nemenyi critical difference" else
      sprintf("it outperforms the baseline (median $\\Delta$RMSE$<0$, paired Wilcoxon signed-rank test, %s $p \\le %.2f$)", mult_label(), ALPHA)
    caption <- sprintf("Mean Friedman ranks of the 16 methods for %s (lower is better; $N=%d$ series per column). Bold entries mark filters for which %s.",
                       model, min(ns), crit)
    tex <- c("% Gerado por gera_graficos.R - requer \\usepackage{booktabs} e \\usepackage{rotating}",
             "\\begin{sidewaystable}", "\\centering", "\\scriptsize", "\\setlength{\\tabcolsep}{3.5pt}",
             sprintf("\\caption{%s}", caption), sprintf("\\label{tab:ranks_%s}", tolower(model)),
             sprintf("\\begin{tabular}{l%s}", paste(rep("r", length(cols)), collapse = "")),
             "\\toprule", paste0(head1, " \\\\"), rules, paste0(head2, " \\\\"), "\\midrule",
             body, "\\bottomrule", "\\end{tabular}", "\\end{sidewaystable}")
    f <- file.path(TAB_DIR, sprintf("tab_ranks_%s.tex", model))
    writeLines(tex, f, useBytes = TRUE)
    message("  tabela: ", f)
    out_all[[model]] <- do.call(rbind, lapply(cols, function(c) data.frame(Model = model, Scenario = c$scen, Horizon = c$h,
                                Filter = names(c$rank), Mean_Rank = unname(c$rank), Bold = unname(c$sig[names(c$rank)]))))
  }
  if (length(out_all)) utils::write.csv(do.call(rbind, out_all), file.path(TAB_DIR, "ranks_friedman.csv"), row.names = FALSE)
}

# ---------------------------------- Main -------------------------------------
main <- function() {
  message("Lendo resultados em ", normalizePath(RESULTS_DIR))
  dat <- load_results()
  res <- dat$res; err <- dat$err; raw <- dat$raw
  dt <- delta_table(res)
  dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(dt, file.path(FIG_DIR, "tabela_delta_rmse_wilcoxon.csv"), row.names = FALSE)
  utils::write.csv(multiplicity_summary(dt), file.path(FIG_DIR, "resumo_multiplicidade.csv"), row.names = FALSE)
  fig_contamination()
  fig_cd(res)
  fig_heatmap(dt)
  fig_pareto(res, dt)
  fig_diagnostic(res, err)
  fig_dm(raw)
  write_rank_tables(res, dt)
  message("Pronto. Figuras em ", normalizePath(FIG_DIR))
}

if (sys.nframe() == 0L || identical(environment(), globalenv())) main()
