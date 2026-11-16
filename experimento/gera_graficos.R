# =============================================================================
# FIGURAS DO ARTIGO - IC M4 HOURLY / RESILIENCIA A CORRUPCAO
# Rodar na mesma pasta dos resultados (<MODELO>_industrial_sensors_results_<N>pts.xlsx)
# e do IC_Common.R. Gera, em FIG_DIR, PNG (300 dpi) e PDF de:
#   1. Diagramas de Diferenca Critica (Demsar / Nemenyi)
#   2. Heatmap de Delta RMSE% mediano vs. baseline (asterisco = Wilcoxon-Holm)
#   3. Fronteira de Pareto: acuracia vs. custo computacional
#   4. Painel diagnostico de degradacao: (a) Error Ratio por passo, (b) 48 h sobrepostas
#   5. Barras empilhadas do teste Diebold-Mariano serie a serie
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
MODELS <- c("ARIMA", "MLP", "ELM")
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
HEAT_LIMIT <- 50                  # Delta RMSE% no extremo da escala de cor (valores acima saturam)

FILTER_ORDER <- c("non_filter", "ma", "ema", "LOESS", "LOWESS", "spline", "winsor", "qes",
                  "fft", "wavelet", "emd", "remd", "hp", "kalman", "seasonal", "gabor")
FILTER_LABEL <- c(non_filter = "Baseline", ma = "MA", ema = "EMA", LOESS = "Smooth (LOESS)",
                  LOWESS = "LOWESS", spline = "Spline", winsor = "Winsor", qes = "QES",
                  fft = "FFT", wavelet = "Wavelet", emd = "EMD", remd = "REMD", hp = "HP",
                  kalman = "Kalman", seasonal = "Seasonal", gabor = "Gabor")
SCENARIO_ORDER <- c("Clean", "Mod_Noise", "High_Noise", "Sev_Noise")

# Textos (traduza aqui se o artigo for em ingles)
TXT <- list(
  cd_title = "Diagrama de Diferenca Critica",
  rank = "Rank medio (Friedman)",
  delta = "Delta RMSE% mediano vs. baseline",
  delta_short = "Delta RMSE% (mediana)",
  heat_caption = "* Wilcoxon pareado vs. baseline com correcao de Holm, p <= 0,05. Azul = reducao do erro.",
  overhead = "Overhead do filtro (ms, escala log)",
  incremental = "Atualizacao",
  frontier = "Fronteira de Pareto",
  step = "Passo dentro da janela de previsao",
  ratio = "Error Ratio (RMSE filtro / RMSE baseline)",
  hours = "Horas (ultimas 48 h do treino)",
  value = "Valor",
  dm_better = "Significativamente melhor",
  dm_equal = "Indistinguivel",
  dm_worse = "Significativamente pior",
  dm_share = "Proporcao das series",
  dm_caption = "Diebold-Mariano (perda quadratica) vs. baseline, p < 0,05 unilateral em cada direcao.",
  clean = "Serie limpa", corrupted = "Serie corrompida", filtered = "Filtrada"
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
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = width, height = height, bg = "white")
  message("  figura: ", file.path(FIG_DIR, paste0(name, ".png")))
}

flab <- function(f) unname(ifelse(f %in% names(FILTER_LABEL), FILTER_LABEL[f], f))
scen_lab <- function(s) gsub("_", " ", s)

# ------------------------------ Leitura dos dados ----------------------------
load_results <- function() {
  pat <- sprintf("^(%s)_industrial_sensors_results_([0-9]+)pts\\.xlsx$", paste(MODELS, collapse = "|"))
  files <- list.files(RESULTS_DIR, pattern = pat, full.names = TRUE)
  if (!length(files)) stop("Nenhum arquivo de resultados encontrado em ", normalizePath(RESULTS_DIR))
  res <- list(); err <- list()
  for (f in files) {
    m <- regmatches(basename(f), regexec(pat, basename(f)))[[1]]
    model <- m[2]; h <- as.integer(m[3]) %/% N_WINDOWS_RESULTS
    r <- as.data.frame(readxl::read_excel(f, sheet = "Resultados_Gerais"))
    r$Model <- model; r$Horizon <- h
    res[[f]] <- r
    if ("Erros_Teste" %in% readxl::excel_sheets(f)) {
      e <- as.data.frame(readxl::read_excel(f, sheet = "Erros_Teste"))
      if (nrow(e)) {
        names(e)[names(e) == "Horizon"] <- "Step"
        e$Model <- model; e$Horizon <- h
        err[[f]] <- e[, c("Model", "Horizon", "Series_ID", "Scenario", "Filter", "Window", "Step", "Error")]
      }
    } else {
      message("  ", basename(f), ": sem aba Erros_Teste (Error Ratio por passo indisponivel para esse arquivo)")
    }
    message("  lido: ", basename(f), " (", model, ", h = ", h, ", ", nrow(r), " linhas)")
  }
  cols <- c("Model", "Horizon", "Series_ID", "Scenario", "Filter", "Filter_Param",
            "Status", "RMSE", "DM_Stat", "DM_Pvalue")
  res <- do.call(rbind, lapply(res, function(r) r[, cols]))
  res <- res[res$Status == "CONVERGED" & is.finite(res$RMSE), ]
  list(res = res, err = if (length(err)) do.call(rbind, err) else NULL)
}

# Matriz series x filtros (so series completas, como no IC_Common.R)
rmse_matrix <- function(res, model, h, scen) {
  d <- res[res$Model == model & res$Horizon == h & res$Scenario == scen, c("Series_ID", "Filter", "RMSE")]
  if (!nrow(d)) return(NULL)
  w <- reshape(d, idvar = "Series_ID", timevar = "Filter", direction = "wide")
  names(w) <- sub("^RMSE\\.", "", names(w))
  fs <- intersect(FILTER_ORDER, names(w))
  m <- as.matrix(w[, fs, drop = FALSE]); rownames(m) <- w$Series_ID
  m[stats::complete.cases(m), , drop = FALSE]
}

combos <- function(res) {
  u <- unique(res[, c("Model", "Horizon", "Scenario")])
  u <- u[order(match(u$Model, MODELS), u$Horizon, match(u$Scenario, SCENARIO_ORDER)), ]
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
    labs(title = title, subtitle = sprintf("N = %d series, k = %d metodos; barras ligam metodos sem diferenca (Nemenyi, alfa = %.2f)",
                                           n, k, ALPHA)) +
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
                           P_Holm = stats::p.adjust(raw_p, "holm"), N = nrow(m))
  }
  d <- do.call(rbind, out); rownames(d) <- NULL; d
}

fig_heatmap <- function(dt) {
  message("2. Heatmap de Delta RMSE%")
  d <- dt
  d$X <- factor(paste0("h=", d$Horizon), levels = paste0("h=", sort(unique(d$Horizon))))
  d$Scenario <- factor(scen_lab(d$Scenario), levels = scen_lab(SCENARIO_ORDER))
  d$Y <- factor(flab(d$Filter), levels = rev(flab(setdiff(FILTER_ORDER, "non_filter"))))
  d$Model <- factor(d$Model, levels = intersect(MODELS, d$Model))
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
                         labels = function(x) ifelse(abs(x) == lim, paste0(ifelse(x > 0, ">= +", "<= "), x, "%"),
                                                     paste0(ifelse(x > 0, "+", ""), x, "%")),
                         name = TXT$delta_short,
                         guide = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
    facet_grid(Model ~ Scenario, scales = "free_x", space = "free_x", switch = "x") +
    labs(title = TXT$delta, x = NULL, y = NULL, caption = TXT$heat_caption) +
    theme_article() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(size = 7.5),
          strip.placement = "outside", panel.spacing.x = unit(0.35, "lines"),
          legend.key.width = unit(1.6, "cm"), legend.key.height = unit(0.3, "cm"))
  save_fig(p, "2_heatmap_delta_rmse", 7.2, 1.2 + 3.1 * nlevels(d$Model))
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
  one <- function(ycol, ylab, ref, file, title) {
    d$Y <- d[[ycol]]
    ref$Model <- factor(ref$Model, levels = levels(d$Model))
    d <- d[is.finite(d$Y), ]
    d$front <- FALSE
    for (mm in levels(d$Model)) { ii <- d$Model == mm; d$front[ii] <- pareto_front(d$Overhead_ms[ii], d$Y[ii]) }
    fr <- d[d$front, ]; fr <- fr[order(fr$Model, fr$Overhead_ms), ]
    p <- ggplot(d, aes(Overhead_ms, Y)) +
      geom_hline(data = ref, aes(yintercept = Y), linetype = "dashed", colour = COL$ink2, linewidth = 0.4) +
      geom_step(data = fr, aes(group = Model), direction = "hv", colour = COL$ink, linewidth = 0.5, alpha = 0.6) +
      geom_point(aes(fill = Update, shape = Update), size = 2.8, colour = "white", stroke = 0.6) +
      geom_text_repel(aes(label = Label, fontface = ifelse(front, "bold", "plain")), size = 2.5,
                      colour = COL$ink, min.segment.length = 0.2, segment.colour = COL$ink2,
                      segment.size = 0.25, max.overlaps = 30, seed = 1) +
      scale_x_log10(labels = function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE, drop0trailing = TRUE)) +
      scale_fill_manual(values = c("O(1)" = COL$o1, "Batch" = COL$batch), name = TXT$incremental, drop = FALSE) +
      scale_shape_manual(values = c("O(1)" = 21, "Batch" = 24), name = TXT$incremental, drop = FALSE) +
      facet_wrap(~ Model, nrow = 1) +
      labs(title = title, subtitle = sprintf("Cenario %s, media sobre h = %s. Tracejado = baseline (sem filtro).\nLinha preta: %s (menor custo e menor valor).",
                                             scen_lab(PARETO_SCENARIO), hs, TXT$frontier),
           x = TXT$overhead, y = ylab) +
      theme_article()
    save_fig(p, file, max(5, 3.2 * nlevels(d$Model) + 0.8), 4.6)
  }
  one("Rank", TXT$rank, setNames(base_rank, c("Model", "Y")), "3_pareto_rank", "Acuracia vs. custo computacional")
  one("Delta", TXT$delta_short, data.frame(Model = levels(d$Model), Y = 0), "3_pareto_delta_rmse",
      "Acuracia vs. custo computacional")
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

overlay_data <- function(res, model) {
  if (!load_ic_common()) return(NULL)
  sl <- load_m4_hourly(N_SERIES)
  n_test <- N_WINDOWS * PRED_HORIZON
  sc <- noise_scenarios[[DIAG_SCENARIO_OVERLAY]]
  ids <- if (is.null(DIAG_SERIES)) names(sl) else DIAG_SERIES
  pick <- NULL
  for (id in ids) {
    cl <- as.numeric(sl[[id]]); cl <- cl[seq_len(length(cl) - n_test)]
    if (is.null(DIAG_SERIES) && forecast::findfrequency(cl) != 24) next
    pick <- id; break
  }
  if (is.null(pick)) pick <- names(sl)[1]
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
        labs(title = sprintf("(a) %s: Error Ratio por passo (h = %d)", model, er$Horizon[1]),
             subtitle = "Mediana entre series do RMSE do filtro / RMSE do baseline em cada passo; > 1 = pior que o baseline (escala log).",
             x = TXT$step, y = TXT$ratio) +
        theme_article()
    } else {
      message("  ", model, ": sem erros por passo (aba Erros_Teste); painel (a) pulado")
    }
    ov <- overlay_data(res, model)
    pb <- NULL
    if (!is.null(ov) && nrow(ov)) {
      ov$F <- factor(sprintf("%s (%s)\nciclo 24 h retido: %.0f%%", flab(ov$Filter), ov$Param, ov$Amp),
                     levels = unique(sprintf("%s (%s)\nciclo 24 h retido: %.0f%%", flab(ov$Filter), ov$Param, ov$Amp)[order(match(ov$Filter, fs))]))
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
        labs(title = sprintf("(b) Ultimas 48 h: serie %s, cenario %s", ov$Series_ID[1], scen_lab(DIAG_SCENARIO_OVERLAY)),
             subtitle = "Linha colorida = filtrada. 'Ciclo 24 h retido' = amplitude da senoide de 24 h na serie filtrada / na serie limpa.",
             x = TXT$hours, y = TXT$value) +
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
  d$Class <- ifelse(!is.finite(d$DM_Pvalue) | identical0, TXT$dm_equal,
                    ifelse(d$DM_Pvalue < ALPHA, TXT$dm_better,
                           ifelse(1 - d$DM_Pvalue < ALPHA, TXT$dm_worse, TXT$dm_equal)))
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
    tab$Scenario <- factor(scen_lab(as.character(tab$Scenario)), levels = scen_lab(SCENARIO_ORDER))
    tab$Model <- factor(tab$Model, levels = intersect(MODELS, unique(tab$Model)))
    p <- ggplot(tab, aes(Share, F, fill = Class)) +
      geom_col(width = 0.75, colour = "white", linewidth = 0.4) +
      scale_fill_manual(values = setNames(c(COL$bad, COL$neutral, COL$good), levels(d$Class)), name = NULL,
                        breaks = c(TXT$dm_better, TXT$dm_equal, TXT$dm_worse)) +
      scale_x_continuous(labels = scales::percent, breaks = c(0.25, 0.5, 0.75), expand = c(0, 0)) +
      facet_grid(Model ~ Scenario) +
      labs(title = sprintf("Teste Diebold-Mariano serie a serie (h = %d)", h),
           subtitle = sprintf("Proporcao das %d series em cada classe; filtros ordenados pela proporcao media de melhora.",
                              max(tab$N)),
           x = TXT$dm_share, y = NULL, caption = TXT$dm_caption) +
      theme_article() + theme(panel.grid.major.y = element_blank(), panel.spacing.x = unit(0.8, "lines"))
    save_fig(p, sprintf("5_dm_series_h%02d", h), 9, 1.4 + 2.9 * nlevels(tab$Model))
  }
}

# ---------------------------------- Main -------------------------------------
main <- function() {
  message("Lendo resultados em ", normalizePath(RESULTS_DIR))
  dat <- load_results()
  res <- dat$res; err <- dat$err
  dt <- delta_table(res)
  dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(dt, file.path(FIG_DIR, "tabela_delta_rmse_wilcoxon.csv"), row.names = FALSE)
  fig_cd(res)
  fig_heatmap(dt)
  fig_pareto(res, dt)
  fig_diagnostic(res, err)
  fig_dm(res)
  message("Pronto. Figuras em ", normalizePath(FIG_DIR))
}

if (sys.nframe() == 0L || identical(environment(), globalenv())) main()
