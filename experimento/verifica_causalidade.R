# =============================================================================
# TESTE EMPIRICO DE CAUSALIDADE DOS FILTROS
# Para cada filtro (todos os parametros candidatos), aplica o filtro ao prefixo
# x[1:n] e ao prefixo estendido x[1:(n + EXTRA)] e compara os n primeiros valores.
#   - Se nada muda: o filtro e causal (a saida em t so depende de x[1:t]).
#   - Se muda: o filtro revisa o passado quando chegam dados novos (nao causal).
# Tambem mede quanto o ULTIMO ponto do prefixo (a origem da previsao) muda
# quando os EXTRA pontos seguintes ficam disponiveis, em unidades do desvio
# padrao da serie. Roda na pasta do IC_Common.R. Saida: causalidade_filtros.csv
# =============================================================================
suppressMessages(source("IC_Common.R"))

CAUSAL_SERIES <- 10L          # series do M4 usadas no teste
CAUSAL_SCENARIO <- "High_Noise"
EXTRA <- 24L                  # observacoes acrescentadas ao prefixo
TOL <- 1e-8                   # diferenca relativa considerada zero

sl <- load_m4_hourly(CAUSAL_SERIES)
sc <- noise_scenarios[[CAUSAL_SCENARIO]]
n_test <- N_WINDOWS * PRED_HORIZON
rows <- list()
for (j in seq_along(filter_names)[-1]) {
  for (cand in get_filter_candidates(j, 24L)) {
    for (i in seq_along(sl)) {
      cl <- as.numeric(sl[[i]])
      x <- corrupt_sensor_signal(cl, sc$noise_lvl, sc$outlier_ratio, sc$spike_mag, BASE_SEED + i)
      n <- length(cl) - n_test - EXTRA
      s <- stats::sd(x[seq_len(n)])
      set.seed(1); a <- apply_filter_candidate(cand, x[seq_len(n)])
      set.seed(1); b <- apply_filter_candidate(cand, x[seq_len(n + EXTRA)])
      if (is.null(a$series) || is.null(b$series)) next
      d <- abs(a$series - b$series[seq_len(n)]) / s
      rows[[length(rows) + 1L]] <- data.frame(
        Filter = filter_names[j], Param = cand$desc, Series_ID = names(sl)[i],
        Max_Change_sd = max(d), Last_Point_Change_sd = d[n],
        Changed_Points = sum(d > TOL), First_Changed_From_End = if (any(d > TOL)) n - min(which(d > TOL)) + 1L else 0L)
    }
  }
}
r <- do.call(rbind, rows)
agg <- do.call(rbind, lapply(split(r, list(r$Filter, r$Param), drop = TRUE), function(d) data.frame(
  Filter = d$Filter[1], Param = d$Param[1],
  Causal = all(d$Max_Change_sd <= TOL),
  Median_Last_Point_Change_sd = stats::median(d$Last_Point_Change_sd),
  Median_Max_Change_sd = stats::median(d$Max_Change_sd),
  Median_Revised_Points = stats::median(d$Changed_Points),
  Median_Revision_Depth = stats::median(d$First_Changed_From_End))))
agg <- agg[order(match(agg$Filter, filter_names), agg$Param), ]
rownames(agg) <- NULL
fam <- aggregate(Causal ~ Filter, agg, all)
fam <- fam[order(match(fam$Filter, filter_names)), ]
utils::write.csv(agg, "causalidade_filtros.csv", row.names = FALSE)
cat("\nPor parametro (mudancas em desvios padrao da serie; profundidade = quantos pontos para tras sao revistos):\n")
print(format(agg, digits = 3), row.names = FALSE)
cat("\nPor familia (causal so se todos os parametros forem causais):\n")
print(fam, row.names = FALSE)
cat("\nSugestao para filter_causality no IC_Common.R:\n")
cat(sprintf('  %s = "%s"', fam$Filter, ifelse(fam$Causal, "causal", "non_causal")), sep = ",\n"); cat("\n")
