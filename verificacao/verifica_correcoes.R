# Verificacao rapida das correcoes (rodar na pasta com IC_Common_2.R e os ts_fil_*.R).
# Nao roda o experimento: so confere wavelet e gabor em series reais do M4.
source("IC_Common_2.R")
sl <- load_m4_hourly(N_SERIES)
n_test <- N_WINDOWS * PRED_HORIZON

cat("\n=== WAVELET: desvio no ULTIMO ponto (filtrado - bruto) ===\n")
for (id in c("H25", "H28", "H30")) {
  x <- as.numeric(sl[[id]]); x <- x[seq_len(length(x) - n_test)]; n <- length(x)
  wt <- wavelets::modwt(x, filter = "haar", boundary = "periodic")
  wt@W[[1]] <- matrix(0, nrow(wt@W[[1]]), 1)
  old <- as.numeric(wavelets::imodwt(wt))
  new <- as.numeric(transform(fit(ts_fil_wavelet(), x), x))
  cat(sprintf("%s  x[1]=%.0f | periodico: %8.1f (formula %8.1f) | reflexao: %6.1f (formula %6.1f) | len=%d/%d\n",
              id, x[1], old[n] - x[n], (x[n-1] + x[1] - 2*x[n]) / 4,
              new[n] - x[n], (x[n-1] - x[n]) / 4, length(new), n))
}

cat("\n=== GABOR: guarda de ganho (esperado: freq 2,3,4,6,8 recusados; 12 e 24 ok) ===\n")
x <- as.numeric(sl[["H48"]]); x <- x[seq_len(length(x) - n_test)]
for (fr in c(2, 3, 4, 6, 8, 12, 24)) {
  r <- apply_filter_candidate(list(filter = ts_fil_gabor(central_frequency = 1 / fr)), x)
  cat(sprintf("freq=%2d  %-13s  max|f|/max|x| = %s  %s\n", fr, r$status,
              if (is.null(r$series)) "NA" else format(max(abs(r$series)) / max(abs(x)), digits = 3),
              if (is.na(r$error)) "" else substr(r$error, 1, 60)))
}

cat("\n=== FREQ: piso aplicado (FREQ_MIN_VALID =", FREQ_MIN_VALID, ", fallback =", FREQ_FALLBACK, ") ===\n")
for (id in c("H48", "H15", "H19", "H12")) {
  i <- match(id, names(sl)); cl <- as.numeric(sl[[id]])
  sc <- noise_scenarios$Sev_Noise
  co <- corrupt_sensor_signal(cl, sc$noise_lvl, sc$outlier_ratio, sc$spike_mag, BASE_SEED + i)
  fr <- forecast::findfrequency(co[seq_len(length(cl) - n_test)])
  cat(sprintf("%s  freq_raw=%d -> usado=%d\n", id, fr, if (fr >= FREQ_MIN_VALID) fr else FREQ_FALLBACK))
}
