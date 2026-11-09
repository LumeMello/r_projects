# =============================================================================
# FUNCOES COMPARTILHADAS - IC M4 HOURLY / RESILIENCIA A CORRUPCAO
# Os scripts de modelo devem ser executados a partir desta pasta.
# =============================================================================

PRED_HORIZON <- 3L
INPUT_LEN <- 10L
N_SERIES <- 50L
N_WINDOWS <- 10L
N_VALIDATION_ORIGINS <- 3L
BASE_SEED <- 20261107L

# Periodo sazonal estimado por forecast::findfrequency na serie corrompida.
# Com ruido severo (spikes de 7 sd) o estimador devolve 1-4 em varias series
# (em Sev_Noise: 18 de 50), e esse valor alimenta os filtros seasonal e gabor.
# Abaixo de FREQ_MIN_VALID o periodo e tratado como "nao detectado" e substituido
# por FREQ_FALLBACK (24 = ciclo diario do M4 Hourly). Decisao metodologica:
# para reproduzir o comportamento anterior (usar o freq estimado sempre),
# defina FREQ_MIN_VALID <- 2L.
FREQ_MIN_VALID <- 12L
FREQ_FALLBACK <- 24L
M4_HOURLY_URL <- "https://raw.githubusercontent.com/M4competition/M4-methods/master/Dataset/Train/Hourly-train.csv"

filter_names <- c("non_filter", "ma", "ema", "LOESS", "LOWESS", "spline",
                  "winsor", "qes", "fft", "wavelet", "emd", "remd",
                  "hp", "kalman", "seasonal", "gabor")

noise_scenarios <- list(
  Clean      = list(noise_lvl = 0.00, outlier_ratio = 0.00, spike_mag = 0.0),
  Mod_Noise  = list(noise_lvl = 0.15, outlier_ratio = 0.03, spike_mag = 3.0),
  High_Noise = list(noise_lvl = 0.35, outlier_ratio = 0.08, spike_mag = 5.0),
  Sev_Noise  = list(noise_lvl = 0.50, outlier_ratio = 0.15, spike_mag = 7.0)
)

# Carrega dependencias e implementacoes locais das familias de filtros.
required_packages <- c("forecast", "writexl", "tspredit", "daltoolbox", "daltoolboxdp")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Instale os pacotes antes de executar: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(forecast)
  library(tspredit)
  library(daltoolbox)
  library(daltoolboxdp)
})

filter_sources <- c("ts_fil_fft.R", "ts_fil_remd.R", "ts_fil_spline.R", "ts_fil_gabor.R",
                    "ts_fil_wavelet.R")  # wavelet local: borda por reflexao (ver arquivo)
missing_filter_sources <- filter_sources[!file.exists(filter_sources)]
if (length(missing_filter_sources)) {
  stop("Arquivos locais de filtros ausentes: ", paste(missing_filter_sources, collapse = ", "))
}
for (f in filter_sources) source(f)
if (!requireNamespace("PMCMRplus", quietly = TRUE)) {
  install.packages("PMCMRplus")
}

get_filter_candidates <- function(j, freq) {
  switch(j,
         list(list(name = "non_filter", filter = NULL, desc = "none")),
         lapply(c(3, 5), function(k) list(name = "ma", filter = ts_fil_ma(k), desc = paste0("k=", k))),
         lapply(c(3, 5), function(k) list(name = "ema", filter = ts_fil_ema(k), desc = paste0("k=", k))),
         list(list(name = "LOESS", filter = ts_fil_smooth(), desc = "default")),
         lapply(c(0.1, 0.2, 0.4), function(f) list(name = "LOWESS", filter = ts_fil_lowess(f = f), desc = paste0("f=", f))),
         lapply(c(0.3, 0.5, 0.7), function(sp) list(name = "spline", filter = ts_fil_spline(spar = sp), desc = paste0("spar=", sp))),
         list(list(name = "winsor", filter = ts_fil_winsor(), desc = "default")),
         list(list(name = "qes", filter = ts_fil_qes(gamma = FALSE), desc = "gamma=F")),
         list(list(name = "fft", filter = ts_fil_fft(), desc = "default")),
         list(list(name = "wavelet", filter = ts_fil_wavelet(), desc = "default")),
         list(list(name = "emd", filter = ts_fil_emd(), desc = "default")),
         list(list(name = "remd", filter = ts_fil_remd(max_imfs = NULL), desc = "default")),
         list(list(name = "hp", filter = ts_fil_hp(), desc = "default")),
         lapply(c(0.01, 0.1, 1.0), function(h) list(name = "kalman", filter = ts_fil_kalman(H = h, Q = 1), desc = paste0("H=", h))),
         list(list(name = "seasonal", filter = ts_fil_seas_adj(frequency = max(2, freq)), desc = paste0("freq=", max(2, freq)))),
         list(list(name = "gabor", filter = ts_fil_gabor(central_frequency = 1 / max(2, freq)), desc = paste0("freq=", max(2, freq))))
  )
}

# Classificacao conservadora: deve ser refinada caso a implementacao do pacote
# comprove causalidade. "Potentially_non_causal" nao significa leakage no teste:
# indica que valores historicos podem ser revistos usando o restante do prefixo.
filter_causality <- c(
  non_filter = "causal_by_design",
  ma = "implementation_dependent",
  ema = "implementation_dependent",
  LOESS = "potentially_non_causal",
  LOWESS = "potentially_non_causal",
  spline = "potentially_non_causal",
  winsor = "implementation_dependent",
  qes = "implementation_dependent",
  fft = "potentially_non_causal",
  wavelet = "potentially_non_causal",
  emd = "potentially_non_causal",
  remd = "potentially_non_causal",
  hp = "potentially_non_causal",
  kalman = "implementation_dependent",
  seasonal = "potentially_non_causal",
  gabor = "potentially_non_causal"
)

apply_filter_candidate <- function(cand, series, max_edge_fill = 10L) {
  series <- as.numeric(series)
  if (is.null(cand$filter)) {
    return(list(series = series, status = "OK_NO_FILTER", error = NA_character_))
  }
  tryCatch({
    fitted <- fit(cand$filter, series)
    values <- as.numeric(transform(fitted, series))
    if (length(values) != length(series)) {
      stop("Saida do filtro com tamanho incorreto: ", length(values), " vs ", length(series))
    }
    bad <- which(!is.finite(values))
    status <- "OK"
    if (length(bad)) {
      leading <- identical(bad, seq_along(bad))   # so NAs contiguos no inicio
      if (!leading || length(bad) > max_edge_fill) {
        stop("Valores nao finitos fora da borda inicial (n=", length(bad), ").")
      }
      values[bad] <- series[bad]                  # fallback causal: valor bruto
      status <- "OK_EDGE_FILLED"
    }
    list(series = values, status = status, error = NA_character_)
  }, error = function(e) {
    list(series = NULL, status = "FILTER_FAILED", error = conditionMessage(e))
  })
}

corrupt_sensor_signal <- function(series, noise_lvl, outlier_ratio, spike_mag, seed) {
  set.seed(seed)
  n <- length(series)
  s <- stats::sd(series, na.rm = TRUE)
  if (!is.finite(s) || s == 0) s <- 0
  out <- series + stats::rnorm(n, mean = 0, sd = noise_lvl * s)
  if (outlier_ratio > 0 && s > 0) {
    n_out <- min(n, max(1L, floor(n * outlier_ratio)))
    idx <- sample(seq_len(n), size = n_out, replace = FALSE)
    spikes <- sample(c(-1, 1), n_out, replace = TRUE) * spike_mag * s
    out[idx] <- out[idx] + spikes
  }
  as.numeric(out)
}

load_m4_hourly <- function(n_series = N_SERIES) {
  message("Baixando M4 Hourly: ", M4_HOURLY_URL)
  dat <- utils::read.csv(M4_HOURLY_URL, stringsAsFactors = FALSE, check.names = FALSE)
  n <- min(as.integer(n_series), nrow(dat))
  ids <- as.character(dat[[1]][seq_len(n)])
  values <- lapply(seq_len(n), function(i) {
    x <- suppressWarnings(as.numeric(unlist(dat[i, -1, drop = FALSE], use.names = FALSE)))
    valid <- which(!is.na(x))
    if (!length(valid)) return(numeric(0))
    last <- max(valid)
    x <- x[seq_len(last)] # remove apenas padding a direita
    if (anyNA(x)) stop("NA interno na serie M4 ", ids[i], "; nao remova pontos e comprima o eixo temporal.")
    x
  })
  names(values) <- ids
  values
}

make_supervised_clean <- function(x_corrupted_filtered, y_clean, input_len) {
  stopifnot(length(x_corrupted_filtered) == length(y_clean))
  n <- length(x_corrupted_filtered)
  if (n <= input_len) stop("Historico insuficiente para janelas supervisionadas.")
  idx <- seq_len(n - input_len)
  X <- t(vapply(idx, function(i) x_corrupted_filtered[i:(i + input_len - 1L)], numeric(input_len)))
  Y <- as.numeric(y_clean[(input_len + 1L):n])
  list(X = X, Y = Y)
}

recursive_forecast <- function(model, last_window, horizon, predict_one) {
  window <- as.numeric(last_window)
  preds <- numeric(horizon)
  for (h in seq_len(horizon)) {
    preds[h] <- as.numeric(predict_one(model, matrix(window, nrow = 1L)))[1]
    window <- c(window[-1], preds[h])
  }
  preds
}

# ---------- ELM: implementacao explicita e reproduzivel ----------
fit_elm <- function(X, y, params, seed) {
  set.seed(seed)
  x_center <- colMeans(X); x_scale <- apply(X, 2, sd); x_scale[!is.finite(x_scale) | x_scale == 0] <- 1
  y_center <- mean(y); y_scale <- sd(y); if (!is.finite(y_scale) || y_scale == 0) y_scale <- 1
  Xs <- sweep(sweep(X, 2, x_center, "-"), 2, x_scale, "/")
  ys <- (y - y_center) / y_scale
  nhid <- as.integer(params$nhid)
  W <- matrix(stats::runif(ncol(Xs) * nhid, -1, 1), nrow = ncol(Xs), ncol = nhid)
  b <- stats::runif(nhid, -1, 1)
  H <- tanh(sweep(Xs %*% W, 2, b, "+"))
  lambda <- 1e-6
  beta <- tryCatch(solve(crossprod(H) + diag(lambda, nhid), crossprod(H, ys)), error = function(e) qr.solve(H, ys))
  structure(list(W = W, b = b, beta = beta, x_center = x_center, x_scale = x_scale,
                 y_center = y_center, y_scale = y_scale), class = "ic_elm")
}
predict_elm_one <- function(model, Xnew) {
  Xs <- sweep(sweep(Xnew, 2, model$x_center, "-"), 2, model$x_scale, "/")
  H <- tanh(sweep(Xs %*% model$W, 2, model$b, "+"))
  as.numeric(H %*% model$beta) * model$y_scale + model$y_center
}

# ---------- MLP: nnet com normalizacao ajustada somente no treino ----------
fit_mlp <- function(X, y, params, seed) {
  if (!requireNamespace("nnet", quietly = TRUE)) stop("Instale o pacote nnet.")
  set.seed(seed)
  x_center <- colMeans(X); x_scale <- apply(X, 2, sd); x_scale[!is.finite(x_scale) | x_scale == 0] <- 1
  y_center <- mean(y); y_scale <- sd(y); if (!is.finite(y_scale) || y_scale == 0) y_scale <- 1
  Xs <- sweep(sweep(X, 2, x_center, "-"), 2, x_scale, "/")
  ys <- (y - y_center) / y_scale
  model <- nnet::nnet(x = Xs, y = ys, size = as.integer(params$size),
                      decay = as.numeric(params$decay), maxit = as.integer(params$maxit),
                      linout = TRUE, trace = FALSE, MaxNWts = 100000)
  structure(list(model = model, x_center = x_center, x_scale = x_scale,
                 y_center = y_center, y_scale = y_scale), class = "ic_mlp")
}
predict_mlp_one <- function(model, Xnew) {
  Xs <- sweep(sweep(Xnew, 2, model$x_center, "-"), 2, model$x_scale, "/")
  as.numeric(stats::predict(model$model, Xs)) * model$y_scale + model$y_center
}

# ---------- LSTM: Keras 3 / backend TensorFlow ----------
fit_lstm <- function(X, y, params, seed) {
  if (!requireNamespace("keras3", quietly = TRUE)) stop("Instale keras3 e configure o backend TensorFlow (ver README).")
  set.seed(seed)
  keras3::set_random_seed(seed)
  x_center <- colMeans(X); x_scale <- apply(X, 2, sd); x_scale[!is.finite(x_scale) | x_scale == 0] <- 1
  y_center <- mean(y); y_scale <- sd(y); if (!is.finite(y_scale) || y_scale == 0) y_scale <- 1
  Xs <- sweep(sweep(X, 2, x_center, "-"), 2, x_scale, "/")
  X3 <- array(as.numeric(Xs), dim = c(nrow(Xs), ncol(Xs), 1L))
  ys <- as.numeric((y - y_center) / y_scale)
  model <- keras3::keras_model_sequential(input_shape = c(ncol(X), 1L)) |>
    keras3::layer_lstm(units = as.integer(params$units), dropout = as.numeric(params$dropout)) |>
    keras3::layer_dense(units = 1L)
  model |> keras3::compile(optimizer = keras3::optimizer_adam(learning_rate = as.numeric(params$learning_rate)),
                           loss = "mse")
  model |> keras3::fit(x = X3, y = ys, epochs = as.integer(params$epochs),
                       batch_size = min(as.integer(params$batch_size), nrow(X3)),
                       verbose = 0L)
  structure(list(model = model, x_center = x_center, x_scale = x_scale,
                 y_center = y_center, y_scale = y_scale), class = "ic_lstm")
}
predict_lstm_one <- function(model, Xnew) {
  Xs <- sweep(sweep(Xnew, 2, model$x_center, "-"), 2, model$x_scale, "/")
  X3 <- array(as.numeric(Xs), dim = c(nrow(Xs), ncol(Xs), 1L))
  as.numeric(stats::predict(model$model, X3, verbose = 0L)) * model$y_scale + model$y_center
}

fit_model_dispatch <- function(model_name, X, y, params, seed) {
  switch(model_name,
         ELM = fit_elm(X, y, params, seed),
         MLP = fit_mlp(X, y, params, seed),
         LSTM = fit_lstm(X, y, params, seed),
         stop("Modelo neural desconhecido: ", model_name)
  )
}
predict_model_dispatch <- function(model_name, model, Xnew) {
  switch(model_name,
         ELM = predict_elm_one(model, Xnew),
         MLP = predict_mlp_one(model, Xnew),
         LSTM = predict_lstm_one(model, Xnew),
         stop("Modelo neural desconhecido: ", model_name)
  )
}

fit_arima_model <- function(x) {
  forecast::auto.arima(stats::ts(as.numeric(x)), seasonal = FALSE,
                       stepwise = TRUE, approximation = FALSE)
}
predict_arima <- function(model, horizon) as.numeric(forecast::forecast(model, h = horizon)$mean)

metric_bundle <- function(actual, predicted) {
  e <- actual - predicted
  mape_ok <- abs(actual) > sqrt(.Machine$double.eps)
  c(RMSE = sqrt(mean(e^2)), MAE = mean(abs(e)),
    MAPE = if (any(mape_ok)) mean(abs(e[mape_ok] / actual[mape_ok])) else NA_real_,
    MAPE_n = sum(mape_ok))
}

# ---------- Teste de Diebold-Mariano ----------
# Os erros de teste sao N_WINDOWS janelas nao sobrepostas de PRED_HORIZON passos;
# dentro de cada janela os erros dos horizontes 1..h sao correlacionados, por isso
# h = PRED_HORIZON. O estimador padrao do forecast ("acf", o do teste original/HLN)
# soma as autocovariancias dos lags 1..h-1 com peso 1 e, com so 30 pontos, pode dar
# variancia NEGATIVA. Nesse caso o forecast avisa "Variance is negative ... Proceeding
# with horizon h=1" e devolve, sem marcar no resultado, um teste com h=1 (que ignora
# a autocorrelacao). Aqui: usa "acf" e, SO quando a variancia e negativa, refaz com
# "bartlett" (pesos 1 - j/h, sempre >= 0), mantendo h. Usar bartlett em todas as
# linhas encolhe a variancia e torna o teste menos conservador (no ARIMA mudou a
# conclusao de 52 de 3000 linhas), por isso fica so como fallback.
DM_VARESTIMATOR <- "acf"
DM_FALLBACK_VARESTIMATOR <- "bartlett"

dm_test_safe <- function(e_base, e_model, h = PRED_HORIZON, varestimator = DM_VARESTIMATOR,
                         fallback = DM_FALLBACK_VARESTIMATOR) {
  out <- list(stat = NA_real_, p = NA_real_, h_used = NA_integer_,
              varestimator = NA_character_, note = NA_character_)
  if (isTRUE(all.equal(e_base, e_model, tolerance = 1e-9))) {
    # Previsoes identicas ao baseline (ex.: LOESS em serie limpa, que nao
    # altera a serie): diferenca nula, sem evidencia de superioridade.
    out$stat <- 0; out$p <- 1; out$h_used <- as.integer(h); out$note <- "IDENTICAL_TO_BASELINE"
    return(out)
  }
  d <- e_base^2 - e_model^2
  if (!is.finite(stats::var(d)) || stats::var(d) <= 1e-12) {
    out$note <- "ZERO_VARIANCE_LOSS_DIFF"
    return(out)
  }
  run <- function(ve) {
    warns <- character(0)
    dm <- tryCatch(
      withCallingHandlers(
        forecast::dm.test(e_base, e_model, alternative = "greater", h = h, power = 2,
                          varestimator = ve),
        warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) structure(list(msg = conditionMessage(e)), class = "dm_error"))
    list(dm = dm, warns = warns)
  }
  r <- run(varestimator); used <- varestimator; note <- NA_character_
  negative <- any(grepl("Variance is negative", r$warns, fixed = TRUE))
  if (negative && !is.null(fallback) && fallback != varestimator) {
    r <- run(fallback); used <- fallback
    note <- paste0("NEGATIVE_VARIANCE_", toupper(varestimator), "_USED_", toupper(fallback))
  }
  if (inherits(r$dm, "dm_error")) { out$note <- paste0("DM_ERROR: ", r$dm$msg); return(out) }
  out$stat <- unname(as.numeric(r$dm$statistic))
  out$p <- as.numeric(r$dm$p.value)
  out$h_used <- as.integer(r$dm$parameter[["Forecast horizon"]])
  out$varestimator <- used
  if (length(r$warns)) note <- paste(c(na.omit(note), paste0("WARNING: ", unique(r$warns))), collapse = " | ")
  out$note <- note
  out
}

select_validation <- function(model_name, candidates, grid, x_train, y_train, freq, series_seed) {
  n <- length(x_train)
  fold_size <- PRED_HORIZON
  offsets <- rev(seq_len(N_VALIDATION_ORIGINS) * fold_size)
  origins <- n - offsets
  origins <- origins[origins > INPUT_LEN + 10L]
  if (!length(origins)) origins <- n - fold_size
  best <- list(score = Inf, candidate = candidates[[1]], params = grid[1, , drop = FALSE],
               filter_status = "NOT_TESTED", filter_error = NA_character_)
  
  for (cand in candidates) {
    for (g in seq_len(nrow(grid))) {
      params <- as.list(grid[g, , drop = FALSE])
      fold_scores <- numeric(0)
      fold_status <- character(0)
      fold_errors <- character(0)
      for (fold_id in seq_along(origins)) {
        cutoff <- origins[fold_id]
        x_prefix <- x_train[seq_len(cutoff)]
        y_prefix <- y_train[seq_len(cutoff)]
        target <- y_train[(cutoff + 1L):(cutoff + fold_size)]
        filtered <- apply_filter_candidate(cand, x_prefix)
        fold_status <- c(fold_status, filtered$status)
        fold_errors <- c(fold_errors, filtered$error)
        if (is.null(filtered$series)) next
        
        pred <- tryCatch({
          if (model_name == "ARIMA") {
            fit <- fit_arima_model(filtered$series)
            predict_arima(fit, fold_size)
          } else {
            sup <- make_supervised_clean(filtered$series, y_prefix, INPUT_LEN)
            fit <- fit_model_dispatch(model_name, sup$X, sup$Y, params,
                                      series_seed + 1000L * fold_id + 10L * g)
            recursive_forecast(fit, tail(filtered$series, INPUT_LEN), fold_size,
                               function(m, x) predict_model_dispatch(model_name, m, x))
          }
        }, error = function(e) NULL)
        if (!is.null(pred) && length(pred) == fold_size && all(is.finite(pred))) {
          fold_scores <- c(fold_scores, sqrt(mean((target - pred)^2)))
        }
      }
      # A configuracao so e elegivel se funcionou em todas as origens de validacao.
      if (length(fold_scores) == length(origins)) {
        score <- mean(fold_scores)
        if (is.finite(score) && score < best$score) {
          best <- list(score = score, candidate = cand, params = grid[g, , drop = FALSE],
                       filter_status = paste(unique(fold_status), collapse = ";"),
                       filter_error = paste(na.omit(unique(fold_errors)), collapse = " | "))
        }
      }
    }
  }
  best
}

# series_idx: indices (1..N_SERIES) a rodar; o padrao roda todas. Permite dividir o
# experimento em processos paralelos sem mudar as sementes (que dependem de i).
run_experiment <- function(model_name, grid, output_prefix, series_idx = NULL) {
  start_time <- Sys.time()
  series_list <- load_m4_hourly(N_SERIES)
  result_rows <- list(); counter <- 1L
  error_rows <- list(); error_counter <- 1L
  total_test_points <- N_WINDOWS * PRED_HORIZON
  if (is.null(series_idx)) series_idx <- seq_along(series_list)

  for (i in series_idx) {
    series_id <- names(series_list)[i]
    clean <- as.numeric(series_list[[i]])
    if (length(clean) < INPUT_LEN + total_test_points + 40L) {
      message("Pulando ", series_id, ": historico insuficiente.")
      next
    }
    for (scenario_name in names(noise_scenarios)) {
      scenario <- noise_scenarios[[scenario_name]]
      corrupted <- corrupt_sensor_signal(clean, scenario$noise_lvl, scenario$outlier_ratio,
                                         scenario$spike_mag, seed = BASE_SEED + i)
      train_end <- length(clean) - total_test_points
      x_initial <- corrupted[seq_len(train_end)]
      y_initial <- clean[seq_len(train_end)]
      freq_raw <- forecast::findfrequency(x_initial)
      freq <- if (is.finite(freq_raw) && freq_raw >= FREQ_MIN_VALID) as.integer(freq_raw) else
        max(2L, as.integer(FREQ_FALLBACK))
      baseline_errors <- NULL
      
      message("\n[", model_name, "] ", series_id, " | ", scenario_name)
      for (j in seq_along(filter_names)) {
        filter_family <- filter_names[j]
        candidates <- get_filter_candidates(j, freq)
        t0 <- Sys.time()
        selection <- select_validation(model_name, candidates, grid, x_initial, y_initial,
                                       freq, BASE_SEED + i * 10000L + j * 100L)
        best_cand <- selection$candidate
        best_params <- selection$params
        
        actual_all <- numeric(0); pred_all <- numeric(0)
        filter_statuses <- character(0); filter_errors <- character(0)
        rolling_ok <- TRUE; failure_reason <- NA_character_
        for (w in 0:(N_WINDOWS - 1L)) {
          cutoff <- train_end + w * PRED_HORIZON
          x_hist <- corrupted[seq_len(cutoff)]
          y_hist <- clean[seq_len(cutoff)]
          target <- clean[(cutoff + 1L):(cutoff + PRED_HORIZON)]
          filtered <- apply_filter_candidate(best_cand, x_hist)
          filter_statuses <- c(filter_statuses, filtered$status)
          filter_errors <- c(filter_errors, filtered$error)
          if (is.null(filtered$series)) {
            rolling_ok <- FALSE; failure_reason <- paste0("FILTER_FAILED: ", filtered$error); break
          }
          pred <- tryCatch({
            if (model_name == "ARIMA") {
              fit <- fit_arima_model(filtered$series)
              predict_arima(fit, PRED_HORIZON)
            } else {
              sup <- make_supervised_clean(filtered$series, y_hist, INPUT_LEN)
              fit <- fit_model_dispatch(model_name, sup$X, sup$Y, as.list(best_params),
                                        BASE_SEED + i * 100000L + j * 1000L + w)
              recursive_forecast(fit, tail(filtered$series, INPUT_LEN), PRED_HORIZON,
                                 function(m, x) predict_model_dispatch(model_name, m, x))
            }
          }, error = function(e) { failure_reason <<- conditionMessage(e); NULL })
          if (is.null(pred) || length(pred) != PRED_HORIZON || any(!is.finite(pred))) {
            rolling_ok <- FALSE
            if (is.na(failure_reason)) failure_reason <- "INVALID_PREDICTION"
            break
          }
          actual_all <- c(actual_all, target); pred_all <- c(pred_all, pred)
        }
        
        if (rolling_ok && length(pred_all) == total_test_points) {
          errors <- actual_all - pred_all
          metrics <- metric_bundle(actual_all, pred_all)
          status <- "CONVERGED"
          if (filter_family == "non_filter") baseline_errors <- errors
          dm_stat <- NA_real_; dm_p <- NA_real_; dm_h <- NA_integer_; dm_note <- NA_character_; dm_ve <- NA_character_
          if (filter_family == "non_filter") {
            dm_stat <- 0; dm_p <- 1; dm_h <- PRED_HORIZON; dm_note <- "BASELINE"; dm_ve <- DM_VARESTIMATOR
          } else if (!is.null(baseline_errors) && length(baseline_errors) == length(errors)) {
            dm <- dm_test_safe(baseline_errors, errors, h = PRED_HORIZON)
            dm_stat <- dm$stat; dm_p <- dm$p; dm_h <- dm$h_used; dm_note <- dm$note; dm_ve <- dm$varestimator
          } else {
            dm_note <- "NO_BASELINE"
          }
          error_rows[[error_counter]] <- data.frame(
            Model = model_name, Series_ID = series_id, Scenario = scenario_name,
            Filter = filter_family, Window = rep(0:(N_WINDOWS - 1L), each = PRED_HORIZON),
            Horizon = rep(seq_len(PRED_HORIZON), N_WINDOWS),
            Actual = actual_all, Predicted = pred_all, Error = errors, stringsAsFactors = FALSE)
          error_counter <- error_counter + 1L
        } else {
          status <- "FAILED_ROLLING"
          metrics <- c(RMSE = NA_real_, MAE = NA_real_, MAPE = NA_real_, MAPE_n = 0)
          dm_stat <- NA_real_; dm_p <- NA_real_; dm_h <- NA_integer_; dm_note <- NA_character_; dm_ve <- NA_character_
        }
        
        params_text <- paste(names(best_params), unlist(best_params), sep = "=", collapse = ";")
        row <- data.frame(
          Model = model_name, Series_ID = series_id, Scenario = scenario_name,
          Noise_Level = scenario$noise_lvl, Outlier_Ratio = scenario$outlier_ratio,
          Filter = filter_family, Filter_Param = best_cand$desc,
          Filter_Causality = unname(filter_causality[filter_family]),
          Freq_Raw = freq_raw, Freq_Used = freq,
          Filter_Status = paste(unique(filter_statuses), collapse = ";"),
          Filter_Error = paste(na.omit(unique(c(filter_errors, selection$filter_error))), collapse = " | "),
          Selection_Validation_RMSE = selection$score, Selected_Params = params_text,
          Status = status, Failure_Reason = failure_reason,
          RMSE = unname(metrics["RMSE"]), MAE = unname(metrics["MAE"]),
          MAPE = unname(metrics["MAPE"]), MAPE_n = unname(metrics["MAPE_n"]),
          DM_Stat = dm_stat, DM_Pvalue = dm_p, DM_H_Used = dm_h,
          DM_Varestimator = dm_ve, DM_Note = dm_note,
          Runtime_Seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
          stringsAsFactors = FALSE
        )
        result_rows[[counter]] <- row; counter <- counter + 1L
        message(sprintf("  %s | RMSE=%s | %s", filter_family,
                        ifelse(is.finite(row$RMSE), round(row$RMSE, 4), "NA"), status))
      }
      if (counter %% 20L == 0L) saveRDS(result_rows, sprintf("checkpoint_%s.rds", model_name))
    }
  }
  
  if (!length(result_rows)) stop("Nenhum resultado foi produzido.")
  results <- do.call(rbind, result_rows)
  summary_rows <- list(); s <- 1L
  for (sc in unique(results$Scenario)) {
    d <- results[results$Scenario == sc & results$Status == "CONVERGED" & is.finite(results$RMSE), ]
    wide <- reshape(d[, c("Series_ID", "Filter", "RMSE")], idvar = "Series_ID",
                    timevar = "Filter", direction = "wide")
    names(wide) <- sub("^RMSE\\.", "", names(wide))
    required <- intersect(filter_names, names(wide))
    rownames(wide) <- as.character(wide$Series_ID)
    mat <- as.matrix(wide[, required, drop = FALSE])
    complete <- stats::complete.cases(mat)
    mat <- mat[complete, , drop = FALSE]
    if (nrow(mat) >= 2L && ncol(mat) >= 3L) {
      ft <- tryCatch(stats::friedman.test(mat), error = function(e) NULL)
      ranks <- if (nrow(mat)) colMeans(t(apply(mat, 1, rank, ties.method = "average"))) else rep(NA_real_, ncol(mat))
      names(ranks) <- colnames(mat)
      nemenyi_p <- rep(NA_real_, ncol(mat)); names(nemenyi_p) <- colnames(mat)
      if (!is.null(ft) && is.finite(ft$p.value) && ft$p.value < 0.05 && requireNamespace("PMCMRplus", quietly = TRUE)) {
        long <- data.frame(Series = rep(rownames(mat), times = ncol(mat)),
                           Filter = rep(colnames(mat), each = nrow(mat)),
                           RMSE = as.vector(mat))
        nm <- tryCatch(PMCMRplus::frdAllPairsNemenyiTest(y = long$RMSE,
                                                         groups = long$Filter, blocks = long$Series), error = function(e) NULL)
        if (!is.null(nm)) {
          pm <- nm$p.value
          # A matriz e triangular inferior: o par (a, b) fica em [a, b] OU em [b, a],
          # e a outra celula e NA. Ler so uma das duas deixava NA para os filtros
          # que sao, ao mesmo tempo, linha e coluna (spline, qes, wavelet, remd, seasonal).
          get_p <- function(a, b) {
            v <- c(if (a %in% rownames(pm) && b %in% colnames(pm)) pm[a, b],
                   if (b %in% rownames(pm) && a %in% colnames(pm)) pm[b, a])
            v <- v[is.finite(v)]
            if (length(v)) v[1L] else NA_real_
          }
          for (f in names(nemenyi_p)) {
            nemenyi_p[f] <- if (f == "non_filter") 1 else get_p("non_filter", f)
          }
        }
      }
      wilcox_p <- rep(NA_real_, ncol(mat)); names(wilcox_p) <- colnames(mat)
      if ("non_filter" %in% colnames(mat)) {
        base <- mat[, "non_filter"]
        raw_p <- sapply(setdiff(colnames(mat), "non_filter"), function(f) {
          tryCatch(stats::wilcox.test(mat[, f], base, paired = TRUE, exact = FALSE)$p.value,
                   error = function(e) NA_real_)
        })
        wilcox_p[names(raw_p)] <- stats::p.adjust(raw_p, method = "holm")
        wilcox_p["non_filter"] <- 1
      }
      for (f in colnames(mat)) {
        summary_rows[[s]] <- data.frame(Scenario = sc, Filter = f,
                                        Mean_Rank = unname(ranks[f]), P_Nemenyi_vs_NoFilter = unname(nemenyi_p[f]),
                                        P_Wilcoxon_vs_NoFilter_Holm = unname(wilcox_p[f]),
                                        P_vs_NoFilter = unname(wilcox_p[f]),
                                        Significante = ifelse(is.finite(wilcox_p[f]) && wilcox_p[f] <= 0.05, "SIM (Holm <= 0.05)",
                                                              ifelse(is.na(wilcox_p[f]), "NA", "NAO")),
                                        N_Complete_Series = nrow(mat), Friedman_Statistic = if (is.null(ft)) NA_real_ else unname(ft$statistic),
                                        Friedman_Pvalue = if (is.null(ft)) NA_real_ else ft$p.value,
                                        Filter_Causality = unname(filter_causality[f]), stringsAsFactors = FALSE)
        s <- s + 1L
      }
    } else {
      for (f in filter_names) {
        summary_rows[[s]] <- data.frame(Scenario = sc, Filter = f, Mean_Rank = NA_real_,
                                        P_Nemenyi_vs_NoFilter = NA_real_, P_Wilcoxon_vs_NoFilter_Holm = NA_real_,
                                        P_vs_NoFilter = if (f == "non_filter") 1 else NA_real_, Significante = "NA",
                                        N_Complete_Series = nrow(mat), Friedman_Statistic = NA_real_, Friedman_Pvalue = NA_real_,
                                        Filter_Causality = unname(filter_causality[f]), stringsAsFactors = FALSE)
        s <- s + 1L
      }
    }
  }
  summary_df <- do.call(rbind, summary_rows)
  out_file <- sprintf("%s_industrial_sensors_results_%dpts.xlsx", output_prefix, total_test_points)
  # Erros ponto a ponto do teste: permitem refazer o DM (ou outro teste) sem rodar de novo.
  errors_df <- if (length(error_rows)) do.call(rbind, error_rows) else data.frame()
  writexl::write_xlsx(list(Resultados_Gerais = results, Resumo_Estatistico = summary_df,
                           Erros_Teste = errors_df), out_file)
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  message("Concluido: ", out_file, " (", round(elapsed, 2), " min)")
  invisible(list(results = results, summary = summary_df, errors = errors_df, file = out_file))
}