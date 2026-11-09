#'@title FFT Filter with Smooth Roll-off and Boundary Protection
#'@description Filtro passa-baixa no dominio da frequencia com padding antissimetrico.
#'@param threshold Fracao acumulada da energia nao-DC a reter, em (0, 1].
#'@param pad_pct Fracao do comprimento usada no padding, em [0, 1].
#'@param taper_pct Largura relativa da banda de transicao, positiva.
#'@return Um objeto `ts_fil_fft`.
#'@importFrom daltoolbox dal_transform
#'@importFrom daltoolbox fit
#'@importFrom daltoolbox transform
#'@export
ts_fil_fft <- function(threshold = 0.90, pad_pct = 0.25, taper_pct = 0.15) {
  if (length(threshold) != 1L || !is.finite(threshold) || threshold <= 0 || threshold > 1)
    stop("threshold deve estar no intervalo (0, 1].")
  if (length(pad_pct) != 1L || !is.finite(pad_pct) || pad_pct < 0 || pad_pct > 1)
    stop("pad_pct deve estar no intervalo [0, 1].")
  if (length(taper_pct) != 1L || !is.finite(taper_pct) || taper_pct <= 0)
    stop("taper_pct deve ser positivo.")
  obj <- daltoolbox::dal_transform()
  obj$threshold <- threshold
  obj$pad_pct <- pad_pct
  obj$taper_pct <- taper_pct
  obj$cutindex <- NULL
  obj$pad_len <- NULL
  class(obj) <- append("ts_fil_fft", class(obj))
  obj
}

compute_cut_index <- function(freqs, threshold = 0.90) {
  energy <- as.numeric(freqs)[-1L] # exclui DC
  total <- sum(energy)
  if (!is.finite(total) || total <= .Machine$double.eps) return(1L)
  idx <- which(cumsum(energy) / total >= threshold)[1L]
  if (is.na(idx)) idx <- length(energy)
  as.integer(idx + 1L)
}

fft_pad_signal <- function(x, pad_len) {
  n <- length(x)
  if (pad_len == 0L) return(x)
  left <- 2 * x[1L] - x[(pad_len + 1L):2L]
  right <- 2 * x[n] - x[(n - 1L):(n - pad_len)]
  c(left, x, right)
}

#'@importFrom daltoolbox fit
#'@exportS3Method fit ts_fil_fft
fit.ts_fil_fft <- function(obj, data, ...) {
  x <- as.numeric(data)
  n <- length(x)
  if (n < 4L) stop("FFT requer pelo menos quatro observacoes.")
  if (any(!is.finite(x))) stop("FFT requer entrada finita e sem NA.")
  pad_len <- min(max(0L, floor(n * obj$pad_pct)), n - 1L)
  padded <- fft_pad_signal(x, pad_len)
  spec <- Mod(stats::fft(padded))^2
  half_len <- floor(length(padded) / 2) + 1L
  obj$cutindex <- compute_cut_index(spec[seq_len(half_len)], obj$threshold)
  obj$pad_len <- pad_len
  obj$fit_length <- n
  obj
}

#'@importFrom daltoolbox transform
#'@exportS3Method transform ts_fil_fft
transform.ts_fil_fft <- function(obj, data, ...) {
  original <- data
  x <- as.numeric(data)
  n <- length(x)
  if (n < 4L) stop("FFT requer pelo menos quatro observacoes.")
  if (any(!is.finite(x))) stop("FFT requer entrada finita e sem NA.")
  if (!is.null(obj$fit_length) && obj$fit_length != n)
    stop("FFT foi ajustado para comprimento diferente da serie transformada.")

  pad_len <- if (is.null(obj$pad_len)) min(max(0L, floor(n * obj$pad_pct)), n - 1L) else obj$pad_len
  padded <- fft_pad_signal(x, pad_len)
  n_pad <- length(padded)
  spectrum <- stats::fft(padded)
  half_len <- floor(n_pad / 2) + 1L
  cut <- if (is.null(obj$cutindex)) compute_cut_index(Mod(spectrum[seq_len(half_len)])^2, obj$threshold) else obj$cutindex
  cut <- max(1L, min(as.integer(cut), half_len))

  taper_len <- max(1L, floor(cut * obj$taper_pct))
  end_taper <- min(half_len, cut + taper_len)
  weights <- rep(0, half_len)
  weights[seq_len(cut)] <- 1
  if (end_taper > cut) {
    ii <- (cut + 1L):end_taper
    phase <- pi * (ii - cut) / (end_taper - cut)
    weights[ii] <- 0.5 * (1 + cos(phase))
  }

  full_weights <- numeric(n_pad)
  full_weights[seq_len(half_len)] <- weights
  neg_start <- half_len + 1L
  if (neg_start <= n_pad) {
    neg_idx <- if (n_pad %% 2L == 0L) (half_len - 1L):2L else half_len:2L
    full_weights[neg_start:n_pad] <- weights[neg_idx]
  }
  reconstructed <- Re(stats::fft(spectrum * full_weights, inverse = TRUE) / n_pad)
  result <- reconstructed[(pad_len + 1L):(pad_len + n)]
  if (length(result) != n || any(!is.finite(result))) stop("FFT produziu saida invalida.")
  names(result) <- names(original)
  result
}
