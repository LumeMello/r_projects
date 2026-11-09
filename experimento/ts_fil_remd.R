#'@title Robust EMD Filter
#'@description CEEMD denoising por remocao dos primeiros IMFs de maior frequencia.
#'@param noise Amplitude do ruido auxiliar do CEEMD (positivo).
#'@param trials Numero de realizacoes do ensemble (inteiro positivo).
#'@param max_imfs Limite superior de IMFs removidos; NULL usa o corte automatico.
#'@return Um objeto `ts_fil_remd`.
#'@references Z. Wu and N. E. Huang (2009). Ensemble Empirical Mode Decomposition.
#'@importFrom daltoolbox dal_transform
#'@importFrom daltoolbox fit
#'@importFrom daltoolbox transform
#'@export
ts_fil_remd <- function(noise = 0.1, trials = 5, max_imfs = NULL) {
  if (length(noise) != 1L || !is.finite(noise) || noise <= 0) stop("noise deve ser positivo.")
  if (length(trials) != 1L || !is.finite(trials) || trials < 1 || trials != as.integer(trials))
    stop("trials deve ser um inteiro positivo.")
  if (!is.null(max_imfs) && (length(max_imfs) != 1L || !is.finite(max_imfs) ||
      max_imfs < 0 || max_imfs != as.integer(max_imfs)))
    stop("max_imfs deve ser NULL ou um inteiro nao negativo.")
  obj <- daltoolbox::dal_transform()
  obj$noise <- noise
  obj$trials <- as.integer(trials)
  obj$max_imfs <- if (is.null(max_imfs)) NULL else as.integer(max_imfs)
  class(obj) <- append("ts_fil_remd", class(obj))
  obj
}

fc_roughness <- function(x) {
  x <- as.numeric(x)
  if (any(!is.finite(x))) stop("IMF contem valores nao finitos.")
  if (length(x) <= 2L) return(0)
  d <- diff(x)
  s <- stats::sd(d)
  if (!is.finite(s) || s <= .Machine$double.eps) return(0)
  mean((diff((d - mean(d)) / s)^2) / 4)
}

#'@importFrom daltoolbox fit
#'@exportS3Method fit ts_fil_remd
fit.ts_fil_remd <- function(obj, data, ...) obj

#'@importFrom daltoolbox transform
#'@exportS3Method transform ts_fil_remd
transform.ts_fil_remd <- function(obj, data, ...) {
  x <- as.numeric(data)
  n <- length(x)
  if (n < 4L) stop("REMD requer pelo menos quatro observacoes.")
  if (any(!is.finite(x))) stop("REMD/CEEMD requer entrada finita e sem NA.")

  decomposition <- tryCatch(
    suppressWarnings(hht::CEEMD(sig = x, tt = seq_len(n), noise.amp = obj$noise,
                                trials = obj$trials, verbose = FALSE)),
    error = function(e) stop("Falha no CEEMD: ", conditionMessage(e))
  )
  if (is.null(decomposition$imf) || is.null(decomposition$nimf))
    stop("CEEMD nao retornou a matriz de IMFs esperada.")
  imfs <- as.matrix(decomposition$imf)
  nimf <- as.integer(decomposition$nimf)
  if (nrow(imfs) != n && ncol(imfs) == n) imfs <- t(imfs)
  if (nrow(imfs) != n || ncol(imfs) < 1L || nimf < 1L)
    stop("Dimensoes inesperadas na decomposicao CEEMD.")
  nimf <- min(nimf, ncol(imfs))
  if (any(!is.finite(imfs[, seq_len(nimf), drop = FALSE])))
    stop("CEEMD produziu IMFs com valores nao finitos.")

  rough <- vapply(seq_len(nimf), function(k) fc_roughness(imfs[, k]), numeric(1))
  cumulative <- cumsum(rough)
  cutoff <- 0L
  if (nimf >= 2L && sum(cumulative) > .Machine$double.eps) {
    curvature <- tryCatch(daltoolbox::transform(daltoolbox::fit_curvature_min(), cumulative),
                          error = function(e) NULL)
    if (!is.null(curvature$x) && length(curvature$x) == 1L && is.finite(curvature$x))
      cutoff <- as.integer(curvature$x)
    else cutoff <- 1L
  } else if (nimf >= 2L) {
    cutoff <- 1L
  }
  if (!is.null(obj$max_imfs)) cutoff <- min(cutoff, obj$max_imfs)
  cutoff <- max(0L, min(cutoff, nimf - 1L))

  result <- x
  if (cutoff > 0L) result <- x - rowSums(imfs[, seq_len(cutoff), drop = FALSE])
  if (length(result) != n || any(!is.finite(result))) stop("REMD produziu saida invalida.")
  names(result) <- names(data)
  result
}
