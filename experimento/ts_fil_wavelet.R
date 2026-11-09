#'@title Wavelet Filter (MODWT, borda por reflexao)
#'@description Versao local do filtro wavelet do tspredit. Remove o primeiro
#'  nivel de detalhe (W[[1]]) da MODWT e reconstrui, como a versao do pacote,
#'  mas usa `boundary = "reflection"` em vez de `"periodic"`.
#'
#'  Motivo: com `"periodic"` o vizinho "futuro" do ultimo ponto e o PRIMEIRO
#'  ponto da serie. Para o Haar o filtro equivale a
#'  x~_t = (x_{t-1} + 2 x_t + x_{t+1}) / 4, e no ultimo ponto o erro vira
#'  (x_{n-1} + x_1 - 2 x_n) / 4: um salto proporcional a diferenca entre o
#'  inicio e o fim da serie, justamente no ponto de onde o ARIMA/rede parte
#'  para prever. Com reflexao, x_{n+1} = x_n e o erro cai para (x_{n-1} - x_n)/4.
#'@param filter Filtros disponiveis: haar, d4, la8, bl14, c6. Se houver mais de
#'  um, `fit` escolhe o de maior R2 em relacao a serie original.
#'@return Um objeto `ts_fil_wavelet`.
#'@importFrom daltoolbox dal_transform
#'@importFrom daltoolbox fit
#'@importFrom daltoolbox transform
#'@export
ts_fil_wavelet <- function(filter = "haar") {
  if (!length(filter) || !all(is.character(filter)))
    stop("filter deve ser um ou mais nomes de filtro wavelet (ex.: 'haar').")
  obj <- daltoolbox::dal_transform()
  obj$filter <- filter
  class(obj) <- append("ts_fil_wavelet", class(obj))
  obj
}

wavelet_smooth_reflect <- function(x, filter) {
  n <- length(x)
  wt <- wavelets::modwt(x, filter = filter, boundary = "reflection")
  wt@W[[1]] <- matrix(0, nrow = nrow(wt@W[[1]]), ncol = ncol(wt@W[[1]]))
  y <- as.numeric(wavelets::imodwt(wt))
  # Com reflexao a serie e estendida para 2n; os n primeiros pontos sao os originais.
  if (length(y) < n) stop("Saida da MODWT menor que a entrada: ", length(y), " vs ", n)
  y[seq_len(n)]
}

#'@importFrom daltoolbox fit
#'@exportS3Method fit ts_fil_wavelet
fit.ts_fil_wavelet <- function(obj, data, ...) {
  if (length(obj$filter) > 1L) {
    x <- as.numeric(data)
    sst <- sum((x - mean(x))^2)
    best <- NULL; best_r2 <- -Inf
    for (f in obj$filter) {
      y <- wavelet_smooth_reflect(x, f)
      r2 <- if (sst > 0) 1 - sum((x - y)^2) / sst else 1
      if (is.finite(r2) && r2 > best_r2) { best <- f; best_r2 <- r2 }
    }
    obj$filter <- if (is.null(best)) obj$filter[1L] else best
  }
  obj
}

#'@importFrom daltoolbox transform
#'@exportS3Method transform ts_fil_wavelet
transform.ts_fil_wavelet <- function(obj, data, ...) {
  original <- data
  x <- as.numeric(data)
  if (length(x) < 4L) stop("Wavelet requer pelo menos quatro observacoes.")
  if (any(!is.finite(x))) stop("Wavelet requer entrada finita e sem NA.")
  filter <- if (length(obj$filter) > 1L) obj$filter[1L] else obj$filter
  result <- wavelet_smooth_reflect(x, filter)
  if (length(result) != length(x) || any(!is.finite(result)))
    stop("Falha numerica ao aplicar o filtro wavelet.")
  names(result) <- names(original)
  result
}
