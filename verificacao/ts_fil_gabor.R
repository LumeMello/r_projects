#'@title One dimensional Gabor Smoothing Filter
#'@description Filtro de suavizacao 1D por kernel de Gabor simetrico.
#'@param window_size Desvio padrao do envelope gaussiano; deve ser positivo.
#'@param central_frequency Frequencia central em ciclos por amostra (0, 0.5].
#'@param max_gain Ganho L1 maximo tolerado do kernel normalizado (padrao 3);
#'  acima disso o filtro recusa os parametros em vez de amplificar o sinal.
#'@return Um objeto `ts_fil_gabor`.
#'@importFrom daltoolbox dal_transform
#'@importFrom daltoolbox fit
#'@importFrom daltoolbox transform
#'@export
ts_fil_gabor <- function(window_size = 3, central_frequency = 0.1, max_gain = 3) {
  if (length(max_gain) != 1L || !is.finite(max_gain) || max_gain < 1)
    stop("max_gain deve ser um numero >= 1.")
  if (length(window_size) != 1L || !is.finite(window_size) || window_size <= 0)
    stop("window_size deve ser um numero positivo.")
  if (length(central_frequency) != 1L || !is.finite(central_frequency) ||
      central_frequency <= 0 || central_frequency > 0.5)
    stop("central_frequency deve estar no intervalo (0, 0.5].")
  obj <- daltoolbox::dal_transform()
  obj$ws <- window_size
  obj$cf <- central_frequency
  obj$max_gain <- max_gain
  class(obj) <- append("ts_fil_gabor", class(obj))
  obj
}

#'@importFrom daltoolbox fit
#'@exportS3Method fit ts_fil_gabor
fit.ts_fil_gabor <- function(obj, data, ...) obj

#'@importFrom daltoolbox transform
#'@exportS3Method transform ts_fil_gabor
transform.ts_fil_gabor <- function(obj, data, ...) {
  original <- data
  x <- as.numeric(data)
  n <- length(x)
  if (n < 4L) stop("Gabor requer pelo menos quatro observacoes.")
  if (any(!is.finite(x))) stop("Gabor requer entrada finita e sem NA.")

  radius <- min(n - 1L, max(1L, ceiling(3 * obj$ws)))
  k <- -radius:radius
  envelope <- exp(-(k^2) / (2 * obj$ws^2))
  kernel <- envelope * cos(2 * pi * obj$cf * k)
  kernel_sum <- sum(kernel)
  if (!is.finite(kernel_sum) || abs(kernel_sum) < 1e-10)
    stop("Kernel Gabor com soma proxima de zero; ajuste window_size/central_frequency.")
  kernel <- kernel / kernel_sum

  # Guarda de ganho: a normalizacao por kernel_sum so e valida se o kernel nao
  # amplifica o sinal. Com central_frequency proxima de 0.5 o kernel_sum fica
  # ~0 (ex.: -0.0059 para freq=2) e a divisao infla os coeficientes (norma L1
  # ~1280 para freq=2, ~1110 para freq=3, ~520 para freq=4, com window_size=3).
  # A saida deixa de ser suavizada e passa a explodir (RMSE de milhoes).
  gain <- sum(abs(kernel))
  if (!is.finite(gain) || gain > obj$max_gain)
    stop("Kernel Gabor amplifica o sinal (ganho L1 = ", signif(gain, 3),
         " > ", obj$max_gain, "); use periodo (1/central_frequency) maior ou window_size menor.")

  # Reflexao nas bordas; sides=2 torna o filtro bilateral e nao causal.
  left <- rev(x[2:(radius + 1L)])
  right <- rev(x[(n - radius):(n - 1L)])
  padded <- c(left, x, right)
  filtered <- as.numeric(stats::filter(padded, kernel, sides = 2))
  result <- filtered[(radius + 1L):(radius + n)]
  if (length(result) != n || any(!is.finite(result)))
    stop("Falha numerica ao aplicar o filtro Gabor.")
  names(result) <- names(original)
  result
}
