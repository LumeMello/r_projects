#'@title One dimensional Gabor Smoothing Filter
#'@description Filtro de suavização 1D baseado em kernel de Gabor integrado ao daltoolbox.
#'
#'@param window_size Tamanho da janela / desvio padrão gaussiano (omega).
#'@param central_frequency Frequência central da onda portadora.
#'@return A `ts_fil_gabor` object.
#'
#'@importFrom daltoolbox dal_transform
#'@importFrom daltoolbox fit
#'@importFrom daltoolbox transform
#'@export
ts_fil_gabor <- function(window_size = 3, central_frequency = 0.1) {
  obj <- dal_transform()
  obj$ws <- window_size
  obj$cf <- central_frequency
  class(obj) <- append("ts_fil_gabor", class(obj))
  return(obj)
}

#'@importFrom daltoolbox fit
#'@exportS3Method fit ts_fil_gabor
fit.ts_fil_gabor <- function(obj, data, ...) {
  return(obj)
}

#'@importFrom daltoolbox transform
#'@importFrom stats filter
#'@exportS3Method transform ts_fil_gabor
transform.ts_fil_gabor <- function(obj, data, ...) {
  data_num <- as.numeric(data)
  n <- length(data_num)
  
  if (n < 4 || sum(!is.na(data_num)) < 4) {
    return(data)
  }
  
  k_radius <- max(1, ceiling(3 * obj$ws))
  k_seq <- -k_radius:k_radius
  
  gauss_env <- (1 / (sqrt(2 * pi) * obj$ws)) * exp(- (k_seq^2) / (2 * obj$ws^2))
  gabor_kernel <- gauss_env * cos(2 * pi * obj$cf * k_seq)
  
  if (sum(gabor_kernel) != 0) {
    gabor_kernel <- gabor_kernel / sum(gabor_kernel)
  }
  
  prefix <- rev(data_num[2:(k_radius + 1)])
  suffix <- rev(data_num[n - k_radius: (n-1)])
  padded_data <- c(prefix, data_num, suffix)
  
  filtered_padded <- stats::filter(padded_data, gabor_kernel, sides = 2)
  
  filtered <- as.numeric(filtered_padded[(k_radius +1):(k_radius+n)])
  
  na_mask <- is.na(filtered)
  if (any(na_mask)) {
    filtered[na_mask] <- data_num[na_mask]
  }
  
  return(filtered)
}
