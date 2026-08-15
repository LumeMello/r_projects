source(url("https://raw.githubusercontent.com/cefet-rj-dal/tspredit/main/examples/seed.R"))
library(tspredit)
library(daltoolbox)
library(daltoolboxdp)
library(ggplot2)
expand_dataset <- function(x) {
  url <- attr(x, "url")
  if (is.null(url) || !nzchar(url)) x else loadfulldata(x)
}

test_values_non_filter <- list()
test_values_ma <- list()
test_values_ema <- list()
test_values_smoothing <- list()
test_values_lowess <- list()
test_values_spline <- list()
test_values_winsor <- list()
test_values_qes <- list()
test_values_ftt <- list()
test_values_wavelet <- list()
test_values_emd <- list()
test_values_remd <- list()
test_values_hp <- list()
test_values_kalman <- list()
test_values_seasonal <- list()

yvalues_non_filter <- list()
yvalues_ma <- list()
yvalues_ema <- list()
yvalues_smoothing <- list()
yvalues_lowess <- list()
yvalues_spline <- list()
yvalues_winsor <- list()
yvalues_qes <- list()
yvalues_ftt <- list()
yvalues_wavelet <- list()
yvalues_emd <- list()
yvalues_remd <- list()
yvalues_hp <- list()
yvalues_kalman <- list()
yvalues_seasonal <- list()

data(m4)
m4 <- expand_dataset(m4)
cat("Dataset: m4\n")
cat("Frequency groups:", paste(names(m4), collapse = ", "), "\n")

first_group <- names(m4)[1]
for (i in 1:10) {
  first_name <- names(m4[[first_group]])[i]
  first_series <- m4[[first_group]][[first_name]]
  
  print(any(is.na(first_series)))
  print(sum(is.na(first_series)))
  length(first_series)
  
  ts.plot(first_series, ylab = "Value", xlab = "Index", main = paste("m4", first_group, first_name))
  
  
  first_series <- na.omit(first_series)
  
  ts.plot(first_series, ylab = "Value", xlab = "Index", main = paste("m4", first_group, first_name))
  
  print(any(is.na(first_series)))
  print(sum(is.na(first_series)))
  length(first_series)
  
  set_example_seed()
  
  #Transformar essa aplicação de filtro para um for simples com if e guardar os resultados em vetores diferentes
  
  for (j in 1:14) {
    #Preciso definir o metodo de escolhas dos hiper parametros
    if(j == 1){
      ts_estructured <- ts_data(first_series,11)
    }else if(j == 2){
      filter <- ts_fil_ma(3)
      filter <- fit(filter, first_series)
      y <- transform(filter, first_series)
      y <- na.omit(y)
      ts_estructured <- ts_data(y,11)
    }else if(j == 3){
      filter <- ts_fil_ema(3)
      filter <- fit(filter, first_series)
      y <- transform(filter, first_series)
      y <- na.omit(y)
      ts_estructured <- ts_data(y,11)
    }else if(j == 4){
      filter <- ts_fil_smooth()
      filter <- fit(filter, first_series)
      y <- transform(filter, first_series)
      ts_estructured <- ts_data(y,11)
    }else if(j == 5){
      filter <- ts_fil_lowess(f = 0.2)
      filter <- fit(filter, first_series)
      y <- transform(filter, first_series)
      ts_estructured <- ts_data(y,11)
    }else if(j == 6){
      filter <- ts_fil_spline(spar = 0.5)
      filter <- fit(filter, first_series)
      y <- transform(filter, first_series)
      ts_estructured <- ts_data(y,11)
    }else if(j == 7){
      filter <- ts_fil_winsor()
      filter <- fit(filter, first_series)
      y <- transform(filter, first_series)
      ts_estructured <- ts_data(y,11)
    }else if(j == 8){
      filter <- ts_fil_qes(gamma = FALSE)
      filter <- fit(filter, first_series)
      y <- transform(filter, first_series)
      y <- na.omit(y)
      ts_estructured <- ts_data(y,11)
    }else if(j == 9){
      filter <- ts_fil_fft()
      filter <- fit(filter, first_series)
      y <- transform(filter, first_series)
      ts_estructured <- ts_data(y,11)
    }else if(j == 10){
      filter <- ts_fil_wavelet()
      filter <- fit(filter, first_series)
      y <- transform(filter, first_series)
      ts_estructured <- ts_data(y,11)
    }else if(j == 11){
      filter <- ts_fil_emd()
      filter <- fit(filter, first_series)
      y <- transform(filter, first_series)
      ts_estructured <- ts_data(y,11)
    }else if(j == 12){
      filter <- ts_fil_remd()
      filter <- fit(filter, first_series)
      y <- transform(filter, first_series)
      ts_estructured <- ts_data(y,11)
    }else if(j == 13){
      filter <- ts_fil_hp()
      filter <- fit(filter, first_series)
      y <- transform(filter, first_series)
      ts_estructured <- ts_data(y,11)
    }else if(j == 14){
      filter <- ts_fil_kalman(H = 0.1, Q = 1)
      filter <- fit(filter, first_series)
      y <- transform(filter, first_series)
      ts_estructured <- ts_data(y,11)
    }else{
      #Precisa resolver o problema de descobrir qual a frequencia da serie temporal
    }
    
    #talvez ao inves de ser 10 passsos pensar em um modelo que calcule x ciclos
    samp <- ts_sample(ts_estructured, test_size = 10)
    io_train <- ts_projection(samp$train)
    io_test <- ts_projection(samp$test)
    
    # MLP
    preproc <- ts_norm_gminmax()
    
    
    model <- ts_mlp(ts_norm_gminmax(), input_size=10, size=10, decay=0)
    model <- fit(model, x=io_train$input, y=io_train$output)
    
    adjust <- predict(model, io_train$input)
    adjust <- as.vector(adjust)
    output <- as.vector(io_train$output)
    ev_adjust <- evaluate(model, output, adjust)
    ev_adjust$mse
    
    prediction <- predict(model, x=io_test$input[1,], steps_ahead=10)
    prediction <- as.vector(prediction)
    output <- as.vector(io_test$output)
    
    
    
    ev_test <- evaluate(model, output, prediction)
    ev_test
    
    #guardando os valores finais para comparação se necessário(não esquecer de fazer um if para colocar na list do filtro certo)
    yvalues <- c(io_train$output, io_test$output)
    
    if(j == 1){
      yvalues_non_filter <- append(yvalues_non_filter, list(yvalues))
      test_values_non_filter <- append(test_values_non_filter, list(ev_test$metrics))
    }else if(j == 2){
      yvalues_ma <- append(yvalues_ma, list(yvalues))
      test_values_ma <- append(test_values_ma, list(ev_test$metrics))
    }else if(j == 3){
      yvalues_ema <- append(yvalues_ema, list(yvalues))
      test_values_ema <- append(test_values_ema, list(ev_test$metrics))
    }else if(j == 4){
      yvalues_smoothing <- append(yvalues_smoothing, list(yvalues))
      test_values_smoothing <- append(test_values_smoothing, list(ev_test$metrics))
    }else if(j == 5){
      yvalues_lowess <- append(yvalues_lowess, list(yvalues))
      test_values_lowess <- append(test_values_lowess, list(ev_test$metrics))
    }else if(j == 6){
      yvalues_spline <- append(yvalues_spline, list(yvalues))
      test_values_spline <- append(test_values_spline, list(ev_test$metrics))
    }else if(j == 7){
      yvalues_winsor <- append(yvalues_winsor, list(yvalues))
      test_values_winsor <- append(test_values_winsor, list(ev_test$metrics))
    }else if(j == 8){
      yvalues_qes <- append(yvalues_qes, list(yvalues))
      test_values_qes <- append(test_values_qes, list(ev_test$metrics))
    }else if(j == 9){
      yvalues_ftt <- append(yvalues_ftt, list(yvalues))
      test_values_ftt <- append(test_values_ftt, list(ev_test$metrics))
    }else if(j == 10){
      yvalues_wavelet <- append(yvalues_wavelet, list(yvalues))
      test_values_wavelet <- append(test_values_wavelet, list(ev_test$metrics))
    }else if(j == 11){
      yvalues_emd <- append(yvalues_emd, list(yvalues))
      test_values_emd <- append(test_values_emd, list(ev_test$metrics))
    }else if(j == 12){
      yvalues_remd <- append(yvalues_remd, list(yvalues))
      test_values_remd <- append(test_values_remd, list(ev_test$metrics))
    }else if(j == 13){
      yvalues_hp <- append(yvalues_hp, list(yvalues))
      test_values_hp <- append(test_values_hp, list(ev_test$metrics))
    }else if(j == 14){
      yvalues_kalman <- append(yvalues_kalman, list(yvalues))
      test_values_kalman <- append(test_values_kalman, list(ev_test$metrics))
    }else{
      yvalues_seasonal <- append(yvalues_seasonal, list(yvalues))
      test_values_seasonal <- append(test_values_seasonal, list(ev_test$metrics))
    }
    
    
    plot_ts_pred(y=yvalues, yadj=adjust, ypre=prediction, color_prediction="orange") + theme(text = element_text(size=16))
  }
  
}
