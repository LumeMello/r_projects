# Experimento ARIMA / filtros - correcao do teste de Diebold-Mariano

Correcao em `IC_Common.R` (os ts_fil_*.R sao os enviados, sem mudanca):

- `dm_test_safe()`: DM com `varestimator = "bartlett"` (constante `DM_VARESTIMATOR`).
  Com o padrao "acf" e h = 3 a variancia pode ficar negativa com 30 erros; o forecast
  avisa "Variance is negative ... Proceeding with horizon h=1" e devolve um teste com h=1
  sem marcar no resultado. Bartlett e sempre >= 0, entao o teste fica com h = 3.
- Novas colunas: `DM_H_Used`, `DM_Varestimator`, `DM_Note`.
- Nova aba `Erros_Teste` (erros ponto a ponto) para refazer o DM sem rodar de novo.
- `run_experiment(..., series_idx = )` opcional para dividir as series entre processos
  (as sementes dependem do indice i, entao os resultados nao mudam).

`ts_fil_spline.R` nao foi enviado: coloque o seu nesta pasta antes de rodar.
