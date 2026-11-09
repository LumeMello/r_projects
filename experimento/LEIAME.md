# Experimento ARIMA / filtros - correcao do teste de Diebold-Mariano

## Problema
`forecast::dm.test(..., h = 3)` usa por padrao o estimador "acf" da variancia, que pode
dar negativo com so 30 erros. O forecast entao avisa "Variance is negative ... Proceeding
with horizon h=1" e devolve um teste com h=1, sem marcar isso no resultado.

## Conferencia (ARIMA, 50 series x 4 cenarios x 15 filtros = 3000 testes)
- Rodei o experimento de novo salvando os erros ponto a ponto e refiz o DM.
- So 2 linhas tinham variancia negativa (os 2 avisos): H21/High_Noise/kalman e
  H44/Sev_Noise/winsor. Nenhuma das duas e significativa nem antes nem depois.
- As outras 2998 linhas reproduzem exatamente o DM da planilha original.
- Usar "bartlett" em TODAS as linhas mudaria a conclusao (p <= 0.05) de 52 linhas,
  porque Bartlett encolhe a variancia e deixa o teste menos conservador. Por isso a
  correcao usa Bartlett so como fallback.

## Correcao em `IC_Common.R` (ts_fil_*.R sao os enviados, sem mudanca)
- `dm_test_safe()`: "acf" (padrao do teste); se a variancia der negativa, refaz com
  "bartlett" mantendo h = 3. Constantes `DM_VARESTIMATOR` e `DM_FALLBACK_VARESTIMATOR`.
- Novas colunas: `DM_H_Used`, `DM_Varestimator`, `DM_Note`.
- Nova aba `Erros_Teste` (erros ponto a ponto) para refazer o DM sem rodar de novo.
- `run_experiment(..., series_idx = )` opcional para dividir as series entre processos
  (as sementes dependem do indice i, entao os resultados nao mudam).

## `conferencia_dm/`
- `ARIMA_industrial_sensors_results_30pts_DM_corrigido.xlsx`: sua planilha com as 2 linhas
  corrigidas e as colunas novas (o resto identico; o resumo estatistico nao usa o DM).
- `comparacao_dm.csv`: DM acf x bartlett para as 3000 linhas.
- `confere_dm.R`, `confere2.R`: scripts usados (caminhos da sessao de verificacao).

Obs.: `ts_fil_spline.R` nao foi enviado; usei o do tspredit (deu os mesmos RMSE).
Na rodada nova, 100 linhas de ema/seasonal tiveram RMSE com diferencas numericas minimas
(mesmos parametros; provavel versao de R/BLAS), sem efeito no DM.
