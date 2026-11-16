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

## Checkpoint (retomar apos interrupcao)
- Cada (serie, cenario, filtro) concluido e salvo em `checkpoints/<output_prefix>/`.
  Se o computador desligar, rode o mesmo script de novo: o que ja foi salvo e pulado.
- Gravacao atomica (arquivo .tmp renomeado), entao desligar no meio nao corrompe nada.
- O estado do gerador aleatorio e salvo junto: a retomada da exatamente os mesmos
  resultados de uma execucao sem parada (testado com interrupcao antes e depois de emd/remd).
- `config.rds` guarda a configuracao (PRED_HORIZON, N_WINDOWS, grade, cenarios...). Se
  mudar algo, a execucao para com mensagem: use outro `output_prefix` ou apague a pasta.
- `collect_experiment("ARIMA")` junta os checkpoints e grava o Excel sem rodar nada
  (util apos rodar processos em paralelo com `series_idx`).
- O CSV do M4 e salvo em `Hourly-train.csv` no primeiro download.
- `PRED_HORIZON <- 5L`, como na versao enviada.

## LSTM
- `predict_lstm_one` chama o modelo direto (`model(X, training = FALSE)`) em vez de
  `predict()`: mesmo resultado, ~2.7x mais rapido por chamada, e evita o aviso
  "triggered tf.function retracing" (predict() monta uma tf.function por modelo).
- `fit_lstm` chama `keras3::clear_session()` antes de cada modelo para nao acumular memoria.
- `TF_CPP_MIN_LOG_LEVEL = "2"` esconde avisos de C++ do TF. As linhas de inicializacao
  ("absl::InitializeLog", "oneDNN custom operations") nao somem com isso e sao informativas.

## Figuras do artigo (`gera_graficos.R`)
Rode `Rscript gera_graficos.R` na pasta dos resultados (`ARIMA_/MLP_/ELM_industrial_sensors_results_<N>pts.xlsx`,
horizonte = N / 10) com o `IC_Common.R` e os `ts_fil_*.R`. Sai em `figuras/` (PNG 300 dpi + PDF):
1. `1_cd_<modelo>_hXX`: Diferenca Critica (Nemenyi; CD = 3,26 para k = 16, N = 50), um painel por cenario.
2. `2_heatmap_delta_rmse`: Delta RMSE% mediano vs. baseline, modelo x (cenario, h); `*` = Wilcoxon-Holm <= 0,05.
   Tabela em `figuras/tabela_delta_rmse_wilcoxon.csv`. Escala satura em +-`HEAT_LIMIT` (50%).
3. `3_pareto_rank` e `3_pareto_delta_rmse`: overhead (ms, log) x rank medio / Delta RMSE% em High_Noise,
   cor = O(1) vs. Batch (tabela `INCREMENTAL` no topo do script). Overhead de `tabela7_overhead.csv`
   (colunas Filter, Overhead_ms) se existir; senao e medido e salvo em `overhead_medido.csv`.
4. `4_diagnostico_<modelo>`: (a) Error Ratio por passo (precisa da aba `Erros_Teste`, usa o maior h),
   (b) ultimas 48 h de uma serie com ciclo de 24 h, com a % da amplitude de 24 h retida.
5. `5_dm_series_hXX`: melhor / indistinguivel / pior que o baseline pelo DM (p < 0,05 em cada cauda).
Textos ficam na lista `TXT` no topo do script (para traduzir).

### Revisao editorial (figuras em ingles, tabelas, Pareto)
- Todos os textos das figuras estao em ingles (lista `TXT`); PDF via cairo (suporta Delta, <=, >=).
  Numeros no padrao ingles (0.1; 1,000). `SHOW_TITLES <- FALSE` tira titulo/subtitulo de dentro
  das figuras (paineis do diagnostico ficam so com "(a)" e "(b)").
- Nova figura `0_contamination_profile`: ultimas 200 observacoes, 4 cenarios, subtitulo
  "Last 200 observations; dual-component contamination applied across the full series prior to windowing."
- Pareto: um painel por modelo (largura total), rotulos com ggrepel sem limite de sobreposicao e eixo
  de Delta RMSE% cortado em `PARETO_DELTA_CAP` (100%): pontos acima ficam na borda com o valor real no rotulo.
- Tabelas LaTeX em `tabelas/tab_ranks_<MODELO>.tex` (sidewaystable, \scriptsize, \tabcolsep 3.5pt,
  booktabs; precisa de \usepackage{booktabs,rotating}). Negrito = filtro melhor que o baseline por
  Wilcoxon-Holm (ou Nemenyi com `TABLE_SIG <- "nemenyi"`). Valores tambem em `tabelas/ranks_friedman.csv`.

## Revisao metodologica
Ver `RESPOSTA_REVISAO.md`: cenarios fatoriais, alvo contaminado, mais series, replicacoes,
multiplicidade, causalidade medida (`verifica_causalidade.R`, resultado em `causalidade_filtros.csv`)
e modelos sazonais (ARIMA_FOURIER, SARIMA, lags sazonais nas redes).
