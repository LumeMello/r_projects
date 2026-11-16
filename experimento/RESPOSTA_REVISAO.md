# Revisão — o que foi implementado e como rodar

Todas as opções novas ficam no topo do `IC_Common.R`. Com os valores padrão, o experimento é
**idêntico ao original** (verificado: contaminação, montagem das entradas e previsão recursiva
reproduzem a versão anterior, e os checkpoints antigos continuam válidos).

Para cada variante, defina as opções **depois** de `source("IC_Common.R")` e use um
`output_prefix` próprio (o checkpoint recusa misturar configurações diferentes):

```r
source("IC_Common.R")
SCENARIO_SET <- "factorial"
run_experiment("ARIMA", data.frame(dummy = 1), "ARIMA_FACT")
```

## Scripts de execução

`roda_experimento.R` roda os modelos ARIMA, ARIMA_FOURIER, SARIMA, MLP e ELM (LSTM desativado por prazo) e
qualquer variante; `roda_paralelo.sh` divide as séries em N processos e junta no fim.
**Confira a grade do MLP** no topo do `roda_experimento.R` (copie a do seu MLP.R).

```bash
# teste de 1-2 min antes de cada rodada longa
IC_TESTE_RAPIDO=1 Rscript roda_experimento.R MLP 5 S24

# rodadas (h = 3, 5, 10), em 4 processos
bash roda_paralelo.sh ARIMA 5 FACT 4           # original + cenarios fatoriais
bash roda_paralelo.sh ARIMA_FOURIER 5 FACT 4   # modelo sazonal
bash roda_paralelo.sh MLP 5 FACT 4
bash roda_paralelo.sh MLP 5 CT+FACT 4          # alvo contaminado  -> MLP_CT
bash roda_paralelo.sh MLP 5 S24+FACT 4         # lag de 24 h       -> MLP_S24
bash roda_paralelo.sh ELM 5 REP3 4             # 3 realizacoes do ruido
```

No RStudio: edite `MODELO`, `H`, `VARIANTES` e `PARTE` no topo do script e rode com Source.
Como FACT, REPx e Nx incluem os resultados originais, rodar já com `FACT` substitui a
rodada base (os 4 cenários originais saem iguais); use `PLOT_SCENARIOS` no `gera_graficos.R`
para mostrar só os 4 originais nas figuras principais.

## Custo (referência medida: ARIMA, 50 séries, 4 cenários, h = 5 ≈ 6,4 h em 1 núcleo)

| Mudança | Multiplica o tempo por | Observação |
|---|---|---|
| `SCENARIO_SET <- "factorial"` | 2,5 (10 cenários em vez de 4) | os 4 originais podem ser reaproveitados se rodar só os 6 novos em outra pasta |
| `N_REPLICATIONS <- R` | R | a réplica 1 é o experimento original (reaproveita checkpoint) |
| `N_SERIES <- N` | N / 50 | com `SERIES_SELECTION = "first"` as 50 já feitas são reaproveitadas |
| `TARGET_MODE <- "corrupted"` | 1 (rodada nova só para MLP/ELM) | ARIMA já não usa alvo limpo no treino |
| `ARIMA_FOURIER` | ~1,5–2 | ARIMA com termos de Fourier de 24 h |
| `SARIMA` | bem maior (lento) | use em um subconjunto, se o revisor insistir |
| `NN_SEASONAL_LAGS <- c(24L)` | ~1 | só acrescenta uma entrada às redes |

Paralelize com `series_idx` (um processo por núcleo, mesmo `output_prefix`) e junte com
`collect_experiment("<prefixo>")`.

## Item a item

**1. Gaussiano e outliers não isolados (Alta)** — `SCENARIO_SET <- "factorial"` acrescenta
`Gauss_Mod/High/Sev` (só ruído) e `Out_Mod/High/Sev` (só spikes), com a mesma intensidade dos
cenários combinados e os **mesmos sorteios** (common random numbers): para cada nível,
combinado = só-ruído + só-spikes − limpa, exatamente. Assim a diferença entre cenários
isola o efeito de cada componente.

**2. Condições de treino diferentes entre modelos (Alta)** — no original, as redes treinam
com alvo **limpo** (oráculo) e o ARIMA só vê a série contaminada. `TARGET_MODE <- "corrupted"`
faz as redes treinarem e validarem com a série contaminada (como o ARIMA); a avaliação final
continua contra a série limpa. Sugestão: rodar MLP e ELM com `output_prefix = "MLP_CT"`,
`"ELM_CT"` e apresentar como experimento de robustez.

**3. Só 50 séries (Alta)** — o M4 Hourly tem 414 séries. `N_SERIES <- 414L` (ou 100/200)
mantém as 50 primeiras já calculadas. Alternativa: `SERIES_SELECTION <- "random"` sorteia
uma amostra (pasta nova). A CD do Nemenyi cai com N (≈ 1,13 para N = 414).

**4. Uma só realização do ruído (Média)** — `N_REPLICATIONS <- 5L` (ou 3). Cada série recebe
R contaminações independentes; ranks, Friedman/Nemenyi e Wilcoxon usam o RMSE médio por série
(as séries continuam sendo os blocos), e o DM usa cada realização. Colunas `Replication` nos
resultados e na aba `Erros_Teste`.

**5. Multiplicidade (Média)** — `gera_graficos.R` agora calcula três correções e grava todas em
`figuras/tabela_delta_rmse_wilcoxon.csv`: Holm por família (15 filtros em cada
modelo × cenário × horizonte, a do artigo), Holm global e Benjamini-Hochberg global.
`MULTIPLICITY` escolhe a usada nos asteriscos e no negrito; `figuras/resumo_multiplicidade.csv`
mostra quantos resultados mudam. `DM_MULTIPLICITY <- "bh"` corrige o DM entre séries.

**6. Filtros causais × não causais (Média)** — medido com `verifica_causalidade.R` (o filtro é
causal se a saída no prefixo não muda quando chegam novos dados). Causais: **MA, EMA, Kalman**.
Não causais: os outros 12 (Winsor e QES revisam pouco: o último ponto muda < 0,02 desvio padrão;
LOESS, Spline, Wavelet, EMD/REMD e Gabor mudam o último ponto em 0,1–0,5 desvio padrão).
`filter_causality` no IC_Common foi atualizado; as figuras e tabelas separam os dois grupos
(`SPLIT_CAUSALITY`), e há diagramas de CD por grupo (`1_cd_*_causal`, `1_cd_*_noncausal`).

**7. Sazonalidade de 24 h (Média)** — novos modelos `ARIMA_FOURIER` (regressão harmônica
dinâmica: Fourier de 24 h, `FOURIER_K = 4`, + erros ARIMA) e `SARIMA` (período 24). Para as
redes, `NN_SEASONAL_LAGS <- c(24L)` (ou `c(24L, 168L)`) acrescenta x[t−24] às entradas.

**8. Generalização excessiva (Baixa)** — texto. Sugestão abaixo.

## Textos sugeridos (inglês)

*Multiple comparisons.* "Pairwise comparisons against the unfiltered baseline were performed
with the paired Wilcoxon signed-rank test, with Holm correction applied within each family of
15 filters (one family per model, contamination scenario and horizon). To assess sensitivity
to the definition of the family, we additionally report Holm and Benjamini–Hochberg corrections
across all N tests (Supplementary Table X); the qualitative conclusions are unchanged/the number
of significant improvements changes from A to B. Friedman tests were followed by the Nemenyi
post-hoc procedure, which already controls the family-wise error across all pairs of methods."

*Causality.* "All filters are applied only to the history available at each forecast origin,
so no test information leaks into training. However, we empirically verified (Section X) that
only MA, EMA and Kalman are causal; the remaining filters revise past values—including the
forecast origin—when new observations arrive. Results are therefore reported separately for the
two groups, and the gains of non-causal filters should be interpreted as an upper bound for
online deployment."

*Scope.* "Our conclusions are restricted to hourly series from the M4 Hourly subset with
synthetic, independently injected Gaussian noise and isolated spikes, short horizons
(h ≤ 10) and the specific model configurations evaluated. They do not cover other sampling
frequencies, real sensor failure modes (drift, stuck values, missing blocks, correlated
noise), multivariate settings, or models with explicit long seasonal memory beyond those
tested in Section X."
