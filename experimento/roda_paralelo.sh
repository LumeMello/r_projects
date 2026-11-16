#!/usr/bin/env bash
# Roda um experimento em N processos paralelos e junta o resultado no fim.
# Uso: bash roda_paralelo.sh <MODELO> <H> [VARIANTES] [N_PROCESSOS]
#   ex.: bash roda_paralelo.sh ARIMA_FOURIER 5 base 4
#        bash roda_paralelo.sh MLP 3 S24+FACT 6
# Logs em logs/<MODELO>_<VARIANTES>_h<H>_parte<k>.log. Se cair no meio, rode de novo:
# o checkpoint retoma de onde parou.
set -u
MODELO="$1"; H="$2"; VAR="${3:-base}"; N="${4:-4}"
mkdir -p logs
TAG="${MODELO}_${VAR//+/-}_h${H}"
pids=()
for k in $(seq 1 "$N"); do
  Rscript roda_experimento.R "$MODELO" "$H" "$VAR" "$k/$N" > "logs/${TAG}_parte${k}.log" 2>&1 &
  pids+=($!)
  echo "parte $k/$N iniciada (pid $!), log: logs/${TAG}_parte${k}.log"
done
falhas=0
for p in "${pids[@]}"; do wait "$p" || falhas=$((falhas + 1)); done
if [ "$falhas" -gt 0 ]; then
  echo "$falhas parte(s) terminaram com erro; veja os logs e rode este comando de novo (retoma do checkpoint)."
  exit 1
fi
Rscript roda_experimento.R "$MODELO" "$H" "$VAR" junta > "logs/${TAG}_junta.log" 2>&1 && \
  tail -n 1 "logs/${TAG}_junta.log"
