# Verificacao das correcoes (wavelet, gabor, freq)

Arquivos exatamente como foram rodados com `Rscript verifica_correcoes.R` (R 4.3.3).
Saida completa em `saida.txt`.

- `ts_fil_gabor.R` = `ts_fil_gabor_2.R` enviado (renomeado porque IC_Common_2.R faz `source("ts_fil_gabor.R")`).
- `ts_fil_fft.R`, `ts_fil_remd.R`, `ts_fil_spline.R` sao **placeholders vazios** (os originais nao foram enviados);
  com eles valem as versoes do pacote tspredit. Substitua pelas suas versoes antes de rodar o experimento completo.
- Pacotes: forecast 8.21.1, writexl 1.5.0, tspredit 2.0.717, daltoolbox 1.3.787, daltoolboxdp 1.3.767,
  PMCMRplus 1.9.12, wavelets 0.3.0.2. Com R < 4.5 foi preciso instalar do arquivo do CRAN
  arules 1.7-10 e arulesSequences 0.2-31 (as versoes atuais exigem R >= 4.5).
