suppressPackageStartupMessages(library(forecast))
fs <- sprintf("ARIMA_w%d_industrial_sensors_results_30pts.xlsx", 1:4)
res <- do.call(rbind, lapply(fs, readxl::read_excel, sheet = "Resultados_Gerais"))
err <- do.call(rbind, lapply(fs, readxl::read_excel, sheet = "Erros_Teste"))
orig <- readxl::read_excel(list.files("/root/.claude/uploads/ad9e2174-8dcd-5c47-bcd8-a7fc94afc86a", pattern = "xlsx$", full.names = TRUE), 1)
key <- function(d) paste(d$Series_ID, d$Scenario, d$Filter)
E <- split(err$Error, key(err))
run_dm <- function(e1, e2, ve) {
  w <- FALSE
  r <- tryCatch(withCallingHandlers(dm.test(e1, e2, alternative = "greater", h = 3, power = 2, varestimator = ve),
         warning = function(x) { w <<- TRUE; invokeRestart("muffleWarning") }), error = function(e) NULL)
  if (is.null(r)) return(c(stat = NA, p = NA, h = NA, warn = w))
  c(stat = unname(r$statistic), p = r$p.value, h = unname(r$parameter[1]), warn = w)
}
rows <- res[res$Filter != "non_filter" & res$Status == "CONVERGED", ]
out <- t(sapply(seq_len(nrow(rows)), function(j) {
  r <- rows[j, ]; e1 <- E[[paste(r$Series_ID, r$Scenario, "non_filter")]]; e2 <- E[[key(r)]]
  if (isTRUE(all.equal(e1, e2, tolerance = 1e-9))) return(c(acf_stat=0, acf_p=1, acf_h=3, acf_warn=0, bt_stat=0, bt_p=1, bt_h=3, bt_warn=0))
  a <- run_dm(e1, e2, "acf"); b <- run_dm(e1, e2, "bartlett")
  setNames(c(a, b), c("acf_stat","acf_p","acf_h","acf_warn","bt_stat","bt_p","bt_h","bt_warn"))
}))
d <- cbind(as.data.frame(rows[, c("Series_ID","Scenario","Filter","DM_Stat","DM_Pvalue","DM_H_Used","DM_Note")]), out)
m <- merge(d, orig[, c("Series_ID","Scenario","Filter","RMSE","DM_Stat","DM_Pvalue")], by = c("Series_ID","Scenario","Filter"), suffixes = c("", "_orig"))
m$RMSE_new <- rows$RMSE[match(key(m), key(rows))]
cat("Linhas comparadas:", nrow(m), "\n")
cat("RMSE igual a planilha original (tol 1e-6):", sum(abs(m$RMSE_new - m$RMSE) < 1e-6), "\n")
cat("Filtros com RMSE diferente:", paste(unique(m$Filter[abs(m$RMSE_new - m$RMSE) >= 1e-6]), collapse = ", "), "\n")
same <- abs(m$RMSE_new - m$RMSE) < 1e-6
cat("DM 'acf' reproduz o DM da planilha original (linhas com RMSE igual):",
    sum(abs(m$acf_stat - m$DM_Stat_orig)[same] < 1e-6, na.rm = TRUE), "de", sum(same), "\n\n")
cat("== Estimador antigo (acf): linhas com aviso 'Variance is negative' -> caiu para h=1 ==\n")
print(table(Cenario = m$Scenario, Aviso = m$acf_warn == 1))
cat("\nPor filtro:\n"); print(sort(table(m$Filter[m$acf_warn == 1]), decreasing = TRUE))
cat("\n== Estimador novo (bartlett): avisos:", sum(m$bt_warn == 1, na.rm = TRUE), "| h usado:", paste(unique(m$bt_h), collapse = ","), "\n")
cat("Codigo corrigido (dm_test_safe) = calculo manual bartlett:", sum(abs(m$DM_Stat - m$bt_stat) < 1e-9, na.rm = TRUE), "de", nrow(m), "\n\n")
cat("== Significativos (p <= 0.05) por cenario ==\n")
print(rbind(acf_antigo = tapply(m$acf_p <= 0.05, m$Scenario, sum, na.rm = TRUE),
            bartlett = tapply(m$bt_p <= 0.05, m$Scenario, sum, na.rm = TRUE)))
w <- m[m$acf_warn == 1, ]
cat("\nNas linhas com aviso: significativos com acf(h=1) =", sum(w$acf_p <= 0.05), "| com bartlett(h=3) =", sum(w$bt_p <= 0.05), "\n")
cat("Mudaram de conclusao (p<=0.05):", sum((m$acf_p <= 0.05) != (m$bt_p <= 0.05), na.rm = TRUE), "\n")
print(head(w[order(w$acf_p), c("Series_ID","Scenario","Filter","acf_stat","acf_p","bt_stat","bt_p")], 15), row.names = FALSE)
write.csv(m, "comparacao_dm.csv", row.names = FALSE)
