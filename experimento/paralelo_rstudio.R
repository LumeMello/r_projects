# =============================================================================
# EXECUCAO PARALELA PELO RSTUDIO SERVER (Linux)
# Uso no console do RStudio, na pasta do experimento:
#   source("paralelo_rstudio.R")
#   lancar("ELM", 5, "FACT")              # divide as series em n processos e junta no fim
#   status("ELM", 5, "FACT")              # progresso (pode chamar quantas vezes quiser)
#   parar()                               # interrompe tudo (o checkpoint guarda o progresso)
#
#   # fila: varios experimentos, um depois do outro, sem precisar ficar olhando
#   lancar_fila(list(c("ARIMA", 3, "FACT"), c("ARIMA", 5, "FACT"), c("ARIMA", 10, "FACT"),
#                    c("ELM", 3, "FACT"), c("MLP", 3, "S24+FACT")))
#   status_fila()
#
# Os processos rodam fora da sessao do RStudio (setsid + nohup): fechar o navegador,
# reiniciar o R ou a sessao ser suspensa NAO interrompe a rodada. Se o servidor
# reiniciar, rode lancar()/lancar_fila() de novo: o checkpoint retoma de onde parou.
# =============================================================================

# Processos: por padrao, todos os nucleos menos 1 (deixa folga para o RStudio).
n_padrao <- function() max(1L, parallel::detectCores(logical = FALSE) - 1L)

# setsid desliga o processo da sessao do RStudio; se nao existir, nohup ja basta na maioria dos casos.
.desliga <- function() if (nzchar(Sys.which("setsid"))) "setsid nohup" else "nohup"

.tag <- function(modelo, h, variantes) sprintf("%s_%s_h%s", modelo, gsub("+", "-", variantes, fixed = TRUE), h)

.rodando <- function(padrao = "roda_experimento.R|roda_paralelo.sh") {
  out <- suppressWarnings(system2("pgrep", c("-fa", shQuote(padrao)), stdout = TRUE, stderr = FALSE))
  out[!grepl("pgrep", out)]
}

lancar <- function(modelo, h, variantes = "base", n = n_padrao()) {
  stopifnot(file.exists("roda_experimento.R"), file.exists("roda_paralelo.sh"))
  tag <- .tag(modelo, h, variantes)
  ja <- .rodando(sprintf("roda_paralelo.sh %s %s %s ", modelo, h, variantes))
  if (length(ja)) { message("Ja esta rodando: ", tag, ". Use status() ou parar()."); return(invisible(FALSE)) }
  dir.create("logs", showWarnings = FALSE)
  unlink(file.path("logs", paste0(tag, "_junta.log")))   # evita mostrar a juncao de uma rodada antiga
  log <- file.path("logs", paste0(tag, "_lancador.log"))
  cmd <- sprintf("%s bash roda_paralelo.sh %s %s %s %d > %s 2>&1 < /dev/null",
                 .desliga(), shQuote(modelo), shQuote(h), shQuote(variantes), as.integer(n), shQuote(log))
  system(cmd, wait = FALSE)
  message(sprintf("Lancado %s em %d processos. Acompanhe com status(\"%s\", %s, \"%s\").", tag, n, modelo, h, variantes))
  invisible(TRUE)
}

# Total esperado de combinacoes serie x cenario x replicacao x filtro para a variante.
.esperado <- function(variantes) {
  v <- toupper(strsplit(variantes, "+", fixed = TRUE)[[1]])
  series <- if ("N414" %in% v) 414L else if ("N200" %in% v) 200L else if ("N100" %in% v) 100L else 50L
  cen <- if ("FACT" %in% v) 10L else 4L
  rep <- if ("REP5" %in% v) 5L else if ("REP3" %in% v) 3L else 1L
  if (nzchar(Sys.getenv("IC_TESTE_RAPIDO"))) return(2L * (if ("FACT" %in% v) 6L else 4L) * rep * 3L)
  series * cen * rep * 16L
}
.prefixo <- function(modelo, variantes) {
  v <- toupper(strsplit(variantes, "+", fixed = TRUE)[[1]]); p <- modelo
  if ("S24" %in% v) p <- paste0(p, "_S24"); if ("S24W" %in% v) p <- paste0(p, "_S24W"); if ("CT" %in% v) p <- paste0(p, "_CT")
  if (nzchar(Sys.getenv("IC_TESTE_RAPIDO"))) p <- paste0("TESTE_", p)
  p
}

status <- function(modelo, h, variantes = "base") {
  tag <- .tag(modelo, h, variantes)
  janelas <- if (nzchar(Sys.getenv("IC_TESTE_RAPIDO"))) 2L else 10L   # N_WINDOWS do IC_Common.R
  ck <- file.path("checkpoints", sprintf("%s_%dpts", .prefixo(modelo, variantes), janelas * as.integer(h)))
  arq <- list.files(ck, pattern = "__[0-9]{2}_.*[.]rds$", full.names = TRUE)
  feitos <- length(arq)
  total <- .esperado(variantes)
  cat(sprintf("%s: %d de %d combinacoes (%.1f%%)\n", tag, feitos, total, 100 * feitos / total))
  # Ritmo medido pelas combinacoes concluidas na ultima hora (vale tambem apos retomar).
  if (feitos > 0L && feitos < total) {
    mt <- file.mtime(arq); agora <- Sys.time()
    janela <- min(60, as.numeric(difftime(agora, min(mt), units = "mins")))
    recentes <- sum(mt > agora - janela * 60)
    if (janela >= 2 && recentes > 0) {
      ritmo <- recentes / janela                      # combinacoes por minuto
      falta <- (total - feitos) / ritmo
      cat(sprintf("  ritmo: %.1f combinacoes/min | previsao de termino: %s (daqui a ~%s)\n", ritmo,
                  format(agora + falta * 60, "%d/%m %H:%M"),
                  if (falta >= 90) sprintf("%.1f h", falta / 60) else sprintf("%.0f min", falta)))
    }
  }
  logs <- list.files("logs", pattern = paste0("^", tag, "_parte[0-9]+[.]log$"), full.names = TRUE)
  logs <- logs[order(as.integer(sub(".*_parte([0-9]+)[.]log$", "\\1", logs)))]
  for (f in logs) {
    l <- tail(readLines(f, warn = FALSE), 1)
    cat(sprintf("  %s: %s\n", sub(".*_(parte[0-9]+)[.]log$", "\\1", f), if (length(l)) substr(l, 1, 90) else "(sem saida ainda)"))
  }
  jl <- file.path("logs", paste0(tag, "_junta.log"))
  if (file.exists(jl)) cat("  junta:", tail(readLines(jl, warn = FALSE), 1), "\n")
  proc <- .rodando(sprintf("roda_experimento.R --args %s %s %s ", modelo, h, variantes))
  cat(sprintf("  processos ativos: %d\n", length(proc)))
  invisible(feitos)
}

# Junta os checkpoints e grava o Excel (o lancar() ja faz isso sozinho no fim).
juntar <- function(modelo, h, variantes = "base") {
  system2("Rscript", c("roda_experimento.R", modelo, h, variantes, "junta"))
}

# Repete ate nao sobrar nada (a fila pode iniciar o proximo experimento no meio da parada).
parar <- function() {
  for (i in 1:10) {
    # pkill chamado direto (sem shell): um "sh -c" com estes nomes na linha de comando
    # seria morto pelo primeiro pkill e os seguintes nao rodariam.
    for (alvo in c("roda_fila.sh", "roda_paralelo.sh", "roda_experimento.R"))
      system2("pkill", c("-f", alvo), stdout = FALSE, stderr = FALSE)
    Sys.sleep(1)
    if (!length(.rodando("roda_fila.sh|roda_paralelo.sh|roda_experimento.R"))) break
  }
  restam <- length(.rodando("roda_fila.sh|roda_paralelo.sh|roda_experimento.R"))
  message(if (restam) paste("Ainda ha", restam, "processo(s); rode parar() de novo.") else
          "Tudo parado. O progresso esta salvo no checkpoint.")
}

# Fila: roda os experimentos em sequencia (cada um em n processos).
lancar_fila <- function(fila, n = n_padrao()) {
  if (length(.rodando("roda_fila.sh"))) { message("Ja existe uma fila rodando. Use status_fila() ou parar()."); return(invisible(FALSE)) }
  dir.create("logs", showWarnings = FALSE)
  linhas <- c("#!/usr/bin/env bash", "set -u")
  for (x in fila) {
    if (length(x) == 2L) x <- c(x, "base")
    linhas <- c(linhas, sprintf("echo \"$(date '+%%F %%T') INICIO %s %s %s\"", x[1], x[2], x[3]),
                sprintf("bash roda_paralelo.sh %s %s %s %d", shQuote(x[1]), shQuote(x[2]), shQuote(x[3]), as.integer(n)),
                sprintf("echo \"$(date '+%%F %%T') FIM %s %s %s (codigo $?)\"", x[1], x[2], x[3]))
  }
  linhas <- c(linhas, "echo \"$(date '+%F %T') FILA CONCLUIDA\"")
  writeLines(linhas, "roda_fila.sh")
  system(paste(.desliga(), "bash roda_fila.sh > logs/fila.log 2>&1 < /dev/null"), wait = FALSE)
  message(sprintf("Fila com %d experimentos lancada em %d processos cada. Acompanhe com status_fila().", length(fila), n))
  invisible(TRUE)
}

status_fila <- function() {
  if (!file.exists("logs/fila.log")) { cat("Nenhuma fila iniciada.\n"); return(invisible()) }
  l <- grep("INICIO|FIM|CONCLUIDA", readLines("logs/fila.log", warn = FALSE), value = TRUE)
  cat(l, sep = "\n")
  ini <- grep("INICIO", l, value = TRUE); fim <- grep("FIM", l, value = TRUE)
  if (length(ini) > length(fim)) {
    x <- strsplit(sub(".*INICIO ", "", tail(ini, 1)), " ")[[1]]
    cat("\nEm andamento:\n"); status(x[1], x[2], x[3])
  }
  invisible()
}
