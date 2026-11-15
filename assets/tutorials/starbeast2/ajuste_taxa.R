# =====================================================================
# Prática 10 - Árvore de espécies datada com o StarBEAST2
# Ajuste de uma lognormal às taxas de substituição do 16S de anuros
# (Lemmon et al. 2007; Fouquet et al. 2012; Gehara et al. 2017;
#  Dufresnes et al. 2022), em substituições por sítio por milhão de anos
# =====================================================================

install.packages("fitdistrplus")   # só na primeira vez
library(fitdistrplus)

taxas <- c(0.006, 0.00555, 0.00277, 0.005)

# Média das taxas: valor inicial da taxa do relógio no BEAUti
mean(taxas)

# Ajuste da lognormal por máxima verossimilhança
ajuste <- fitdist(taxas, "lnorm")
ajuste

# Parâmetros para a prior strictClockRate.c:16S no BEAUti (M = meanlog, S = sdlog)
ajuste$estimate

# Gráfico da distribuição ajustada, para comparar com o do BEAUti
x <- seq(0, 0.01, length = 100)
plot(x, dlnorm(x, meanlog = -5.374295, sdlog = 0.304072),
     type = "l", xlim = c(0.001, 0.01),
     xlab = "Taxa de substituição (subst./sítio/Ma)", ylab = "Densidade")

# Prior larga para o relógio nuclear (Tyr, RAG1 e POMC): mediana 0,001
# M = log(0.001) = -6.907755 e S = 1. Mediana e intervalo de 95%:
log(0.001)
qlnorm(c(0.025, 0.5, 0.975), meanlog = log(0.001), sdlog = 1)
