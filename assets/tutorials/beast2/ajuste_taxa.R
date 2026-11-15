# =====================================================================
# Prática 6 - Árvore datada com o BEAST 2.7
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

# Parâmetros para a prior ORCucldMean no BEAUti (M = meanlog, S = sdlog)
ajuste$estimate

# Gráfico da distribuição ajustada, para comparar com o do BEAUti
x <- seq(0, 0.01, length = 100)
plot(x, dlnorm(x, meanlog = -5.374295, sdlog = 0.304072),
     type = "l", xlim = c(0.001, 0.01),
     xlab = "Taxa de substituição (subst./sítio/Ma)", ylab = "Densidade")
