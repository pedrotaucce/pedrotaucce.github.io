# =====================================================================
# Prática 10 - DESAFIO: árvore de espécies datada com a escala geológica
# Script de Pedro P. G. Taucci, adaptado para a prática.
#
# Coloque este script na mesma pasta do Phyzelaphryninae.tre (a árvore
# MCC do TreeAnnotator) e, no RStudio, use
#   Session > Set Working Directory > To Source File Location
#
# Pacotes (só na primeira vez):
#   install.packages(c("phytools", "remotes", "XML", "colorspace"))
#   remotes::install_github("fmichonneau/phyloch")
# =====================================================================

library(phyloch)    # read.beast() e HPDbars()
library(phytools)   # plotTree() e make.transparent()

# O TreeAnnotator do BEAST 2.7 chama a árvore de "TREE_MCC_mean", mas o
# read.beast() do phyloch só encontra as estatísticas dos nós (posterior,
# HPD...) se ela se chamar "TREE1". Por isso, renomeamos numa cópia temporária.
arvore <- readLines("Phyzelaphryninae.tre")
arvore <- sub("^[[:space:]]*tree [^ ]+ = ", "tree TREE1 = [&R] ", arvore)
temporario <- tempfile(fileext = ".tre")
writeLines(arvore, temporario)

sptree <- read.beast(temporario)
posterior <- sptree$posterior

# Nós com PP > 0,95 em preto, entre 0,5 e 0,95 em cinza, o resto sem ponto
bg  <- ifelse(posterior > 0.95, "black",
              ifelse(posterior > 0.5, "darkgray", "transparent"))
col <- ifelse(posterior > 0.5, "black", "transparent")

# Cores e limites (em milhões de anos) das épocas geológicas
cols <- setNames(c(rgb(166, 216, 74, 0.5, max = 255), rgb(253, 167, 95, 0.5, max = 255),
                   rgb(253, 180, 108, 0.5, max = 255), rgb(253, 192, 122, 0.5, max = 255),
                   rgb(255, 255, 0, 0.5, max = 255), rgb(255, 255, 153, 0.5, max = 255),
                   rgb(255, 242, 174, 0.5, max = 255), rgb(254, 242, 224, 0.5, max = 255)),
                 c("Cretaceous", "Paleocene", "Eocene", "Oligocene", "Miocene",
                   "Pliocene", "Pleistocene", ""))

epochs <- matrix(c(100.5, 66.0,
                   66.0, 56.0,
                   56.0, 33.9,
                   33.9, 23.03,
                   23.03, 5.333,
                   5.333, 2.58,
                   2.58, 0.0117,
                   0.0117, 0),
                 8, 2, byrow = TRUE, dimnames = list(names(cols), c("start", "end")))

# Limite do eixo x: a maior idade do intervalo HPD, arredondada para cima
xmax <- ceiling(max(sptree$"height_95%_HPD_MAX", na.rm = TRUE) / 10) * 10
xlim <- c(xmax, -0.3 * xmax)   # à direita do zero fica o espaço para os nomes
ylim <- c(0.5, 1.2 * Ntip(sptree))

desenha <- function() {
  par(mar = c(3.1, 2.1, 1.1, 1.1))
  plot(NA, xlim = xlim, ylim = ylim, bty = "n", axes = FALSE, xlab = "", ylab = "")
  for (i in seq_along(cols)) {
    if (epochs[i, "end"] > xmax) next          # épocas mais antigas que o gráfico
    polygon(x = c(epochs[i, ], epochs[i, 2:1]),
            y = c(rep(ylim[1], 2), rep(Ntip(sptree) + 0.5, 2)),
            col = make.transparent(cols[i], 0.5), border = "gray")
    text(min(mean(epochs[i, ]), xmax), Ntip(sptree) + 0.6, names(cols)[i],
         adj = c(0, 0.25), srt = 45, cex = 1)
  }
  axis(1, at = seq(0, xmax, by = 10), cex.axis = 0.8)
  plotTree(ladderize(sptree), ftype = "reg", direction = "leftwards",
           xlim = xlim, ylim = ylim, mar = par()$mar, lwd = 3, add = TRUE)
  # O HPDbars avisa que as barras "excedem a região do gráfico" sempre que o eixo
  # está invertido (direction = "leftwards"); as barras são desenhadas certas,
  # então o aviso é silenciado.
  suppressWarnings(HPDbars(sptree, label = "height_95%_HPD",
                           col = adjustcolor("blue", alpha.f = 0.4), lwd = 8, broken = FALSE))
  nodelabels(frame = "none", pch = 21, col = col, bg = bg, cex = 1)
  mtext("milhões de anos", side = 1, at = xmax / 2, line = 2, cex = 0.9)
}

pdf("Phyzelaphryninae.pdf", height = 8, width = 10)
desenha()
dev.off()

svg("Phyzelaphryninae.svg", height = 8, width = 10)
desenha()
dev.off()

desenha()   # mostra também na tela
