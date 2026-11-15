# =====================================================================
# Prática 4 - Alinhamento múltiplo
# Script 4: concatenar os alinhamentos de 16S e COI pelo nome da espécie
# =====================================================================
# Rode na pasta de trabalho que contém a pasta "alinhamentos/"
# (confira com getwd(); mude com setwd() se precisar).

library(ape)
getwd()

# ---------------------------------------------------------------------
# 1. Ler os alinhamentos como matrizes (linhas = sequências,
#    colunas = posições do alinhamento)
# ---------------------------------------------------------------------
s16 <- as.matrix(read.FASTA("alinhamentos/16S_alinhado.fasta"))
coi <- as.matrix(read.FASTA("alinhamentos/COI_alinhado.fasta"))
dim(s16)
dim(coi)

# ---------------------------------------------------------------------
# 2. Separar o nome da espécie do número de acesso
#    Crossodactylus_schmidti_AY843579  ->  Crossodactylus_schmidti
#                                     ->  AY843579
# ---------------------------------------------------------------------
especie <- function(x) sub("_[A-Z]+_?[0-9]+$", "", x)
acesso  <- function(x) sub("^.*_([A-Z]+_?[0-9]+)$", "\\1", x)

# Tabela com o acesso usado para cada espécie em cada marcador
# (all = TRUE mantém também as espécies que só têm um dos marcadores)
tabela <- merge(data.frame(especie = especie(rownames(s16)),
                           acesso_16S = acesso(rownames(s16))),
                data.frame(especie = especie(rownames(coi)),
                           acesso_COI = acesso(rownames(coi))),
                by = "especie", all = TRUE)
tabela

# A partir daqui, cada linha passa a se chamar só pelo nome da espécie
rownames(s16) <- especie(rownames(s16))
rownames(coi) <- especie(rownames(coi))

# ---------------------------------------------------------------------
# 3. Completar cada matriz com as espécies que faltam (dado faltante "?")
# ---------------------------------------------------------------------
nomes <- sort(union(rownames(s16), rownames(coi)))
length(nomes)

setdiff(nomes, rownames(s16))   # espécies sem 16S
setdiff(nomes, rownames(coi))   # espécies sem COI

completar <- function(m, nomes) {
  faltam <- setdiff(nomes, rownames(m))
  if (length(faltam) == 0) return(m[nomes, ])
  vazio <- matrix("?", nrow = length(faltam), ncol = ncol(m),
                  dimnames = list(faltam, NULL))
  rbind(m, as.DNAbin(vazio))[nomes, ]
}

s16_completo <- completar(s16, nomes)
coi_completo <- completar(coi, nomes)

# ---------------------------------------------------------------------
# 4. Concatenar: o 16S ocupa as primeiras colunas e o COI vem em seguida
# ---------------------------------------------------------------------
concatenada <- cbind(s16_completo, coi_completo)
dim(concatenada)

# ---------------------------------------------------------------------
# 5. Gravar a matriz em três formatos, os limites de cada fragmento
#    (partições) e a tabela de acessos
# ---------------------------------------------------------------------
dir.create("matriz", showWarnings = FALSE)

write.FASTA(concatenada, "matriz/hylodidae_concatenado.fasta")
write.dna(concatenada, "matriz/hylodidae_concatenado.phy",
          format = "sequential", nbcol = -1, colsep = "")
write.nexus.data(concatenada, "matriz/hylodidae_concatenado.nex",
                 format = "dna", interleaved = FALSE)

particoes <- c(paste0("DNA, 16S = 1-", ncol(s16)),
               paste0("DNA, COI = ", ncol(s16) + 1, "-", ncol(s16) + ncol(coi)))
writeLines(particoes, "matriz/particoes.txt")
particoes

write.csv(tabela, "matriz/tabela_acessos.csv", row.names = FALSE, na = "")

list.files("matriz")
