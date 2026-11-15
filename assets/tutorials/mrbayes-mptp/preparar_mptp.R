# =====================================================================
# preparar_mptp.R
# Prepara a árvore do MrBayes para o mPTP
#
# O mPTP precisa de uma árvore binária (sem politomias), enraizada e em
# formato Newick. Este script:
#   1. lê o consenso "allcompat" do MrBayes (anolis_ND2_allcompat.con.tre);
#   2. enraíza no grupo externo e retira os grupos externos;
#   3. resolve as politomias que sobrarem com ramos de comprimento zero;
#   4. grava a árvore (Newick) e o alinhamento sem os grupos externos
#      (FASTA), que o mPTP usa para calcular o comprimento mínimo de ramo.
#
# No RStudio: Session > Set Working Directory > To Source File Location
# Pacote (só na primeira vez): install.packages("ape")
# =====================================================================

library(ape)

arvore_mrbayes <- "anolis_ND2_allcompat.con.tre"
alinhamento    <- "anolis_ND2.fas"
externos       <- c("BR_ACB.CHUNB14548", "PL_BPN.96")

# 1. Lê o consenso (formato NEXUS). O arquivo do MrBayes traz a mesma árvore
#    duas vezes (com e sem as probabilidades posteriores); ficamos com a primeira
arvore <- read.nexus(arvore_mrbayes)
if (inherits(arvore, "multiPhylo")) arvore <- arvore[[1]]

# 2. Enraíza no primeiro grupo externo e retira os dois
arvore <- root(arvore, outgroup = externos[1], resolve.root = TRUE)
arvore <- drop.tip(arvore, externos)
arvore$node.label <- NULL                    # o mPTP não usa os valores de suporte

# 3. Resolve politomias (sequências idênticas costumam formar politomias)
set.seed(2026)
arvore <- multi2di(arvore)

cat("Terminais:", Ntip(arvore),
    "| binária:", is.binary(arvore),
    "| enraizada:", is.rooted(arvore), "\n")

# 4. Grava a árvore e o alinhamento só com o grupo interno
write.tree(arvore, "anolis_mptp.tre")

seqs <- read.dna(alinhamento, format = "fasta")
seqs <- seqs[!rownames(seqs) %in% externos, ]
write.dna(seqs, "anolis_mptp.fas", format = "fasta", colsep = "", nbcol = -1)

cat("Arquivos gravados: anolis_mptp.tre e anolis_mptp.fas\n")
