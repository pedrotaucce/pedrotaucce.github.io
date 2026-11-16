# =====================================================================
# Prática 3 - Aquisição de sequências no GenBank
# Script 3: baixar sequências a partir de uma lista de números de acesso
# =====================================================================

library(ape)

dir.create("sequencias", showWarnings = FALSE)

# Números de acesso usados em um trabalho publicado (tabela do material
# suplementar). Inclui Limnomedusa macroglossa (Alsodidae) como grupo externo.
ids_16S <- c("KY007117.1", "KU215902.1", "KM390791.1", "AY843579.1", "AY263235.1",
             "MF624225.1", "KY627916.1", "KJ961575.1", "KJ961577.1", "MF624238.1",
             "AY843689.1", "KJ961586.1", "KJ961578.1", "DQ283072.1", "KJ961583.1",
             "KM390793.1", "KM390794.1", "KU495250.1", "KY002953.1", "KY627911.1",
             "KY627913.1", "KY627914.1")

ids_COI <- c("KU494387.1", "DQ502738.1", "KU494457.1", "KC593345.1", "KJ961562.1",
             "KJ961561.1", "KJ961558.1", "KJ961557.1", "KJ961556.1", "DQ502873.1",
             "KJ961565.1", "KU494594.1", "KJ961566.1", "DQ502839.1")

# read.GenBank() baixa as sequências e, no atributo "species", traz o nome
# do organismo segundo a taxonomia ATUAL do NCBI (campo ORGANISM).
s16 <- read.GenBank(ids_16S)
coi <- read.GenBank(ids_COI)

# Compare: a descrição foi escrita quando a sequência foi depositada;
# o campo ORGANISM acompanha as mudanças taxonômicas.
data.frame(acesso = names(coi),
           organismo = attr(coi, "species"),
           descricao = substr(attr(coi, "description"), 1, 60))

# Nomes no formato Genero_especie_acesso (sem o número da versão)
names(s16) <- paste(attr(s16, "species"), sub("\\.[0-9]+$", "", names(s16)), sep = "_")
names(coi) <- paste(attr(coi, "species"), sub("\\.[0-9]+$", "", names(coi)), sep = "_")
names(coi)

write.FASTA(s16, "sequencias/16S_ids.fasta")
write.FASTA(coi, "sequencias/COI_ids.fasta")

# ---------------------------------------------------------------------
# Alternativa com o rentrez (como na versão anterior do tutorial):
# grava o FASTA exatamente como o NCBI o fornece, sem renomear.
# ---------------------------------------------------------------------
# library(rentrez)
# fasta <- entrez_fetch(db = "nuccore", id = ids_16S, rettype = "fasta")
# cat(fasta, file = "sequencias/16S_ids_ncbi.fasta")
