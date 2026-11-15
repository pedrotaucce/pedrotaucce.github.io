# =====================================================================
# Prática 3 - Aquisição de sequências no GenBank
# Script 2: baixar todas as sequências de 16S e COI de Hylodidae
# =====================================================================

library(ape)
library(rentrez)

# Opcional: chave pessoal da API do NCBI (aumenta o limite de 3 para 10
# requisições por segundo). Crie a sua em https://www.ncbi.nlm.nih.gov/account/
# e NUNCA a publique junto com o script.
# set_entrez_key("COLE_SUA_CHAVE_AQUI")

# Os arquivos serão gravados na subpasta "sequencias" da pasta de trabalho
dir.create("sequencias", showWarnings = FALSE)

# --- 1. O que vamos buscar ------------------------------------------------

# Os quatro gêneros de Hylodidae. Desde 2021, Megaelosia contém apenas
# M. goeldii; as demais espécies estão em Phantasmarana.
generos <- c("Crossodactylus", "Hylodes", "Megaelosia", "Phantasmarana")

# Termos de busca de cada marcador. O COI aparece no GenBank com vários nomes.
termo_coi <- paste("(COI[All Fields] OR COX1[All Fields]",
                   'OR "cytochrome oxidase subunit 1"[All Fields]',
                   'OR "cytochrome oxidase subunit I"[All Fields]',
                   'OR "cytochrome c oxidase subunit I"[All Fields])')
marcadores <- c("16S" = "16S[All Fields]", "COI" = termo_coi)

# Só sequências de 100 a 5000 pares de bases: exclui genomas mitocondriais
# completos (cerca de 17 mil pb), que contêm os dois genes e atrapalhariam
# o alinhamento.
filtro_tamanho <- "100:5000[SLEN]"

# "(Crossodactylus[Organism] OR Hylodes[Organism] OR ...)"
termo_generos <- paste0("(", paste0(generos, "[Organism]", collapse = " OR "), ")")

# --- 2. Função que busca e baixa em lotes --------------------------------

# A busca fica guardada no servidor do NCBI (web_history) e as sequências
# são baixadas em lotes, o que funciona mesmo com milhares de resultados.
baixar_fasta <- function(termo, arquivo, lote = 200) {
  busca <- entrez_search(db = "nuccore", term = termo, use_history = TRUE)
  if (busca$count == 0) return(0)
  if (file.exists(arquivo)) file.remove(arquivo)   # evita duplicar ao rodar de novo
  for (inicio in seq(0, busca$count - 1, by = lote)) {
    fasta <- entrez_fetch(db = "nuccore", web_history = busca$web_history,
                          rettype = "fasta", retstart = inicio, retmax = lote)
    cat(fasta, file = arquivo, append = TRUE)
    Sys.sleep(0.5)                                  # respeita o limite do NCBI
  }
  busca$count
}

# --- 3. Baixar: um arquivo por marcador -----------------------------------

for (m in names(marcadores)) {
  termo <- paste(termo_generos, "AND", marcadores[[m]], "AND", filtro_tamanho)
  arquivo <- file.path("sequencias", paste0("hylodidae_", m, ".fasta"))
  n <- baixar_fasta(termo, arquivo)
  message(m, ": ", n, " sequências gravadas em ", arquivo)
}

# --- 4. Quantas sequências há de cada gênero? -----------------------------

contagem <- sapply(generos, function(g) {
  sapply(names(marcadores), function(m) {
    Sys.sleep(0.5)
    termo <- paste(paste0(g, "[Organism]"), "AND", marcadores[[m]], "AND", filtro_tamanho)
    entrez_search(db = "nuccore", term = termo, retmax = 0)$count
  })
})
contagem

# --- 5. Conferir os arquivos ----------------------------------------------

s16 <- read.FASTA("sequencias/hylodidae_16S.fasta")
coi <- read.FASTA("sequencias/hylodidae_COI.fasta")
length(s16)                 # número de sequências
length(coi)
head(names(s16))            # primeiros cabeçalhos
summary(lengths(s16))       # comprimentos (pb)

# --- 6. Nomes curtos: Genero_especie_acesso -------------------------------

# A mesma expressão regular da Prática 1 (no R, cada \ é escrita \\).
# (?:UNVERIFIED: )? ignora o aviso que o GenBank põe em alguns registros.
renomear <- function(nomes) {
  sub("^([A-Z]+_?[0-9]+)\\.[0-9]+ (?:UNVERIFIED: )?([A-Z][a-z]+) ([a-z]+).*$",
      "\\2_\\3_\\1", nomes, perl = TRUE)
}
names(s16) <- renomear(names(s16))
names(coi) <- renomear(names(coi))
head(names(s16))

# Algum nome ficou repetido? (deve mostrar character(0))
names(s16)[duplicated(names(s16))]

write.FASTA(s16, "sequencias/hylodidae_16S_nomes.fasta")
write.FASTA(coi, "sequencias/hylodidae_COI_nomes.fasta")
