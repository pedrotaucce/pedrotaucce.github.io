---
layout: page
title: Aquisição de sequências no GenBank
description: "BLAST, buscas no site do GenBank e download de sequências com o R"
description_en: "Searching and downloading sequences from GenBank with BLAST and R."
order: 1
toc:
  sidebar: left
---

Tutorial da disciplina de pós-graduação [Métodos Computacionais em Sistemática e Filogeografia]({{ '/teaching/metodos_computacionais/' | relative_url }}) (Unesp Bauru).

<div class="alert alert-secondary" role="alert">
<strong>Arquivos da prática:</strong>
<a href="{{ '/assets/tutorials/genbank/pratica_genbank.zip' | relative_url }}">pratica_genbank.zip</a>
(scripts em R e a sequência do BLAST). Individualmente:
<a href="{{ '/assets/tutorials/genbank/01_preparar_R.R' | relative_url }}">01_preparar_R.R</a>,
<a href="{{ '/assets/tutorials/genbank/02_baixar_por_nome.R' | relative_url }}">02_baixar_por_nome.R</a>,
<a href="{{ '/assets/tutorials/genbank/03_baixar_por_acesso.R' | relative_url }}">03_baixar_por_acesso.R</a>,
<a href="{{ '/assets/tutorials/genbank/panthera_onca.fas' | relative_url }}">panthera_onca.fas</a>.
</div>

## O que é o GenBank?

O GenBank (<https://www.ncbi.nlm.nih.gov/genbank/>), mantido pelo NCBI (*National Center for Biotechnology Information*, EUA), é a base de dados que reúne as sequências de DNA disponíveis publicamente. Ele faz parte de um consórcio internacional (INSDC) com o ENA, europeu, e o DDBJ, japonês: as três bases trocam dados diariamente e compartilham os mesmos registros e números de acesso.

Quando você elabora um trabalho com sequências de DNA, provavelmente visitará o GenBank para adquirir sequências antes de começar as análises e, depois, para depositar as sequências que produziu: o depósito em um repositório público é exigido para a publicação de artigos. Cada registro recebe um **número de acesso** (por exemplo, `KJ961566`) e uma versão (`KJ961566.1`), que identificam a sequência de modo permanente.

Uma ferramenta muito útil associada ao GenBank é o **BLAST** (*Basic Local Alignment Search Tool*), que compara uma sequência com as sequências da base do NCBI e lista as mais parecidas. Ele serve, por exemplo, para verificar se uma sequência recém-produzida é de fato do organismo esperado ou se está contaminada.

## Usando o BLAST {#blast}

1.  Abra o arquivo `dados/panthera_onca.fas` em um editor de texto (Bloco de Notas, Sublime Text, Notepad++). Ele contém uma sequência parcial do gene mitocondrial 16S rRNA identificada como de onça-pintada (*Panthera onca*) do Pantanal.

2.  Acesse o BLAST de nucleotídeos: <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn>.

3.  Copie a sequência inteira, incluindo a linha de cabeçalho (`>`), e cole-a no campo *Enter accession number(s), gi(s), or FASTA sequence(s)*. Outra possibilidade é enviar o arquivo inteiro pelo botão *Escolher arquivo* (*Or, upload file*).

4.  Em *Choose Search Set*, mantenha a base padrão, *Core nucleotide database (core_nt)*; em *Program Selection*, mantenha *Highly similar sequences (megablast)*.

5.  Clique no botão *BLAST* e aguarde (de alguns segundos a alguns minutos).

Os resultados aparecem na aba *Descriptions*, ordenados da sequência mais parecida para a menos parecida. As colunas mais importantes são:

- *Scientific Name*: o organismo da sequência encontrada;

- *Query Cover*: a porcentagem da sua sequência que foi alinhada;

- *Per. Ident*: a porcentagem de bases idênticas no trecho alinhado;

- *E value*: o número de resultados com essa qualidade esperados ao acaso (quanto menor, melhor; valores como `0.0` ou `1e-150` indicam semelhança que não se explica pelo acaso);

- *Accession*: o número de acesso da sequência encontrada.

**Para refletir**

Qual foi a sequência mais próxima da nossa? Se você tivesse acabado de produzir esta sequência, o que isto poderia significar?

## Buscando sequências no site do GenBank

Se todas as sequências de DNA públicas estão no GenBank, é lá que procuramos quando precisamos de alguma. Acesse <https://www.ncbi.nlm.nih.gov/nuccore/> e busque pela sequência de algum organismo e de um fragmento. Se estiver sem ideia, vejamos as sequências do marcador mitocondrial 16S rRNA da rã-de-unhas-africana (*Xenopus laevis*). Digite no campo de busca:

    Xenopus laevis 16S

e clique em *Search*. A página mostra as sequências de *X. laevis* que mencionam o marcador: são centenas. Repare nos filtros da coluna da esquerda (por exemplo, *Sequence length*) e no quadro *Search details*, à direita, que mostra como o NCBI interpretou a sua busca, com *campos* entre colchetes, como `[Organism]`. Os mesmos campos podem ser usados diretamente:

    Xenopus laevis[Organism] AND 16S[All Fields] AND 100:5000[SLEN]

`[SLEN]` é o comprimento da sequência; aqui, de 100 a 5000 pares de bases, o que exclui os genomas mitocondriais completos (cerca de 17 mil pares de bases).

Baixar centenas de sequências clicando uma a uma daria um trabalho danado! Por isso, foram criados pacotes para a plataforma R que interagem com o GenBank.

Vamos supor que você encontrou uma espécie nova de rã que parece ser do gênero *Hylodes* e quer fazer uma árvore com dois marcadores mitocondriais, 16S rRNA e COI (citocromo c oxidase subunidade I), da família Hylodidae, para posicionar a nova espécie. Baixar no clique todas as sequências disponíveis para a família seria trabalhoso; vamos fazer isso no R.

## Um parêntese taxonômico: *Megaelosia* e *Phantasmarana*

Até 2021, a família Hylodidae era dividida em três gêneros: *Crossodactylus*, *Hylodes* e *Megaelosia*. Vittorazzi et al. (2021) mostraram que *Megaelosia* não era monofilético e descreveram o gênero *Phantasmarana* (“rãs fantasma”, pela raridade com que são encontradas). Desde então, *Megaelosia* contém uma única espécie, *M. goeldii*, e as demais foram transferidas para *Phantasmarana*, que hoje tem oito espécies (Frost, 2026), duas delas descritas depois da criação do gênero (de Sá et al., 2022):

| *Megaelosia*                 | *Phantasmarana*                                          |
|:-----------------------------|:---------------------------------------------------------|
| *M. goeldii* (Baumann, 1912) | *P. apuana* (Pombal, Prado, and Canedo, 2003)            |
|                              | *P. bocainensis* (Giaretta, Bokermann, and Haddad, 1993) |
|                              | *P. boticariana* (Giaretta and Aguiar, 1998)             |
|                              | *P. curucutuensis* de Sá et al., 2022                    |
|                              | *P. jordanensis* (Heyer, 1983)                           |
|                              | *P. lutzae* (Izecksohn and Gouvêa, 1987)                 |
|                              | *P. massarti* (De Witte, 1930)                           |
|                              | *P. tamuia* de Sá et al., 2022                           |

Isso afeta a busca de sequências de duas maneiras:

1.  A família passa a ter quatro gêneros, e *Phantasmarana* precisa estar na busca.

2.  As sequências depositadas antes de 2021 foram descritas como *Megaelosia*. A Taxonomia do NCBI atualiza o campo `ORGANISM` dos registros, que é o usado nas buscas com `[Organism]`; mas o texto descritivo do registro (a linha que aparece no cabeçalho do FASTA) foi escrito por quem depositou a sequência e pode continuar com o nome antigo. Ao montar uma matriz, confira sempre os nomes.

## Usando o R para adquirir sequências

### Preparando o R

Usaremos o R (<https://www.r-project.org/>) por meio do RStudio. Recomenda-se criar um *projeto* do RStudio para a prática (*File \> New Project \> New Directory*): o RStudio passa a usar a pasta do projeto como pasta de trabalho, e todos os arquivos são gravados ali, sem a necessidade de `setwd()` com caminhos que só existem no seu computador. Baixe os arquivos da prática (veja o quadro no início) e copie para a pasta do projeto as pastas `R` e `dados`.

Precisamos de dois pacotes: `ape` (Paradis & Schliep, 2019), para ler e gravar sequências, e `rentrez` (Winter, 2017), para conversar com as bases do NCBI. Abra o script `R/01_preparar_R.R` e execute-o (no RStudio, **Ctrl+Enter** executa a linha atual; **Ctrl+Shift+Enter**, o script inteiro):

    install.packages(c("ape", "rentrez"))
    library(ape)
    library(rentrez)
    packageVersion("ape")
    packageVersion("rentrez")

`library()` carrega o pacote; se ele não estiver instalado ou não funcionar, aparece uma mensagem de erro (por exemplo, `there is no package called ’ape’`). Se nenhum erro aparecer e as duas últimas linhas mostrarem as versões instaladas (algo como `[1] '5.8.1'`), você pode seguir adiante. Anote as versões: elas fazem parte da documentação de um trabalho. A instalação só precisa ser feita uma vez em cada computador; nas sessões seguintes basta carregar os pacotes com `library()`.

**Limite de acessos do NCBI**

O NCBI aceita até 3 requisições por segundo de cada endereço de internet. Em uma sala em que todos usam a mesma rede, esse limite pode ser atingido, e aparecerá um erro como `HTTP 429 Too Many Requests`: espere alguns segundos e execute o comando de novo. Quem cria uma conta no NCBI pode gerar uma chave pessoal da API (*API key*), que eleva o limite para 10 requisições por segundo, e informá-la ao R com `set_entrez_key("...")`. A chave é pessoal: não a publique junto com scripts.

### Baixando as sequências por nome

Abra o script `R/02_baixar_por_nome.R`. Vamos executá-lo por partes, entendendo cada uma.

#### Parte 1: o que buscar

```
library(ape)
library(rentrez)
dir.create("sequencias", showWarnings = FALSE)

generos <- c("Crossodactylus", "Hylodes", "Megaelosia", "Phantasmarana")

termo_coi <- paste("(COI[All Fields] OR COX1[All Fields]",
                   'OR "cytochrome oxidase subunit 1"[All Fields]',
                   'OR "cytochrome oxidase subunit I"[All Fields]',
                   'OR "cytochrome c oxidase subunit I"[All Fields])')
marcadores <- c("16S" = "16S[All Fields]", "COI" = termo_coi)
filtro_tamanho <- "100:5000[SLEN]"
termo_generos <- paste0("(", paste0(generos, "[Organism]", collapse = " OR "), ")")
```

O COI aparece no GenBank com vários nomes; por isso o termo junta as variantes com `OR` (`paste()` une os pedaços de texto, separados por espaço). A última linha monta o texto `(Crossodactylus[Organism] OR Hylodes[Organism] OR ...)`: digite `termo_generos` no console para vê-lo.

#### Parte 2: uma função para buscar e baixar

```
baixar_fasta <- function(termo, arquivo, lote = 200) {
  busca <- entrez_search(db = "nuccore", term = termo, use_history = TRUE)
  if (busca$count == 0) return(0)
  if (file.exists(arquivo)) file.remove(arquivo)
  for (inicio in seq(0, busca$count - 1, by = lote)) {
    fasta <- entrez_fetch(db = "nuccore", web_history = busca$web_history,
                          rettype = "fasta", retstart = inicio, retmax = lote)
    cat(fasta, file = arquivo, append = TRUE)
    Sys.sleep(0.5)
  }
  busca$count
}
```

`entrez_search()` faz a busca na base `nuccore` (nucleotídeos) e, com `use_history = TRUE`, deixa a lista de resultados guardada no servidor do NCBI. `entrez_fetch()` baixa as sequências em formato FASTA, em lotes de 200 (`retstart` indica a partir de qual resultado), o que funciona mesmo quando há milhares de resultados. O arquivo antigo é apagado antes, para que executar o script de novo não duplique as sequências, e `Sys.sleep(0.5)` faz uma pausa entre os lotes, respeitando o limite do NCBI.

#### Parte 3: baixar um arquivo por marcador

```
for (m in names(marcadores)) {
  termo <- paste(termo_generos, "AND", marcadores[[m]], "AND", filtro_tamanho)
  arquivo <- file.path("sequencias", paste0("hylodidae_", m, ".fasta"))
  n <- baixar_fasta(termo, arquivo)
  message(m, ": ", n, " sequências gravadas em ", arquivo)
}
```

Ao final, a pasta `sequencias` contém `hylodidae_16S.fasta` e `hylodidae_COI.fasta`, cada um com todas as sequências de Hylodidae do marcador. Anote quantas sequências foram gravadas: o número cresce à medida que novas sequências são depositadas.

#### Parte 4: quantas sequências há de cada gênero?

O script repete a busca para cada gênero e marcador, pedindo só a contagem (`retmax = 0`), e monta uma tabela:

```
contagem <- sapply(generos, function(g) {
  sapply(names(marcadores), function(m) {
    Sys.sleep(0.5)
    termo <- paste(paste0(g, "[Organism]"), "AND", marcadores[[m]], "AND", filtro_tamanho)
    entrez_search(db = "nuccore", term = termo, retmax = 0)$count
  })
})
contagem
```

#### Parte 5: conferir os arquivos

    s16 <- read.FASTA("sequencias/hylodidae_16S.fasta")
    coi <- read.FASTA("sequencias/hylodidae_COI.fasta")
    length(s16)
    head(names(s16))
    summary(lengths(s16))

`length()` informa o número de sequências; `names()`, os cabeçalhos; e `summary(lengths())`, o resumo dos comprimentos em pares de bases.

#### Parte 6: nomes curtos

Os cabeçalhos do GenBank são longos (`KJ961566.1 Megaelosia boticariana voucher CFBHT00425 cytochrome oxidase…`). Para transformá-los em `Genero_especie_acesso`, usamos o mesmo tipo de expressão regular visto na disciplina (no R, cada `\` é escrita `\\`, e `perl = TRUE` permite o grupo de não captura `(?:...)`, que aqui ignora o aviso `UNVERIFIED:` presente em alguns registros):

```
renomear <- function(nomes) {
  sub("^([A-Z]+_?[0-9]+)\\.[0-9]+ (?:UNVERIFIED: )?([A-Z][a-z]+) ([a-z]+).*$",
      "\\2_\\3_\\1", nomes, perl = TRUE)
}
names(s16) <- renomear(names(s16))
names(coi) <- renomear(names(coi))
names(s16)[duplicated(names(s16))]
write.FASTA(s16, "sequencias/hylodidae_16S_nomes.fasta")
write.FASTA(coi, "sequencias/hylodidae_COI_nomes.fasta")
```

A penúltima linha verifica se há nomes repetidos (deve mostrar `character(0)`).

### Baixando as sequências por número de acesso {#acesso}

Suponhamos que acaba de sair a descrição de uma nova espécie de *Phantasmarana*, como *P. curucutuensis* e *P. tamuia*, descritas em 2022 (de Sá et al., 2022), e que os autores fizeram uma ótima árvore, à qual você gostaria de acrescentar a sua espécie nova de *Hylodes*. Infelizmente, os autores não disponibilizaram a matriz usada no artigo, mas o material suplementar traz a lista de números de acesso de todas as sequências. Nesse caso, você pode baixar exatamente as mesmas sequências. Abra o script `R/03_baixar_por_acesso.R`:

```
library(ape)
ids_16S <- c("KY007117.1", "KU215902.1", "KM390791.1", "AY843579.1", "AY263235.1",
             "MF624225.1", "KY627916.1", "KJ961575.1", "KJ961577.1", "MF624238.1",
             "AY843689.1", "KJ961586.1", "KJ961578.1", "DQ283072.1", "KJ961583.1",
             "KM390793.1", "KM390794.1", "KU495250.1", "KY002953.1", "KY627911.1",
             "KY627913.1", "KY627914.1")
ids_COI <- c("KU494387.1", "DQ502738.1", "KU494457.1", "KC593345.1", "KJ961562.1",
             "KJ961561.1", "KJ961558.1", "KJ961557.1", "KJ961556.1", "DQ502873.1",
             "KJ961565.1", "KU494594.1", "KJ961566.1", "DQ502839.1")
s16 <- read.GenBank(ids_16S)
coi <- read.GenBank(ids_COI)
```

A função `read.GenBank()`, do `ape`, baixa as sequências e guarda, no atributo `species`, o nome do organismo segundo a taxonomia *atual* do NCBI. Compare esse nome com a descrição original de cada registro:

```
data.frame(acesso = names(coi),
           organismo = attr(coi, "species"),
           descricao = substr(attr(coi, "description"), 1, 60))
```

Por fim, damos às sequências nomes no formato `Genero_especie_acesso`, usando o nome atual do organismo, e gravamos os arquivos:

```
names(s16) <- paste(attr(s16, "species"), sub("\\.[0-9]+$", "", names(s16)), sep = "_")
names(coi) <- paste(attr(coi, "species"), sub("\\.[0-9]+$", "", names(coi)), sep = "_")
write.FASTA(s16, "sequencias/16S_ids.fasta")
write.FASTA(coi, "sequencias/COI_ids.fasta")
```

Agora você tem no seu computador os arquivos `16S_ids.fasta` e `COI_ids.fasta`, com as sequências desejadas. O script traz, comentada, a versão com o `rentrez` (`entrez_fetch()`), que grava o FASTA exatamente como o NCBI o fornece.

**Sem acesso ao NCBI?**

Se a rede da sala bloquear o NCBI ou o serviço estiver fora do ar, baixe os arquivos prontos [`16S_ids.fasta`]({{ '/assets/tutorials/genbank/16S_ids.fasta' | relative_url }}) e [`COI_ids.fasta`]({{ '/assets/tutorials/genbank/COI_ids.fasta' | relative_url }}): são as mesmas sequências da lista acima, já com os nomes no formato `Genero_especie_acesso`.

Com as sequências em mãos, o próximo passo é o alinhamento.

## Para refletir

1.  Na tabela `contagem`, quantas sequências de 16S e de COI há para *Megaelosia* e para *Phantasmarana*? Que espécie(s) aparece(m) quando se busca `Megaelosia[Organism]`?

2.  Procure, nos cabeçalhos de `hylodidae_COI.fasta` (ou nos nomes gerados pela função `renomear()`), sequências de *Phantasmarana*. Elas aparecem com esse nome? O que isso indica sobre a confiabilidade das linhas de descrição dos registros?

3.  Na a seção [Baixando as sequências por número de acesso](#acesso), compare as colunas `organismo` e `descricao` da tabela. Quais registros mudaram de nome?

4.  A lista de números de acesso inclui *Limnomedusa macroglossa*, que não pertence a Hylodidae. Por que ela foi incluída?

5.  Por que excluímos sequências com mais de 5000 pares de bases? O que aconteceria se um genoma mitocondrial completo entrasse no arquivo do 16S?

6.  Executar o script `02_baixar_por_nome.R` hoje e daqui a um ano produziria o mesmo arquivo? Como registrar, em um trabalho, quais sequências foram usadas?

## Problemas comuns

- **`HTTP 429` ou `Too Many Requests`:** limite de requisições do NCBI atingido (comum quando a turma toda usa a mesma rede). Aguarde alguns segundos e repita; se possível, use uma chave da API (`set_entrez_key()`).
- **`HTTP 500`, `HTTP 502` ou tempo esgotado:** instabilidade do NCBI. Repita mais tarde ou use os arquivos prontos da seção anterior.
- **`cannot open file`:** o arquivo não está na pasta de trabalho. Confira com `getwd()` e `list.files()`; no RStudio, abra o projeto da prática.
- **`there is no package called`:** o pacote não foi instalado. Execute `01_preparar_R.R`.

## Referências

- de Sá, F. P., Condez, T. H., Lyra, M. L., Haddad, C. F. B. & Malagoli, L. R. 2022. Unveiling the diversity of giant neotropical torrent frogs (Hylodidae): phylogenetic relationships, morphology, and the description of two new species. *Systematics and Biodiversity* 20(1): 2039318.
- Frost, D. R. 2026. *Amphibian Species of the World: an Online Reference*, versão 6.2. American Museum of Natural History, Nova York. <https://amphibiansoftheworld.amnh.org/>. Acesso em 9 nov. 2026.
- Paradis, E. & Schliep, K. 2019. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics* 35: 526–528.
- Vittorazzi, S. E., Augusto-Alves, G., Neves-da-Silva, D., Carvalho-e-Silva, A. M. P. T., Recco-Pimentel, S. M., Toledo, L. F., Lourenço, L. B. & Bruschi, D. P. 2021. Paraphyly in the giant torrent-frogs (Anura: Hylodidae: *Megaelosia*) and the description of a new genus. *Salamandra* 57(2): 274–284.
- Winter, D. J. 2017. rentrez: an R package for the NCBI eUtils API. *The R Journal* 9(2): 520–526.

