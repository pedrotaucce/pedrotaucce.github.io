---
layout: page
title: "Árvore de espécies com o StarBEAST2"
description: "Árvore de espécies datada a partir de vários locos"
description_en: "Dated species tree from multiple loci with StarBEAST2."
order: 6
toc:
  sidebar: left
---

Tutorial da disciplina de pós-graduação [Métodos Computacionais em Sistemática e Filogeografia]({{ '/teaching/metodos_computacionais/' | relative_url }}).

<div class="alert alert-secondary" role="alert">
<strong>Arquivos da prática:</strong> <a href="{{ '/assets/tutorials/starbeast2/pratica_starbeast.zip' | relative_url }}">pratica_starbeast.zip</a>. Individualmente: <a href="{{ '/assets/tutorials/starbeast2/16S.nex' | relative_url }}">16S.nex</a>, <a href="{{ '/assets/tutorials/starbeast2/Tyr.nex' | relative_url }}">Tyr.nex</a>, <a href="{{ '/assets/tutorials/starbeast2/RAG1.nex' | relative_url }}">RAG1.nex</a>, <a href="{{ '/assets/tutorials/starbeast2/POMC.nex' | relative_url }}">POMC.nex</a>, <a href="{{ '/assets/tutorials/starbeast2/especies.dat' | relative_url }}">especies.dat</a>, <a href="{{ '/assets/tutorials/starbeast2/ajuste_taxa.R' | relative_url }}">ajuste_taxa.R</a>, <a href="{{ '/assets/tutorials/starbeast2/desafio_arvore.R' | relative_url }}">desafio_arvore.R</a>.
</div>

**Um aviso**

Esta é uma análise de tutorial, montada para rodar em menos de meia hora. Ao longo do texto, os quadros “**Na vida real**” indicam o que faríamos diferente numa análise para publicação.

## Introdução

Até aqui, todas as árvores da disciplina foram estimadas a partir de uma matriz concatenada, como se todos os genes tivessem a mesma história. Mas cada loco tem a sua própria genealogia, a *árvore de genes*, e ela nem sempre coincide com a história das espécies, a *árvore de espécies*. A principal razão é a separação incompleta de linhagens (*incomplete lineage sorting*): quando duas espécies se separam há pouco tempo, ou quando as populações ancestrais eram grandes, as cópias de um gene podem continuar misturadas por muitas gerações, e a árvore daquele gene pode mostrar um agrupamento diferente do das espécies (Degnan & Rosenberg 2009).

O modelo que descreve essa relação é o *coalescente multiespécies*: as árvores de genes evoluem “dentro” da árvore de espécies, como linhas dentro de tubos, e coalescem (se juntam num ancestral comum) mais depressa ou mais devagar dependendo do tamanho das populações. O StarBEAST2 (Ogilvie et al. 2017), um pacote do BEAST 2 (Bouckaert et al. 2019), estima ao mesmo tempo a árvore de cada loco e a árvore de espécies que melhor explica todas elas. Para isso, ele precisa saber a qual espécie pertence cada sequência, e é bom ter mais de um indivíduo por espécie.

Os dados desta prática são de Phyzelaphryninae, uma subfamília de rãzinhas de folhiço minúsculas com dois gêneros: *Phyzelaphryne*, da Amazônia, e *Adelophryne*, com espécies na Amazônia, no Escudo das Guianas e na Mata Atlântica (Fouquet et al. 2012). A filogenia do grupo mostrou três clados em *Adelophryne*, um no norte da Amazônia e dois na Mata Atlântica (norte e sul), e sugeriu que a subfamília surgiu na Amazônia e chegou à Mata Atlântica no início do Mioceno (Fouquet et al. 2012). São 14 espécies (12 de *Adelophryne* e 2 de *Phyzelaphryne*), de 1 a 4 indivíduos por espécie, e quatro locos: o 16S, mitocondrial, e três genes nucleares, tirosinase (Tyr), RAG1 e POMC. As sequências foram reunidas do GenBank; os números de acesso estão no seção [Números de acesso do GenBank](#ap-acessos).

## Instalação

1.  **BEAST 2, Tracer e FigTree.** Os mesmos do tutorial [Árvore datada com o BEAST 2.7]({{ '/tutorials/beast2/' | relative_url }}). O BEAST traz também o **DensiTree**, na mesma pasta do BEAUti.

2.  **StarBEAST2.** Abra o BEAUti, vá em *File \> Manage Packages*, selecione *starbeast2* e clique em *Install/Upgrade*. Feche e abra o BEAUti de novo para que o pacote seja carregado.

3.  **R**, com o pacote `fitdistrplus` (tutorial [Árvore datada com o BEAST 2.7]({{ '/tutorials/beast2/' | relative_url }})) e, para o desafio final, os pacotes `phytools` e `phyloch` (instruções no script `desafio_arvore.R`).

## Os dados

Desta vez, cada loco fica num arquivo NEXUS separado, porque cada um terá a sua árvore. Os nomes das sequências seguem o padrão

    A_baturitensis_JX298277        (16S)
    A_baturitensis_JX298096_a      (genes nucleares, cópia a)
    A_baturitensis_JX298096_b      (genes nucleares, cópia b)

com a inicial do gênero (`A` para *Adelophryne*, `P` para *Phyzelaphryne*), a espécie e um número de acesso do GenBank, que identifica o indivíduo. Nos genes nucleares, cada indivíduo aparece duas vezes, com `_a` e `_b`: são as duas cópias do gene (os dois alelos), separadas computacionalmente. O arquivo `especies.dat` é uma tabela com duas colunas, separadas por tabulação: o nome da sequência e a sua espécie.

|            | 16S | Tyr | RAG1 | POMC |
|:-----------|:---:|:---:|:----:|:----:|
| Sequências | 46  | 64  |  62  |  66  |
| Sítios     | 583 | 471 | 504  | 495  |

Nem todas as espécies têm os quatro locos: falta a Tyr de *A. amapaensis* e o RAG1 de *A. nordestina*. O StarBEAST2 aceita isso, desde que cada espécie esteja em pelo menos um loco.

## Montando a análise no BEAUti

### O modelo StarBEAST2

Abra o BEAUti e vá em *File \> Template \> StarBeast2*. As abas mudam: aparecem *Taxon sets*, *Gene Ploidy* e *Population Model*, que não existem no modelo padrão.

### Partitions

1.  Clique no botão *+*, no canto inferior esquerdo. Uma janela pergunta o que você quer adicionar: troque *Add Morphological Data* (a opção que vem marcada) por *Import Alignment* e clique em *OK*.

2.  Na pasta `dados`, selecione os quatro arquivos `.nex` (clique no primeiro e, com **Ctrl**, nos outros) e clique em *Open*. Não selecione o `especies.dat`: ele não é uma matriz.

3.  O BEAUti mostra uma linha para cada loco: `16S`, `POMC`, `RAG1` e `Tyr`, com o número de sequências (*Taxa*) e de sítios (*Sites*) da tabela acima.

Agora, os relógios. Para simplificar, os três genes nucleares vão compartilhar um único relógio; o 16S, que evolui bem mais rápido, terá o seu. Selecione as linhas `POMC`, `RAG1` e `Tyr` (clique na primeira e, com **Shift**, na última) e clique em *Link Clock Models*. Para facilitar a leitura, renomeie esse relógio: clique no nome na coluna *Clock Model* de uma dessas linhas, apague, escreva `nuclear` e tecle **Enter**.

Não clique em *Link Trees*: cada loco deve ter a sua árvore de genes. Os modelos de substituição (*Site Model*) também ficam separados.

> **Na vida real.** Os relógios não deveriam ser ligados: cada gene tem a sua própria taxa de evolução, e o ideal é deixar um relógio para cada loco. Ligamos os três genes nucleares só para a análise ter menos parâmetros e terminar no tempo da aula.

### Taxon sets

Aqui dizemos ao StarBEAST2 a qual espécie pertence cada sequência. Clique em *Guess*, marque *read from file*, clique em *Browse*, escolha o `especies.dat` e clique em *OK*. A coluna *Species/Population* passa a mostrar a espécie de cada sequência. Role a lista para conferir que nenhuma linha ficou vazia: devem ser 14 espécies.

### Tip Dates

Deixe como está.

### Gene Ploidy

A *ploidia* indica quantas cópias de cada loco há em cada indivíduo, comparado com os outros locos. Ela importa porque o tempo até as cópias de um gene coalescerem depende do número de cópias na população. Nos genes nucleares, cada indivíduo tem duas cópias (ploidia 2, o padrão). O DNA mitocondrial é herdado só da mãe e tem uma cópia por indivíduo: o número efetivo de cópias na população é um quarto do de um gene nuclear. Por isso, troque a ploidia de `geneTree.t:16S` para `0.5`; deixe os três genes nucleares em `2.0`.

### Population Model

Deixe *Analytical Population Size Integration*, o padrão. Com essa opção, o StarBEAST2 não estima o tamanho de cada população separadamente: ele “integra” sobre todos os valores possíveis, o que deixa a análise bem mais rápida. O único parâmetro de população que aparece no Tracer é a média, `popMean`.

### Site Model {#sec-sitemodel}

Desta vez, vamos usar o mesmo modelo nos quatro locos: HKY com variação de taxas entre sítios (Gama). Na lista da esquerda, clique em `16S` e, no *Gamma Site Model*:

1.  **desmarque** *estimate* em *Substitution Rate*, que vem marcado, e deixe o valor em 1,0;

2.  em *Gamma Category Count*, digite `4` e marque *estimate* em *Shape*;

3.  em *Subst Model*, escolha *HKY*. Deixe *Kappa* com *estimate* marcado e *Frequencies* em *Estimated*.

Repita para `POMC`, `RAG1` e `Tyr`. Se a caixa *estimate* de *Substitution Rate* voltar a se marcar sozinha, desmarque *Automatic set fix mean substitution rate flag* no menu *Mode* e desmarque a caixa de novo.

> **Por que não estimar a Substitution Rate?** Ela é a taxa relativa de cada loco, e quem cuida das taxas aqui são os relógios (seção [Clock Model](#sec-clock)). Se as duas forem estimadas ao mesmo tempo, os dados só “enxergam” o produto entre elas, e a cadeia não converge.

> **Na vida real.** O ideal seria não fixar nada disso. Primeiro, deixar o bModelTest (Bouckaert & Drummond 2017) escolher o modelo de cada loco, como no tutorial [Árvore datada com o BEAST 2.7]({{ '/tutorials/beast2/' | relative_url }}). Segundo, estimar as taxas relativas dos genes nucleares: marcar *estimate* em *Substitution Rate* só em `POMC`, `RAG1` e `Tyr`, deixando o 16S fixo em 1, para que cada gene nuclear tenha a sua taxa dentro do relógio compartilhado. Tiramos o bModelTest porque, nos nossos testes, vários parâmetros dele não convergiram em 10<sup>7</sup> gerações; e deixamos as taxas fixas para simplificar o passo a passo.

### Clock Model {#sec-clock}

Usaremos o relógio estrito (*Strict Clock*), o padrão. Como no tutorial [Árvore datada com o BEAST 2.7]({{ '/tutorials/beast2/' | relative_url }}), o tempo vem de uma *prior* para a taxa do 16S, ajustada a quatro estimativas publicadas para anuros (Lemmon et al. 2007; Fouquet et al. 2012; Gehara et al. 2017; Dufresnes et al. 2022) pelo script `ajuste_taxa.R`: média 0,00483 substituições por sítio por milhão de anos, com `meanlog` = -5,374295 e `sdlog` = 0,304072. A taxa dos genes nucleares será estimada a partir dos dados.

1.  Na lista da esquerda, clique em `16S`, digite `0.0048` em *Clock.rate* e marque *estimate*.

2.  Clique em `nuclear`, digite `0.001` em *Clock.rate* e marque *estimate*.

Se a caixa *estimate* estiver cinza, desmarque *Automatic set clock rate* no menu *Mode* e tente de novo.

> **Na vida real.** Como no tutorial [Árvore datada com o BEAST 2.7]({{ '/tutorials/beast2/' | relative_url }}), um relógio relaxado seria mais realista numa árvore com dezenas de milhões de anos. O StarBEAST2 tem relógios relaxados próprios para a árvore de espécies (na lista de *Clock Model*), mas eles têm mais parâmetros e demoram mais para convergir.

### Priors {#sec-priors}

1.  Deixe a *prior* da árvore de espécies no *Yule Model*, o padrão, e as *priors* de `popMean.Species` (*1/X*) e `speciationRate.t:Species` (*Uniform*) como estão.

2.  **Taxa do 16S.** Na linha `strictClockRate.c:16S`, troque a distribuição por *Log Normal* e clique no triângulo à esquerda do nome para abrir os detalhes. **Desmarque** *Mean In Real Space*, que vem marcado, e digite `M` = -5.374295 e `S` = 0.304072. O gráfico ao lado deve mostrar mediana 0,00463 e média 0,00485. Compare com o gráfico feito pelo script `ajuste_taxa.R`.

3.  **Taxa nuclear.** Na linha `strictClockRate.c:nuclear`, troque por *Log Normal*, desmarque *Mean In Real Space* e digite `M` = -6.907755 e `S` = 1.0. Esse *M* é o logaritmo de 0,001: a mediana da *prior* é 0,001, cerca de um quinto da taxa do 16S, porque genes nucleares codificantes evoluem bem mais devagar que o DNA mitocondrial. O desvio é grande (S = 1, intervalo de 95% de 0,00014 a 0,0071) porque não temos uma boa estimativa para esses genes: a *prior* só delimita a região, e quem decide são os dados. O fim do script `ajuste_taxa.R` calcula esses valores.

> **Por que desmarcar Mean In Real Space?** Com a caixa marcada, o BEAUti interpreta *M* como a média da própria taxa, e não como a média do seu logaritmo, que é o que o `fitdist` calcula. Um *M* negativo, como o nosso, nem faz sentido como taxa.

### MCMC

Deixe *Chain Length* em `10000000` (dez milhões de gerações), o padrão. Agora, o intervalo dos registros: clique no triângulo de cada registro (*tracelog*, *speciesTreeLogger*, *screenlog* e os quatro *treelog*, um para a árvore de cada loco) e troque *Log Every* de 5000 para `1000`. Assim, cada arquivo terá 10.000 amostras. Salve o arquivo com *File \> Save As*, com o nome `Phyzelaphryninae.xml`, na pasta `pratica_starbeast`.

> **Na vida real.** Uma corrida só não basta. O certo é rodar pelo menos duas réplicas independentes, conferir no Tracer se elas chegam aos mesmos valores e se todos os parâmetros têm ESS acima de 200, e aumentar a cadeia se for preciso. Análises com mais espécies e locos podem precisar de centenas de milhões de gerações.

## Rodando o BEAST

Abra o BEAST, escolha o arquivo `Phyzelaphryninae.xml` em *Input file* e clique em *Run*. A análise leva de 15 a 30 minutos, dependendo do computador. Ao terminar, a pasta terá seis arquivos novos:

- `starbeast.log`: os valores dos parâmetros a cada 1000 gerações;

- `species.trees`: as árvores de espécies, a cada 1000 gerações;

- `16S.trees`, `Tyr.trees`, `RAG1.trees` e `POMC.trees`: as árvores de genes, também a cada 1000 gerações.

## Resultados

**Tracer**

Carregue o `starbeast.log` e olhe a coluna *ESS*: praticamente todos os valores devem passar de 200. A exceção mais comum são as verossimilhanças de alguns genes (`treeLikelihood.16S`, por exemplo), que podem ficar um pouco abaixo sem afetar a árvore de espécies. Veja também:

- `TreeHeight.Species`: a idade da raiz da árvore de espécies, em milhões de anos;

- `strictClockRate.c:16S` e `strictClockRate.c:nuclear`: as taxas dos dois relógios. Quantas vezes o 16S é mais rápido que os genes nucleares?

- `popMean.Species`: o tamanho médio das populações, numa escala que combina tamanho populacional e tempo de geração.

**DensiTree**

Abra o DensiTree, vá em *File \> Load* e escolha o `species.trees`. Ele desenha todas as árvores amostradas umas sobre as outras: onde elas concordam, as linhas se sobrepõem e formam um traço forte; onde discordam, se espalham numa “nuvem”. É uma forma rápida de ver em que partes da árvore há incerteza. Se quiser, faça o mesmo com as árvores de genes (por exemplo, `16S.trees`) e compare.

**TreeAnnotator**

Como no tutorial [Árvore datada com o BEAST 2.7]({{ '/tutorials/beast2/' | relative_url }}): *Burn in percentage* 10, *Maximum clade credibility tree* e *Mean heights*. Em *Input Tree File*, escolha `species.trees`; em *Output File*, escreva `Phyzelaphryninae.tre`. Clique em *Run*.

**FigTree**

Abra o `Phyzelaphryninae.tre` e mostre, como no tutorial [Árvore datada com o BEAST 2.7]({{ '/tutorials/beast2/' | relative_url }}), a probabilidade posterior nos nós (*Node Labels*, *posterior*), as barras de 95% HPD das idades (*Node Bars*, *height_95%\_HPD*) e o eixo de tempo invertido (*Scale Axis*, *Reverse axis*). Repare que, ao contrário das práticas anteriores, os terminais agora são espécies, e não sequências.

## Desafio: a árvore no R {#sec-desafio}

O FigTree é ótimo para explorar a árvore, mas, para uma figura de artigo, muita gente prefere o R, em que tudo fica registrado num script. O `desafio_arvore.R` desenha a árvore de espécies datada sobre as épocas geológicas, com as barras de 95% HPD das idades e um círculo em cada nó, preto para probabilidade posterior acima de 0,95 e cinza entre 0,5 e 0,95.

1.  Copie o `desafio_arvore.R` para a pasta onde está o `Phyzelaphryninae.tre`.

2.  Abra o script no RStudio e vá em *Session \> Set Working Directory \> To Source File Location*.

3.  Instale os pacotes (só na primeira vez), com as linhas indicadas no começo do script. O `phyloch` não está no CRAN: ele é instalado do GitHub, pelo pacote `remotes`.

4.  Rode o script. Ele mostra a árvore na tela e salva os arquivos `Phyzelaphryninae.pdf` e `Phyzelaphryninae.svg`.

Conseguiu? Então tente:

1.  trocar a cor das barras de HPD (`col = adjustcolor("blue", ...)`) para cinza;

2.  mudar o limite de 0,95 para 0,85 e ver quais nós mudam de cor;

3.  descobrir em que época geológica cada um dos três clados de *Adelophryne* começou a se diversificar.

## Problemas comuns {#ap-problemas}

- ***File \> Template* não mostra *StarBeast2*.** O pacote não foi instalado, ou o BEAUti não foi reaberto depois da instalação.

- **O BEAUti pergunta sobre dados morfológicos (Mkv).** A opção *Add Morphological Data* ficou marcada ao clicar em *+*. Cancele e escolha *Import Alignment*. Se o BEAUti travar, feche e comece de novo.

- **Erro *Unsupported sequence file*.** Algum arquivo que não é matriz (como o `especies.dat`) foi selecionado junto com os `.nex`. Feche a mensagem; as quatro matrizes são importadas mesmo assim.

- **A coluna *Species/Population* tem linhas vazias.** O `especies.dat` foi editado e perdeu a tabulação, ou algum nome de sequência foi alterado.

- **A taxa do relógio não pode ser estimada.** Desmarque *Automatic set clock rate* no menu *Mode*.

- **Aparece `mutationRate` no Tracer, com ESS baixo.** A caixa *estimate* de *Substitution Rate* ficou marcada em algum loco (seção [Site Model](#sec-sitemodel)).

- **Não aparece `kappa` no Tracer.** O modelo ficou no padrão, JC69: a escolha de *HKY* em *Subst Model* não foi feita em algum loco.

- **As idades saem enormes ou minúsculas.** *Mean In Real Space* ficou marcado em alguma *prior* de taxa.

- **O BEAST avisa que o arquivo `.log` já existe.** Rode a nova réplica em outra pasta, ou escolha a opção de sobrescrever.

## Números de acesso do GenBank {#ap-acessos}

Nos genes nucleares, cada número corresponde a um indivíduo, representado na matriz pelas suas duas cópias (`_a` e `_b`).

| Espécie             | 16S                                    | Tyr                                    | RAG1                                   | POMC                                   |
|:--------------------|:---------------------------------------|:---------------------------------------|:---------------------------------------|:---------------------------------------|
| *A. adiastola*      | JX298299                               | JX298221                               | JX298167                               | JX298121                               |
| *A. amapaensis*     | JX298295                               | –                                      | JX298162                               | JX298116                               |
| *A. baturitensis*   | JX298277, JX298278, JX298279, JX298280 | JX298197, JX298199, JX298201, JX298202 | JX298145, JX298147, JX298149, JX298150 | JX298096, JX298098, JX298100, JX298101 |
| *A. glandulata*     | JX298287, JX298288, JX298289, MH304347 | JX298208, JX298209, JX298210, JX298211 | JX298154, JX298155                     | JX298107, JX298108, JX298109, JX298110 |
| *A. gutturosa*      | JX298300, JX298301, JX298302           | JX298222, JX298223, JX298224           | JX298168, JX298169, JX298170           | JX298122, JX298123, JX298124           |
| *A. maranguapensis* | JX298285, JX298286, KU495097, KU495098 | JX298206, JX298207                     | JX298152, JX298153                     | JX298105                               |
| *A. mucronata*      | JX298291, JX298292, MH304341, MH304342 | JX298213, JX298214                     | JX298158, JX298159                     | JX298112, JX298113                     |
| *A. nordestina*     | JX298284, MT731282, MT731283, MT731284 | JX298205                               | –                                      | JX298104                               |
| *A. pachydactyla*   | JX298283, MH304343, MH304344, MH304345 | JX298204                               | JX298151                               | JX298103                               |
| *A. patamona*       | JX298296, JX298297, JX298298           | JX298218, JX298219, JX298220           | JX298163, JX298164, JX298165, JX298166 | JX298117, JX298118, JX298119, JX298120 |
| *A.* sp. 4          | JX298290, MH304348                     | JX298212                               | JX298157                               | JX298111                               |
| *A.* sp. 8          | JX298293, JX298294, MH304333, MH304334 | JX298215, JX298216                     | JX298160, JX298161                     | JX298114, JX298115                     |
| *P. miriamae*       | EU186689, JX298303, JX298304, JX298305 | JX298225, JX298226, JX298228, JX298229 | JX298171, JX298172, JX298174, JX298175 | JX298125, JX298126, JX298128, JX298129 |
| *P. nimio*          | JX298308, JX298309, JX298310, JX298311 | JX298230, JX298231, JX298233, JX298235 | JX298176, JX298177, JX298179, JX298181 | JX298130, JX298131, JX298133, JX298134 |

</div>

## Referências

- Bouckaert, R. *et al.* 2019. BEAST 2.5: an advanced software platform for Bayesian evolutionary analysis. *PLoS Computational Biology* 15(4): e1006650.
- Bouckaert, R. R. & Drummond, A. J. 2017. bModelTest: Bayesian phylogenetic site model averaging and model comparison. *BMC Evolutionary Biology* 17: 42.
- Degnan, J. H. & Rosenberg, N. A. 2009. Gene tree discordance, phylogenetic inference and the multispecies coalescent. *Trends in Ecology & Evolution* 24(6): 332–340.
- Dufresnes, C. *et al.* 2022. Shedding light on taxonomic chaos: diversity and distribution of South Asian skipper frogs (Anura, Dicroglossidae, *Euphlyctis*). *Systematics and Biodiversity* 20(1): 2102686.
- Fouquet, A., Loebmann, D., Castroviejo-Fisher, S., Padial, J. M., Orrico, V. G. D., Lyra, M. L., Roberto, I. J., Kok, P. J. R., Haddad, C. F. B. & Rodrigues, M. T. 2012. From Amazonia to the Atlantic forest: molecular phylogeny of Phyzelaphryninae frogs reveals unexpected diversity and a striking biogeographic pattern emphasizing conservation challenges. *Molecular Phylogenetics and Evolution* 65(2): 547–561.
- Fouquet, A., Noonan, B. P., Rodrigues, M. T., Pech, N., Gilles, A. & Gemmell, N. J. 2012. Multiple Quaternary refugia in the eastern Guiana Shield revealed by comparative phylogeography of 12 frog species. *Systematic Biology* 61(3): 461–489.
- Gehara, M., Barth, A., Oliveira, E. F., Costa, M. A., Haddad, C. F. B. & Vences, M. 2017. Model-based analyses reveal insular population diversification and cryptic frog species in the *Ischnocnema parva* complex in the Atlantic forest of Brazil. *Molecular Phylogenetics and Evolution* 112: 68–78.
- Lemmon, E. M., Lemmon, A. R. & Cannatella, D. C. 2007. Geological and climatic forces driving speciation in the continentally distributed trilling chorus frogs (*Pseudacris*). *Evolution* 61(9): 2086–2103.
- Ogilvie, H. A., Bouckaert, R. R. & Drummond, A. J. 2017. StarBEAST2 brings faster species tree inference and accurate estimates of substitution rates. *Molecular Biology and Evolution* 34(8): 2101–2114.
- Werneck, F. P. *et al.* 2023. *Tutoriais: Filogeografia Estatística e Adaptação Ambiental*. Repositório Institucional do INPA. <https://repositorio.inpa.gov.br/handle/1/40376>.
