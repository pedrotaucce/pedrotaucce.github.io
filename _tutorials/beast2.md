---
layout: page
title: "Árvore datada com o BEAST 2.7"
description: "Árvore bayesiana datada com relógio molecular"
description_en: "Time-calibrated Bayesian tree with a molecular clock in BEAST 2.7."
order: 5
toc:
  sidebar: left
---

Tutorial da disciplina de pós-graduação [Métodos Computacionais em Sistemática e Filogeografia]({{ '/teaching/metodos_computacionais/' | relative_url }}).

<div class="alert alert-secondary" role="alert">
<strong>Arquivos da prática:</strong> <a href="{{ '/assets/tutorials/beast2/pratica_beast.zip' | relative_url }}">pratica_beast.zip</a>. Individualmente: <a href="{{ '/assets/tutorials/beast2/hylodidae_concatenado.nex' | relative_url }}">hylodidae_concatenado.nex</a>, <a href="{{ '/assets/tutorials/beast2/particoes.txt' | relative_url }}">particoes.txt</a>, <a href="{{ '/assets/tutorials/beast2/ajuste_taxa.R' | relative_url }}">ajuste_taxa.R</a>.
</div>

**Um aviso**

Esta é uma análise de tutorial, montada para rodar em poucos minutos. Ao longo do texto, os quadros “**Na vida real**” indicam o que faríamos diferente numa análise para publicação.

## Introdução

O BEAST (*Bayesian Evolutionary Analysis Sampling Trees*) (Bouckaert et al. 2019) é um programa de análise filogenética bayesiana que estima árvores com o tempo nos ramos (cronogramas). Há muitos tutoriais no site do BEAST 2 (<https://www.beast2.org>) e no *Taming the BEAST* (<https://taming-the-beast.org>). Uma análise bayesiana combina três ingredientes: a matriz de dados, o modelo de substituição de nucleotídeos e as distribuições *a priori* (*priors*) de cada parâmetro. O resultado não é uma única árvore, mas uma amostra de milhares de árvores e valores de parâmetros, a distribuição *a posteriori*.

Vamos usar a matriz concatenada de Hylodidae (16S + COI) e, para datar a árvore, uma *prior* para a taxa de substituição do 16S baseada em estimativas publicadas para anuros.

## Instalação

1.  **BEAST 2.** Baixe e instale a versão para o seu sistema em <https://www.beast2.org>. O pacote traz o BEAUti (que monta a análise), o BEAST (que roda a análise) e o TreeAnnotator (que resume as árvores).

2.  **bModelTest.** Abra o BEAUti, vá em *File \> Manage Packages*, selecione *bModelTest* (Bouckaert & Drummond 2017) e clique em *Install/Upgrade*. Feche e abra o BEAUti de novo para que o pacote seja carregado.

3.  **Tracer.** Baixe em <https://github.com/beast-dev/tracer/releases> (Rambaut et al. 2018).

4.  **FigTree.** O mesmo do tutorial [Árvore de máxima verossimilhança com o IQ-TREE 3]({{ '/tutorials/iqtree3/' | relative_url }}).

## A matriz com as partições

O BEAUti lê as partições de um bloco `sets` no fim do arquivo NEXUS. Abra o `hylodidae_concatenado.nex` no Sublime Text ou no Notepad++ e cole no fim do arquivo, depois do `END;` que fecha a matriz, o bloco abaixo, usando os números do seu `particoes.txt`:

      ;
    END;

    begin sets;
      charset p16S = 1-2475;
      charset COI  = 2476-3187;
    end;

Salve com o nome `hylodidae_particionado.nex`. O nome `p16S` começa com uma letra porque, no formato NEXUS, um nome não pode começar com número.

> **Na vida real.** No tutorial [Árvore de máxima verossimilhança com o IQ-TREE 3]({{ '/tutorials/iqtree3/' | relative_url }}), o COI foi dividido nas três posições do códon. Aqui usamos só duas partições para ter menos parâmetros e uma análise mais rápida.

## Montando a análise no BEAUti

### Partitions

Vá em *File \> Import Alignment* e abra o `hylodidae_particionado.nex`. O BEAUti mostra uma linha para cada partição: `p16S` (2475 sítios) e `COI` (712).

Os dois fragmentos vêm do DNA mitocondrial, que é herdado como uma unidade só: eles compartilham a mesma árvore. Selecione as duas linhas (clique na primeira e, com **Shift**, na segunda) e clique em *Link Trees*. A coluna *Tree* passa a mostrar `p16S` nas duas linhas. Os modelos de substituição (*Site Model*) e os relógios (*Clock Model*) continuam separados: cada marcador terá o seu.

### Tip Dates

Deixe como está.

### Site Model

Em vez de escolher um modelo de substituição, vamos usar o pacote bModelTest (Bouckaert & Drummond 2017), que estima o modelo junto com a árvore: durante a análise, ele passa por vários modelos e calcula a probabilidade de cada um. Na lista da esquerda, clique em `p16S` e troque *Gamma Site Model* por *BEAST Model Test*. Faça o mesmo para `COI`. Deixe *Mutation Rate* em 1,0, sem marcar *estimate*: as diferenças de taxa entre os dois marcadores ficarão por conta dos relógios.

### Clock Model {#sec-clock}

Usaremos o relógio estrito (*Strict Clock*), o padrão do BEAUti, que supõe a mesma taxa de substituição em todos os ramos da árvore. Cada marcador terá o seu relógio, porque o COI evolui mais rápido que o 16S.

> **Na vida real.** Numa árvore com dezenas de milhões de anos, como esta, dificilmente a taxa é a mesma em todos os ramos. Um relógio relaxado, em que cada ramo tem a sua taxa (por exemplo, o *Optimised Relaxed Clock*, do pacote ORC), seria provavelmente mais adequado, mas tem mais parâmetros e demora mais para convergir.

Para transformar distâncias genéticas em tempo, é preciso informar ao BEAST algo sobre o tempo: a idade de um nó (por um fóssil ou um evento geológico) ou a taxa de substituição. Vamos usar uma *prior* para a taxa do 16S, ajustada a quatro estimativas publicadas para anuros (Lemmon et al. 2007; Fouquet et al. 2012; Gehara et al. 2017; Dufresnes et al. 2022): 0,006; 0,00555; 0,00277 e 0,005 substituições por sítio por milhão de anos. O script `ajuste_taxa.R` ajusta uma distribuição lognormal a esses valores:

    install.packages("fitdistrplus")
    library(fitdistrplus)
    taxas <- c(0.006, 0.00555, 0.00277, 0.005)
    mean(taxas)
    ajuste <- fitdist(taxas, "lnorm")
    ajuste

A média das taxas (0,00483) será o valor inicial da taxa do 16S, e os parâmetros da lognormal, `meanlog` = -5,374295 e `sdlog` = 0,304072, serão usados na *prior* (seção [Priors](#sec-priors)). Como as taxas estão em milhões de anos, as idades da árvore também sairão em milhões de anos.

O 16S calibra o tempo; a taxa do COI será estimada a partir dos dados. No BEAUti:

1.  No menu *Mode*, desmarque *Automatic set clock rate*.

2.  Na aba *Clock Model*, clique em `p16S`, digite `0.00483` em *Mean clock rate* e marque *estimate*.

3.  Clique em `COI`, digite `0.01` em *Mean clock rate* e marque *estimate*.

### Priors {#sec-priors}

1.  Na primeira linha, `Tree.t:p16S`, deixe o *Yule Model*, a *prior* de árvore mais simples: as linhagens só se originam, com uma única taxa de especiação. (O *Birth Death Model* também permite a extinção, mas tem um parâmetro a mais.)

2.  **Taxa do 16S.** Na linha `clockRate.c:p16S`, troque *Uniform* por *Log Normal* e clique no triângulo à esquerda do nome para abrir os detalhes. Digite `M` = -5.374295 e `S` = 0.304072 (sem marcar *Mean In Real Space*). Ao lado, o BEAUti desenha a distribuição: ela deve ter mediana 0,00463 e média 0,00485. Compare com o gráfico feito no R pelo fim do script `ajuste_taxa.R`.

3.  **Taxa do COI.** Na linha `clockRate.c:COI`, troque *Uniform* por *Log Normal* e digite `M` = -4.68 e `S` = 1.0. Em anfíbios, o COI acumula diferenças cerca de duas vezes mais depressa que o 16S (Vences et al. 2005). Por isso a mediana dessa *prior* é o dobro da do 16S (-5,374 + ln 2 ≈ -4,68, mediana 0,00928), mas o desvio é bem maior (S = 1, contra 0,30 do 16S): o intervalo de 95% vai de 0,0013 a 0,066, ou seja, de cerca de 0,3 a 14 vezes a taxa do 16S. A *prior* do 16S é estreita porque é ela que calibra o tempo; a do COI é larga porque só ajuda, e quem decide são os dados.

4.  **Monofilia de Hylodidae.** Clique em *+ Add Prior*, escolha *MRCA prior* e clique em *OK*. Em *Taxon set label*, escreva `Hylodidae`. Selecione todos os táxons da lista da esquerda (**Ctrl+A**) e passe-os para a direita com o botão `>>`; depois, selecione *Limnomedusa macroglossa* na lista da direita e devolva-a com `<<`. Clique em *OK* e, na nova linha (a que começa com o botão *Hylodidae*), marque *monophyletic*.

    Temos um único grupo externo, e só um relógio não basta para posicionar bem a raiz da árvore. Por isso forçamos a monofilia de Hylodidae, que é bem sustentada por estudos anteriores (Vittorazzi et al. 2021; de Sá et al. 2022). Bônus: o Tracer passará a mostrar a idade do ancestral comum mais recente de Hylodidae (`mrca.age(Hylodidae)`).

5.  Deixe as outras *priors* como estão.

> **Na vida real.** O ideal é não forçar a monofilia e incluir vários grupos externos, de modo que os dados mostrem se o grupo é monofilético.

### MCMC

Em *Chain Length*, digite `1000000` (um milhão de gerações). Mantenha o registro a cada 1000 gerações (*Log Every*, dentro de *tracelog* e *treelog*), o valor padrão. Salve o arquivo com *File \> Save As*, com o nome `hylodidae.xml`, na pasta `pratica_beast`.

> **Na vida real.** Um milhão de gerações é pouco: o ideal é rodar pelo menos 10 milhões (10<sup>7</sup>) e conferir no Tracer se todos os parâmetros convergiram.

## Rodando o BEAST

Abra o BEAST, clique em *Choose file...* ao lado de *Input file*, escolha o arquivo `hylodidae.xml` e clique em *Run*. A janela mostra o andamento da análise, que leva poucos minutos. Ao terminar, a pasta terá dois arquivos novos:

- `hylodidae.log`: os valores dos parâmetros a cada 1000 gerações;

- `hylodidae-p16S.trees`: as árvores amostradas a cada 1000 gerações. O nome vem da árvore, que o BEAUti batizou com o nome da primeira partição.

É recomendável rodar pelo menos mais uma réplica com outra semente (*Random number seed:*, na janela do BEAST), em outra pasta, para conferir se as duas chegam ao mesmo resultado.

## Resultados

**Tracer**

Abra o Tracer e carregue o `hylodidae.log` com o botão *+* no canto superior esquerdo. Na coluna *ESS* (*effective sample size*), procure valores acima de 200; valores menores aparecem em destaque e indicam que a cadeia precisa rodar por mais tempo. Se você rodou duas réplicas, carregue as duas e compare as médias.

**TreeAnnotator**

O TreeAnnotator resume as milhares de árvores do arquivo `.trees` em uma só, a árvore de máxima credibilidade de clados (*Maximum clade credibility tree*).

1.  Em *Burn in percentage*, deixe 10 (as primeiras 10% das árvores são descartadas).

2.  Em *Target tree type*, deixe *Maximum clade credibility tree*.

3.  Em *Node heights*, troque *Common Ancestor heights* por *Mean heights*.

4.  Em *Input Tree File*, escolha `hylodidae-p16S.trees` (*Choose file...*); em *Output File*, escreva `hylodidae_mcc.tree`. Clique em *Run*.

**FigTree**

Abra o `hylodidae_mcc.tree` no FigTree. No painel da esquerda:

- marque *Node Labels* e, em *Display*, escolha *posterior*: é a probabilidade posterior de cada clado;

- marque *Node Bars* e, em *Display*, escolha *height_95%\_HPD*: as barras mostram o intervalo de 95% de maior densidade da idade de cada nó;

- marque *Scale Axis* e, dentro dele, *Reverse axis*, para ler as idades em milhões de anos antes do presente.

## Perguntas

1.  Todos os parâmetros chegaram a ESS acima de 200? Quais ficaram abaixo? O que você faria para resolver?

2.  Qual a idade estimada (média e intervalo de 95%) do ancestral comum mais recente de Hylodidae?

3.  No Tracer, compare as médias de `clockRate.COI` e `clockRate.p16S`. Quantas vezes o COI evolui mais rápido que o 16S? O resultado está dentro do que esperávamos pela *prior*?

## Problemas comuns {#ap-problemas}

- ***BEAST Model Test* não aparece.** O pacote não foi instalado, ou o BEAUti não foi reaberto depois da instalação.

- **O BEAUti mostra uma partição só.** O bloco `begin sets; ... end;` não foi colado no fim do arquivo, ou falta algum ponto e vírgula.

- **A taxa do relógio não pode ser estimada.** Desmarque *Automatic set clock rate* no menu *Mode*.

- **O BEAST avisa que o arquivo `.log` já existe.** Já há uma análise com o mesmo nome na pasta. Rode a nova réplica em outra pasta, ou escolha a opção de sobrescrever.

## Referências

- Bouckaert, R. *et al.* 2019. BEAST 2.5: an advanced software platform for Bayesian evolutionary analysis. *PLoS Computational Biology* 15(4): e1006650.
- Bouckaert, R. R. & Drummond, A. J. 2017. bModelTest: Bayesian phylogenetic site model averaging and model comparison. *BMC Evolutionary Biology* 17: 42.
- de Sá, F. P., Condez, T. H., Lyra, M. L., Haddad, C. F. B. & Malagoli, L. R. 2022. Unveiling the diversity of Giant Neotropical Torrent frogs (Hylodidae): phylogenetic relationships, morphology, and the description of two new species. *Systematics and Biodiversity* 20(1): 2039318.
- Dufresnes, C. *et al.* 2022. Shedding light on taxonomic chaos: diversity and distribution of South Asian skipper frogs (Anura, Dicroglossidae, *Euphlyctis*). *Systematics and Biodiversity* 20(1): 2102686.
- Fouquet, A., Noonan, B. P., Rodrigues, M. T., Pech, N., Gilles, A. & Gemmell, N. J. 2012. Multiple Quaternary refugia in the eastern Guiana Shield revealed by comparative phylogeography of 12 frog species. *Systematic Biology* 61(3): 461–489.
- Gehara, M., Barth, A., Oliveira, E. F., Costa, M. A., Haddad, C. F. B. & Vences, M. 2017. Model-based analyses reveal insular population diversification and cryptic frog species in the *Ischnocnema parva* complex in the Atlantic forest of Brazil. *Molecular Phylogenetics and Evolution* 112: 68–78.
- Lemmon, E. M., Lemmon, A. R. & Cannatella, D. C. 2007. Geological and climatic forces driving speciation in the continentally distributed trilling chorus frogs (*Pseudacris*). *Evolution* 61(9): 2086–2103.
- Rambaut, A., Drummond, A. J., Xie, D., Baele, G. & Suchard, M. A. 2018. Posterior summarization in Bayesian phylogenetics using Tracer 1.7. *Systematic Biology* 67(5): 901–904.
- Vences, M., Thomas, M., Bonett, R. M. & Vieites, D. R. 2005. Deciphering amphibian diversity through DNA barcoding: chances and challenges. *Philosophical Transactions of the Royal Society B* 360(1462): 1859–1868.
- Vittorazzi, S. E., Augusto-Alves, G., Neves-da-Silva, D., Carvalho-e-Silva, A. M. P. T., Recco-Pimentel, S. M., Toledo, L. F., Lourenço, L. B. & Bruschi, D. P. 2021. Paraphyly in the giant torrent-frogs (Anura: Hylodidae: *Megaelosia*) and the description of a new genus. *Salamandra* 57(2): 274–284.
