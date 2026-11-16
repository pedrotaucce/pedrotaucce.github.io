---
layout: page
title: "Máxima verossimilhança com dados de RAD-seq no IQ-TREE 3"
description: "Matriz completa, SNPs e fatores de concordância em <em>Pedicularis</em>"
description_en: "Maximum likelihood with RAD-seq data in IQ-TREE 3: full matrix, SNPs and concordance factors."
order: 4
toc:
  sidebar: left
---

Tutorial da disciplina de pós-graduação [Métodos Computacionais em Sistemática e Filogeografia]({{ '/teaching/metodos_computacionais/' | relative_url }}).

<div class="alert alert-secondary" role="alert">
<strong>Arquivos da prática:</strong> <a href="{{ '/assets/tutorials/iqtree3-radseq/pratica_iqtree_radseq.zip' | relative_url }}">pratica_iqtree_radseq.zip</a>. Individualmente: <a href="{{ '/assets/tutorials/iqtree3-radseq/pedicularis.snps' | relative_url }}">pedicularis.snps</a>, <a href="{{ '/assets/tutorials/iqtree3-radseq/pedicularis.usnps' | relative_url }}">pedicularis.usnps</a>.<br>O arquivo <code>pedicularis.phy</code> (40 MB) é grande demais para este site: <a href="https://drive.google.com/file/d/19NrBcR2-rtG0y2-I4DS_we8FFbLwWgJn/view?usp=sharing">baixe-o pelo Google Drive</a>.
</div>

**O roteiro**

Vamos estimar a filogenia das 13 amostras de *Pedicularis* por máxima verossimilhança no IQ-TREE 3 (Wong et al. 2025), com três matrizes que o ipyrad gerou a partir dos mesmos dados, e comparar o que cada uma mostra. No fim, vamos medir o conflito entre os sítios com os fatores de concordância de sítio.

**Programas**

IQ-TREE 3 e FigTree.

## Antes de começar: o tutorial do ipyrad

Esta prática parte das saídas do tutorial do ipyrad com os dados de RAD-seq de *Pedicularis* (13 amostras, Eaton & Ree (2013)), feito antes desta prática: <https://ipyrad.readthedocs.io/en/master/pedicularis_.html>. Três avisos sobre aquele tutorial:

1.  **O download pode falhar.** O `ipyrad --download` usa o `fastq-dump` do sra-tools, e o ambiente do ipyrad pode vir com uma versão de 2019 (2.9.6), que não lê o formato atual dos arquivos do NCBI. O erro aparece como `rcColumn, rcNotFound`. Para consertar, dentro do ambiente do ipyrad:

        conda install -c conda-forge -c bioconda "sra-tools>=3"
        fastq-dump --version

    A versão deve começar com 3. Apague os arquivos `.sra` que tiverem sobrado na pasta e rode o `ipyrad --download` de novo.

2.  **O arquivo de SNPs mudou de nome.** Quando o tutorial manda olhar o `pedicularis.snps.phy`, na versão atual do ipyrad o arquivo é o `pedicularis.snps` (o mesmo formato PHYLIP, só sem o `.phy` no nome).

3.  **Cuidado com os parâmetros 25 e 26.** O tutorial muda o 26 (`trim_loci`) para `0, 5, 0, 0`. O 25 (`trim_reads`), logo acima, tem o mesmo formato: se o valor for parar nele, cada sequência é cortada para 5 pb e o passo 2 descarta todas (`No reads passed filtering`).

## Introdução

**As três matrizes**

O ipyrad escreve a mesma informação de vários jeitos. Vamos usar três:

| Arquivo             | O que tem                                           | Correção +ASC? |
|:--------------------|:----------------------------------------------------|:--------------:|
| `pedicularis.phy`   | todos os sítios de todos os loci, concatenados      |      não       |
| `pedicularis.snps`  | só os sítios variáveis: todos os SNPs de cada locus |      sim       |
| `pedicularis.usnps` | um SNP sorteado por locus (*unlinked*)              |      sim       |

**Por que corrigir**

Uma matriz só de SNPs não tem sítios invariáveis, porque eles foram tirados. Se o modelo não souber disso, ele interpreta “todo sítio variou” como uma taxa de substituição altíssima, e os ramos saem longos demais. A correção do viés de averiguação (*ascertainment bias correction*, `+ASC` no IQ-TREE (Lewis 2001)) calcula a verossimilhança de cada sítio *sabendo que ele é variável*: divide a verossimilhança do sítio pela probabilidade de um sítio qualquer ser variável,

$$L_i^{\text{ASC}} = \frac{L_i}{1 - \left[L(\text{AAA}\ldots) + L(\text{CCC}\ldots) +
L(\text{GGG}\ldots) + L(\text{TTT}\ldots)\right]},$$

em que L(AAA…) é a verossimilhança de uma coluna com A em todas as amostras.

**Por que olhar além do *bootstrap***

Com centenas de milhares de sítios, o *bootstrap* ultrarrápido (UFBoot) (Hoang et al. 2018) tende a dar 100 em quase todo nó, mesmo onde há conflito real entre os sítios (por separação incompleta de linhagens, por exemplo). O fator de concordância de sítio (sCF) (Minh et al. 2020; Mo et al. 2023) mede outra coisa: a porcentagem dos sítios informativos que apoiam cada ramo, contra as duas alternativas. Vamos calcular os dois.

## Instalação e arquivos

**O IQ-TREE 3**

Os executáveis estão em <https://iqtree.github.io>, e o programa se chama `iqtree3`. Use a versão 3.1.4 ou mais nova: foi a versão testada neste tutorial. Como as saídas do ipyrad estão no Linux (ou no WSL) ou no Mac, instale-o no mesmo lugar: baixe o arquivo do seu sistema (por exemplo, `iqtree-3.1.4-Linux-intel.tar.gz` no Linux e no WSL) e descompacte com `tar xzf` seguido do nome do arquivo. O programa é o `iqtree3` da pasta `bin`.

**O `PATH`**

Até aqui, rodamos programas baixados com `./`, de dentro da pasta onde eles estão. O jeito mais prático é colocá-los numa pasta onde o terminal sempre procura comandos. Essa lista de pastas é a variável de ambiente `PATH`; veja quais são com

    echo $PATH

A pasta `/usr/local/bin` está nessa lista no Linux, no WSL e no Mac. Mova o `iqtree3` para lá (o `sudo` pede a senha de administrador; no WSL, é a senha que você criou ao instalar o Ubuntu):

    sudo mkdir -p /usr/local/bin
    sudo mv iqtree-3.1.4-Linux-intel/bin/iqtree3 /usr/local/bin/
    iqtree3 --version

O primeiro comando só cria a pasta, se ela ainda não existir (nos Macs com chip Apple, ela às vezes não existe). No segundo, troque `iqtree-3.1.4-Linux-intel` pelo nome da pasta que saiu do arquivo do seu sistema. No Mac, antes do `mv`, rode `xattr -d com.apple.quarantine` seguido do caminho do `iqtree3`, para tirar o bloqueio que o Mac põe em programas baixados. Agora o `iqtree3` funciona de qualquer pasta.

**Os arquivos**

Crie uma pasta para a prática e copie para ela as três matrizes, que ficam na pasta de saídas do ipyrad (ou use as da pasta da prática):

    mkdir iqtree
    cp analysis-ipyrad/pedicularis_outfiles/pedicularis.phy iqtree/
    cp analysis-ipyrad/pedicularis_outfiles/pedicularis.snps iqtree/
    cp analysis-ipyrad/pedicularis_outfiles/pedicularis.usnps iqtree/
    cd iqtree

**O grupo externo**

As duas amostras de *P. przewalskii* enraízam as árvores. Os nomes vêm com o número do SRA, e o IQ-TREE recebe os dois separados por vírgula, sem espaço:

    -o 32082_przewalskii_SRR1754729,33588_przewalskii_SRR1754727

Para conferir os nomes, olhe o começo de cada linha do `pedicularis.snps` com `cut -c1-40 pedicularis.snps`.

**Os comandos**

Os comandos abaixo terminam algumas linhas com `\`: é o jeito de quebrar um comando longo em várias linhas no terminal. Pode copiar assim mesmo, ou digitar tudo numa linha só, sem as barras.

## Análise 1: a matriz completa

A primeira árvore usa todos os sítios e serve de referência para as outras duas. O IQ-TREE busca a árvore de máxima verossimilhança e calcula o UFBoot com 1.000 réplicas (`-B 1000`):

    iqtree3 -s pedicularis.phy -st DNA -m GTR+I+G4 -B 1000 -T AUTO \
      -o 32082_przewalskii_SRR1754729,33588_przewalskii_SRR1754727 \
      --prefix completo --seed 2026

- `-st DNA`: avisa que os dados são DNA. Sem isso, o IQ-TREE para com o erro `cannot confidently detect the sequence data type`: a matriz tem tantos `N` (dados faltantes) e códigos de ambiguidade (os heterozigotos) que ele não tem certeza de que é DNA.

- `-m GTR+I+G4`: o modelo de substituição, fixado. Nas Análises 2 e 3, o IQ-TREE vai escolher o modelo sozinho, com o ModelFinder (Kalyaanamoorthy et al. 2017); aqui, com 2,8 milhões de colunas, a escolha levaria horas. Usamos o modelo mais geral, o GTR, com variação de taxas entre sítios (`+G4`) e uma proporção de sítios invariáveis (`+I`), já que quase todos os sítios da matriz são constantes.

- `-T AUTO`: o IQ-TREE escolhe quantos núcleos do processador usar.

- `--prefix`: o começo do nome de todos os arquivos de saída.

- `--seed`: a semente dos números aleatórios, para todos chegarem ao mesmo resultado.

Esta é a análise mais demorada da prática: deixe-a rodando e, enquanto isso, faça as Análises 2 e 3 em outra janela do terminal. Os principais arquivos de saída:

- `completo.iqtree`: o relatório, com o modelo escolhido, o logaritmo da verossimilhança, o comprimento total da árvore e a árvore desenhada em texto;

- `completo.treefile`: a árvore de máxima verossimilhança, com o UFBoot em cada nó;

- `completo.log`: tudo o que apareceu na tela.

Abra o `completo.treefile` no FigTree, enraíze nos dois *P. przewalskii* e mostre os rótulos dos nós (*Node Labels*).

**Questão 1.** No `completo.iqtree`, quantos sítios a matriz tem? Quantos são constantes, quantos são informativos para parcimônia (*parsimony-informative*) e quantos são *singletons*? Qual proporção de sítios invariáveis o modelo estimou (`Proportion of invariable sites`)? Compare com a proporção de sítios constantes.

**Questão 2.** As amostras da mesma espécie formam clados? Quantos nós têm UFBoot abaixo de 95 (o limite que a documentação do IQ-TREE sugere)?

## Análise 2: só os SNPs, com +ASC

Agora a matriz só de SNPs. O `+ASC` no fim do `-m MFP` faz o ModelFinder testar só modelos com a correção:

    iqtree3 -s pedicularis.snps -st DNA -m MFP+ASC -B 1000 -T AUTO \
      -o 32082_przewalskii_SRR1754729,33588_przewalskii_SRR1754727 \
      --prefix snps --seed 2026

**O erro esperado**

O IQ-TREE vai parar logo no começo, com esta mensagem:

```
ERROR: For your convenience alignment with variable sites printed to snps.varsites.phy
ERROR: Invalid use of +ASC because of 33197 invariant sites in the alignment
```

Mas o ipyrad não tinha tirado os sítios invariáveis? Tinha, só que “variável” não quer dizer a mesma coisa para os dois programas. Para o ipyrad, um sítio onde uma amostra é heterozigota (`R`, por exemplo, que é A ou G) já tem um SNP. Para o +ASC, a coluna só é variável se não puder ser explicada por uma única base: uma coluna só com `A`, `R` e `N` pode ser toda de A, e conta como invariável.

O IQ-TREE já resolve: junto com o erro, ele grava um alinhamento só com as colunas realmente variáveis, no arquivo `snps.varsites.phy` (o nome vem do `--prefix`). Rode de novo com ele:

    iqtree3 -s snps.varsites.phy -st DNA -m MFP+ASC -B 1000 -T AUTO \
      -o 32082_przewalskii_SRR1754729,33588_przewalskii_SRR1754727 \
      --prefix snps --seed 2026 -redo

O `-redo` é importante. A tentativa que deu erro deixou um arquivo de checkpoint (`snps.model.gz`), e, sem o `-redo`, o IQ-TREE tenta aproveitá-lo, mesmo ele sendo de outro alinhamento. Dependendo da versão, isso faz o programa travar no meio do ModelFinder, com um erro `Assertion ... failed` e `IQ-TREE CRASHES WITH SIGNAL ABORTED`. O `-redo` manda ignorar o checkpoint e começar do zero.

**Questão 3.** Quantas colunas o IQ-TREE tirou? Que fração das colunas do `pedicularis.snps` isso representa? Compare o número de colunas que sobraram com o número de sítios variáveis da Análise 1 (informativos mais *singletons*).

**Questão 4.** Qual modelo o ModelFinder escolheu, e por qual critério? Compare o comprimento total da árvore (no `.iqtree`) com o da Análise 1. Por que eles são tão diferentes, se os dados são os mesmos?

**Questão 5.** (Opcional.) Rode a Análise 2 *sem* a correção, com o modelo que o ModelFinder escolheu, mas sem o `+ASC` (para não esperar o ModelFinder de novo):

    iqtree3 -s snps.varsites.phy -st DNA -m TVM+G4 -B 1000 -T AUTO \
      -o 32082_przewalskii_SRR1754729,33588_przewalskii_SRR1754727 \
      --prefix snps_sem_asc --seed 2026

(Se o seu ModelFinder tiver escolhido outro modelo, use o seu, sem o `+ASC`.) O que acontece com o comprimento total da árvore? E com a topologia?

## Análise 3: um SNP por locus

Os SNPs de um mesmo locus estão a poucas dezenas de pares de bases uns dos outros e, quase sempre, compartilham a mesma história genealógica. Na matriz de todos os SNPs, um locus com dez SNPs pesa dez vezes mais que um locus com um. O `pedicularis.usnps` guarda um SNP sorteado por locus, mais perto da ideia de sítios independentes que o modelo assume:

    iqtree3 -s pedicularis.usnps -st DNA -m MFP+ASC -B 1000 -T AUTO \
      -o 32082_przewalskii_SRR1754729,33588_przewalskii_SRR1754727 \
      --prefix usnps --seed 2026

O mesmo erro do +ASC vai aparecer: rode de novo com o `usnps.varsites.phy`, como na Análise 2, sem esquecer o `-redo`:

    iqtree3 -s usnps.varsites.phy -st DNA -m MFP+ASC -B 1000 -T AUTO \
      -o 32082_przewalskii_SRR1754729,33588_przewalskii_SRR1754727 \
      --prefix usnps --seed 2026 -redo

**Questão 6.** Quantas colunas tem a matriz final, comparada à da Análise 2? O que aconteceu com os valores de UFBoot? A topologia mudou?

## Comparando as árvores

**As topologias**

A distância de Robinson-Foulds (RF) conta quantos clados estão numa árvore e não na outra; 0 quer dizer a mesma topologia. Junte as três árvores num arquivo e peça ao IQ-TREE as distâncias entre todos os pares:

    cat completo.treefile snps.treefile usnps.treefile > tres_arvores.tre
    iqtree3 -t tres_arvores.tre -rf_all

As distâncias ficam no `tres_arvores.tre.rfdist`, numa matriz em que as árvores aparecem na ordem do arquivo (`Tree0` é a da matriz completa, `Tree1` a dos SNPs e `Tree2` a de um SNP por locus).

**Os números**

Preencha com os relatórios (`.iqtree`) das três análises:

|                             | Completa | Todos os SNPs | Um SNP por locus |
|:----------------------------|:--------:|:-------------:|:----------------:|
| Colunas analisadas          |          |               |                  |
| Modelo                      | GTR+I+G4 |               |                  |
| Comprimento total da árvore |          |               |                  |
| Nós com UFBoot &lt; 95       |          |               |                  |
| Distância RF até a completa |    0     |               |                  |

**Questão 7.** Onde estão as diferenças entre as três árvores? Elas aparecem em nós com UFBoot alto ou baixo?

## O fator de concordância de sítio

Na árvore da matriz completa, vamos calcular o sCF de cada ramo. Para cada ramo, o IQ-TREE sorteia 100 quartetos de amostras em volta dele e, em cada quarteto, vê quais sítios apoiam o ramo e quais apoiam cada uma das duas alternativas. A árvore (`-te`) e o modelo são os da Análise 1:

    iqtree3 -te completo.treefile -s pedicularis.phy -st DNA -m GTR+I+G4 \
      --scfl 100 -T AUTO --prefix concord --seed 2026

- `concord.cf.tree`: a árvore com o UFBoot e o sCF em cada nó, separados por uma barra (por exemplo, `100/44.6`);

- `concord.cf.stat`: uma tabela com uma linha por ramo, com o sCF e os dois fatores de discordância, `sDF1` e `sDF2` (a porcentagem dos sítios que apoiam cada alternativa), e o número médio de sítios informativos por quarteto (`sN`);

- `concord.cf.branch`: a mesma árvore, com o número de cada ramo (a coluna `ID` do `concord.cf.stat`) no lugar dos suportes.

Abra o `concord.cf.tree` e o `concord.cf.branch` no FigTree, lado a lado, para saber a que ramo corresponde cada linha da tabela. (A árvore sai sem raiz: enraíze nos dois *P. przewalskii*, como antes.)

**Como ler**

Em volta de cada ramo, há três arranjos possíveis dos quatro grupos de amostras. Um sCF perto de 33% quer dizer que os sítios não preferem o ramo da árvore a nenhuma das alternativas; um sCF bem acima disso quer dizer que a maioria dos sítios concorda com ele. Um ramo pode ter UFBoot de 100 e sCF de 40%: é o mais apoiado dos três arranjos, mas com muito conflito. Os dois fatores de discordância também contam uma história: se sDF1 e sDF2 são parecidos, o conflito é o esperado pela separação incompleta de linhagens; se um é bem maior que o outro, um dos arranjos alternativos tem apoio demais, e isso pode indicar fluxo gênico entre as linhagens desse arranjo.

**Questão 8.** Quais ramos têm UFBoot alto (95 ou mais) e sCF abaixo de 50%? Eles estão entre espécies ou entre amostras da mesma espécie?

**Questão 9.** Para esses ramos, sDF1 e sDF2 são parecidos ou um é bem maior que o outro? O que isso sugere?

## Na vida real

**Um segundo suporte**

O UFBoot é rápido e, com modelos adequados, bem calibrado: um clado com 95 é verdadeiro em cerca de 95% dos casos. Mas ele pode superestimar o apoio quando o modelo está errado. Por isso, a documentação do IQ-TREE recomenda calcular junto o SH-aLRT (`-alrt 1000`) (Guindon et al. 2010), um teste que compara cada ramo com os dois arranjos alternativos em volta dele, reamostrando os sítios. Os dois saem no mesmo rótulo do nó (SH-aLRT/UFBoot), e um clado é considerado confiável quando o SH-aLRT é pelo menos 80 e o UFBoot pelo menos 95. Numa análise para publicação, rodaríamos sempre os dois.

**A escolha do modelo**

Na Análise 1, fixamos o GTR+I+G4 porque o ModelFinder levaria horas numa matriz de 2,8 milhões de colunas. Numa análise para publicação, escolheríamos o modelo mesmo assim, deixando o computador trabalhar ou limitando os modelos testados (por exemplo, só variação de taxas por gama, com `-mrate G,I+G`).

**Além da concatenação**

Uma árvore de máxima verossimilhança de dados concatenados assume que todos os loci têm a mesma história. Com RAD-seq de espécies próximas, isso raramente é verdade, e o sCF mostra o tamanho do conflito. Numa análise para publicação, complementaríamos com métodos de árvore de espécies que lidam com essa discordância, como o SVDquartets ou o tetrad (das ferramentas de análise do ipyrad). O fator de concordância de gene (gCF) exigiria uma árvore por locus, e os loci de RAD-seq (cerca de 70 pb aqui) são curtos demais para dar árvores informativas.

## Referências

- Eaton, D. A. R. & Ree, R. H. 2013. Inferring phylogeny and introgression using RADseq data: an example from flowering plants (*Pedicularis*: Orobanchaceae). *Systematic Biology* 62: 689–706.
- Guindon, S., Dufayard, J.-F., Lefort, V., Anisimova, M., Hordijk, W. & Gascuel, O. 2010. New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. *Systematic Biology* 59: 307–321.
- Hoang, D. T., Chernomor, O., von Haeseler, A., Minh, B. Q. & Vinh, L. S. 2018. UFBoot2: improving the ultrafast bootstrap approximation. *Molecular Biology and Evolution* 35: 518–522.
- Kalyaanamoorthy, S., Minh, B. Q., Wong, T. K. F., von Haeseler, A. & Jermiin, L. S. 2017. ModelFinder: fast model selection for accurate phylogenetic estimates. *Nature Methods* 14: 587–589.
- Lewis, P. O. 2001. A likelihood approach to estimating phylogeny from discrete morphological character data. *Systematic Biology* 50: 913–925.
- Minh, B. Q., Hahn, M. W. & Lanfear, R. 2020. New methods to calculate concordance factors for phylogenomic datasets. *Molecular Biology and Evolution* 37: 2727–2733.
- Mo, Y. K., Lanfear, R., Hahn, M. W. & Minh, B. Q. 2023. Updated site concordance factors minimize effects of homoplasy and taxon sampling. *Bioinformatics* 39: btac741.
- Wong, T. K. F., Ly-Trong, N., Ren, H. *et al.* 2025. IQ-TREE 3: phylogenomic inference software using complex evolutionary models. Preprint, EcoEvoRxiv, doi:10.32942/X2P62N.
