---
layout: page
title: "Árvore de máxima verossimilhança com o IQ-TREE 3"
description: "Análise particionada no terminal"
description_en: "Partitioned maximum likelihood analysis with IQ-TREE 3 on the command line."
order: 3
toc:
  sidebar: left
---

Tutorial da disciplina de pós-graduação [Métodos Computacionais em Sistemática e Filogeografia]({{ '/teaching/metodos_computacionais/' | relative_url }}).

<div class="alert alert-secondary" role="alert">
<strong>Arquivos da prática:</strong> <a href="{{ '/assets/tutorials/iqtree3/pratica_iqtree.zip' | relative_url }}">pratica_iqtree.zip</a>. Individualmente: <a href="{{ '/assets/tutorials/iqtree3/hylodidae_concatenado.fasta' | relative_url }}">hylodidae_concatenado.fasta</a>, <a href="{{ '/assets/tutorials/iqtree3/particoes.txt' | relative_url }}">particoes.txt</a>.
</div>

## Máxima verossimilhança, partições e grupo externo

Na **máxima verossimilhança**, cada árvore possível recebe uma nota: a probabilidade de observarmos o nosso alinhamento se aquela árvore, com aqueles comprimentos de ramo, fosse a verdadeira. A árvore escolhida é a de maior verossimilhança. Para calcular essa probabilidade, é preciso um **modelo de substituição de nucleotídeos**, que descreve como as bases mudam ao longo do tempo. O IQ-TREE 3 (Wong et al. 2026) escolhe o modelo automaticamente com o ModelFinder (Kalyaanamoorthy et al. 2017) e, logo em seguida, procura a melhor árvore. O apoio a cada clado é medido pelo *bootstrap* ultrarrápido (Hoang et al. 2018).

Nossos dois marcadores não evoluem do mesmo jeito, por isso vamos dividir a matriz em **partições**, cada uma com o seu próprio modelo:

- o 16S, um gene de RNA ribossomal;

- a 1ª, a 2ª e a 3ª posições dos códons do COI, separadas porque mudanças na 3ª posição muitas vezes não alteram o aminoácido e se acumulam bem mais depressa que nas outras.

Separar as posições só é possível porque o alinhamento do COI começa na posição 1 do códon, como você conferiu no MEGA no tutorial [Alinhamento múltiplo de sequências]({{ '/tutorials/alinhamento/' | relative_url }}).

Por fim, a árvore de máxima verossimilhança não tem raiz. Para enraizá-la, indicamos um **grupo externo**: *Limnomedusa macroglossa*, que não pertence a Hylodidae.

## Instalação

Vamos baixar o IQ-TREE 3 (versão 3.1.4) direto do GitHub pelo terminal e descompactá-lo na sua pasta pessoal (`~`). Escolha os comandos do seu sistema.

**Windows (Ubuntu no WSL)**

```
$ cd ~
$ curl -LO \
  https://github.com/iqtree/iqtree3/releases/download/v3.1.4/iqtree-3.1.4-Linux-intel.tar.gz
$ tar -xzf iqtree-3.1.4-Linux-intel.tar.gz
```

**macOS**

```
$ cd ~
$ curl -LO \
  https://github.com/iqtree/iqtree3/releases/download/v3.1.4/iqtree-3.1.4-macOS.zip
$ unzip iqtree-3.1.4-macOS.zip
```

O `curl -LO` baixa o arquivo do endereço indicado, e o `tar -xzf` (ou o `unzip`) o descompacta. Agora existe, na sua pasta pessoal, a pasta `iqtree-3.1.4-Linux-intel` (ou `iqtree-3.1.4-macOS`), com o programa dentro de `bin`.

## Preparando a pasta de trabalho

O programa e os dados ficarão na mesma pasta. Vá até a pasta `pratica_iqtree` (veja o tutorial [Alinhamento múltiplo de sequências]({{ '/tutorials/alinhamento/' | relative_url }}): no Windows, ela fica dentro de `/mnt/c/Users/...`; no macOS, digite `cd`, um espaço, e arraste a pasta para o Terminal) e copie o programa para ela. O ponto no fim do comando significa “para a pasta atual”:

    $ cp ~/iqtree-3.1.4-Linux-intel/bin/iqtree3 .

No macOS, o caminho é `~/iqtree-3.1.4-macOS/bin/iqtree3`. Confira:

    $ ./iqtree3 --version

A primeira linha da resposta deve mostrar a versão: `IQ-TREE version 3.1.4`. O `./` diz ao terminal que o programa está na pasta atual.

## O arquivo de partições

Abra o `particoes.txt` do tutorial [Alinhamento múltiplo de sequências]({{ '/tutorials/alinhamento/' | relative_url }}) no Sublime Text ou no Notepad++. Ele diz quais colunas da matriz são do 16S e quais são do COI. Com o MAFFT v7.471, por exemplo:

    DNA, 16S = 1-2475
    DNA, COI = 2476-3187

Crie no editor um arquivo novo, escreva o texto abaixo e salve-o na pasta `pratica_iqtree` com o nome `particoes.nex`:

    #nexus
    begin sets;
      charset p16S     = 1-2475;
      charset COI_pos1 = 2476-3187\3;
      charset COI_pos2 = 2477-3187\3;
      charset COI_pos3 = 2478-3187\3;
    end;

Cada `charset` dá um nome a um conjunto de colunas. O `\3` significa “de 3 em 3”: `2476-3187\3` são as colunas 2476, 2479, 2482…, ou seja, todas as 1ª posições dos códons do COI; começando em 2477 temos as 2ª posições e, em 2478, as 3ª. O nome `p16S` começa com uma letra porque, no formato NEXUS, um nome não pode começar com número.

**Se os números do seu `particoes.txt` forem outros**, use os seus: o 16S vai de 1 até o último número dele, e as três partições do COI começam no primeiro número do COI, no seguinte e no outro, todas terminando no último número do COI. Não esqueça os pontos e vírgulas.

## Rodando a análise

Na pasta `pratica_iqtree`, digite:

    $ ./iqtree3 -s hylodidae_concatenado.fasta -p particoes.nex \
          -o Limnomedusa_macroglossa -B 1000 --prefix hylodidae

Como no tutorial [Alinhamento múltiplo de sequências]({{ '/tutorials/alinhamento/' | relative_url }}), a barra invertida indica que o comando continua na linha de baixo. Por partes:

- `-s hylodidae_concatenado.fasta`: o alinhamento;

- `-p particoes.nex`: as partições. O IQ-TREE escolhe um modelo para cada uma e estima uma única árvore para todas;

- `-o Limnomedusa_macroglossa`: o grupo externo, escrito exatamente como na matriz;

- `-B 1000`: *bootstrap* ultrarrápido com 1000 réplicas;

- `--prefix hylodidae`: o começo do nome de todos os arquivos de resultado.

A análise leva poucos segundos (cerca de 12 s em um computador comum). Enquanto ela roda, aparecem várias mensagens; algumas começam com `WARNING` (por exemplo, sobre sequências com muitos dados faltantes) e, nesta prática, podem ser ignoradas. No Ubuntu, também aparece no início a mensagem `libarcher.so: cannot open shared object file`, que não atrapalha nada. A análise terminou quando aparecer a lista de arquivos de resultado e o *prompt* voltar.

## Os resultados

O IQ-TREE cria vários arquivos começando com `hylodidae`. Os três principais são:

- `hylodidae.treefile`: a árvore de máxima verossimilhança, com os valores de *bootstrap* nos nós, no formato Newick;

- `hylodidae.iqtree`: o relatório da análise, em texto. Na seção `ModelFinder`, a linha `Best-fit model according to BIC` traz o modelo escolhido para cada partição;

- `hylodidae.log`: tudo o que apareceu na tela durante a análise.

**Vendo a árvore**

Use o FigTree (<https://github.com/rambaut/figtree/releases>; no Windows, ele precisa do Java instalado). No Windows, o FigTree é um programa comum, fora do Ubuntu, e enxerga a pasta `pratica_iqtree` normalmente.

1.  Abra o arquivo `hylodidae.treefile` (*File \> Open*). Quando o FigTree pedir um nome para os rótulos dos nós, digite `bootstrap` e clique em *OK*.

2.  Clique no ramo de *Limnomedusa macroglossa* e no botão *Reroot*, para que o grupo externo fique de um lado da raiz e Hylodidae do outro.

3.  No painel da esquerda, marque *Node Labels* e, em *Display*, escolha *bootstrap*.

Compare a sua árvore com a de um colega: os valores de *bootstrap* podem diferir um pouco, porque as réplicas são sorteadas.

## Perguntas

1.  Qual modelo foi escolhido para cada uma das quatro partições?

2.  Os gêneros *Crossodactylus*, *Hylodes* e *Phantasmarana* saíram monofiléticos? Com que valor de *bootstrap*?

## Problemas comuns {#ap-problemas}

- **`iqtree3: command not found`.** Faltou o `./` antes do nome, ou o programa não foi copiado para a pasta atual (confira com `ls`).

- **macOS: o sistema diz que não pode verificar o desenvolvedor.** Acontece quando o programa é baixado pelo navegador, e não pelo `curl`. Rode `xattr -d com.apple.quarantine iqtree3` na pasta `pratica_iqtree`, ou libere o programa em *Ajustes do Sistema \> Privacidade e Segurança*.

- **`ERROR: Too large site ID`.** Algum número do `particoes.nex` é maior que o número de colunas da matriz. Confira os números com o seu `particoes.txt`.

- **`ERROR: Alignment does not have specified outgroup taxon`.** O nome depois do `-o` está diferente do nome na matriz.

- **O IQ-TREE não encontra `particoes.nex`.** No Windows, o editor pode ter salvo o arquivo como `particoes.nex.txt`. Confira com `ls` e, se for o caso, renomeie com `mv particoes.nex.txt particoes.nex`.

- **`Use -redo option if you really want to redo the analysis`.** Já existe uma análise com o mesmo `--prefix`. Acrescente `-redo` ao comando para refazê-la, ou use outro prefixo.

## Referências

- Hoang, D. T., Chernomor, O., von Haeseler, A., Minh, B. Q. & Vinh, L. S. 2018. UFBoot2: improving the ultrafast bootstrap approximation. *Molecular Biology and Evolution* 35(2): 518–522.
- Kalyaanamoorthy, S., Minh, B. Q., Wong, T. K. F., von Haeseler, A. & Jermiin, L. S. 2017. ModelFinder: fast model selection for accurate phylogenetic estimates. *Nature Methods* 14: 587–589.
- Wong, T. K. F., Ly-Trong, N., Ren, H., Demotte, P., Banos, H., Roger, A. J., Susko, E., Bielow, C., De Maio, N., Goldman, N., Hahn, M. W., dos Reis, M., Vinh, L. S., Huttley, G., Lanfear, R. & Minh, B. Q. 2026. IQ-TREE 3: phylogenomic inference software using complex evolutionary models. *Molecular Biology and Evolution*, msag117.
