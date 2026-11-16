---
layout: page
title: "Inferência Bayesiana com MrBayes e delimitação de linhagens com mPTP"
description: "O ND2 de <em>Anolis chrysolepis</em> e <em>A. tandai</em>"
description_en: "Bayesian inference with MrBayes and lineage delimitation with mPTP."
order: 7
toc:
  sidebar: left
---

Tutorial da disciplina de pós-graduação [Métodos Computacionais em Sistemática e Filogeografia]({{ '/teaching/metodos_computacionais/' | relative_url }}).

<div class="alert alert-secondary" role="alert">
<strong>Arquivos da prática:</strong> <a href="{{ '/assets/tutorials/mrbayes-mptp/pratica_mrbayes_mptp.zip' | relative_url }}">pratica_mrbayes_mptp.zip</a>. Individualmente: <a href="{{ '/assets/tutorials/mrbayes-mptp/anolis_ND2.nex' | relative_url }}">anolis_ND2.nex</a>, <a href="{{ '/assets/tutorials/mrbayes-mptp/anolis_ND2.fas' | relative_url }}">anolis_ND2.fas</a>, <a href="{{ '/assets/tutorials/mrbayes-mptp/preparar_mptp.R' | relative_url }}">preparar_mptp.R</a>, <a href="{{ '/assets/tutorials/mrbayes-mptp/anolis_mptp.tre' | relative_url }}">anolis_mptp.tre</a>, <a href="{{ '/assets/tutorials/mrbayes-mptp/anolis_mptp.fas' | relative_url }}">anolis_mptp.fas</a>.
</div>

**Os dados**

Um subconjunto das sequências do gene mitocondrial ND2 de lagartos do grupo *Anolis chrysolepis*, publicadas por Vieira de Oliveira *et al.* (2025): 110 sequências de 924 pares de bases, sendo 58 de *A. chrysolepis* (nomes que começam com `CR_`), 50 de *A. tandai* (`TA_`) e uma de cada um dos dois grupos externos, *A. brasiliensis* (`BR_`) e *A. planiceps* (`PL_`).

**O roteiro**

Duas análises encadeadas, a segunda usando a saída da primeira:

1.  **MrBayes**: a filogenia Bayesiana do ND2 (Parte 1);

2.  **mPTP**: quantas linhagens a árvore sugere (Parte 2).

**Programas**

MrBayes 3.2.9 (Ronquist et al. 2012), Tracer (Rambaut et al. 2018), FigTree e mPTP 0.2.5 (Kapli et al. 2017). Opcional: R com o pacote `ape` (Paradis & Schliep 2019), para quem quiser preparar a própria árvore para o mPTP. As instruções de instalação estão em cada parte.

**Um aviso**

Esta é uma análise de tutorial, montada para rodar numa aula. Ao longo do texto, os quadros “**Na vida real**” indicam o que faríamos diferente numa análise para publicação.

## Introdução

### Inferência Bayesiana de filogenias

Na inferência Bayesiana, o que queremos é a probabilidade posterior de cada árvore (e de cada valor dos parâmetros do modelo), dados os dados. Pelo teorema de Bayes,

$$P(\text{árvore, parâmetros} \mid \text{dados}) =
\frac{P(\text{dados} \mid \text{árvore, parâmetros}) \times P(\text{árvore, parâmetros})}
{P(\text{dados})}.$$

No numerador estão a verossimilhança (a mesma da máxima verossimilhança) e o *prior*. O denominador, a verossimilhança marginal, é a soma de verossimilhança × *prior* sobre todas as árvores e todos os valores dos parâmetros. É por causa dele que não dá para calcular o posterior diretamente: para mais do que poucos terminais, há árvores demais. O MCMC (*Markov chain Monte Carlo*) contorna o problema porque só compara o posterior de um estado com o de outro. Nessa razão, P(dados) aparece em cima e embaixo e se cancela, e basta saber que o posterior é *proporcional* ao numerador.

**A analogia do andarilho**

Imagine todas as combinações possíveis de árvore e parâmetros como uma paisagem de morros, em que a altura de cada ponto é o posterior. Um andarilho começa num ponto qualquer e, a cada passo, sorteia um ponto vizinho: se o novo ponto é mais alto, ele vai; se é mais baixo, ele vai só às vezes, com probabilidade igual à razão das alturas (um ponto com metade da altura é aceito metade das vezes). Ele nunca precisa saber a altura de um ponto, só quantas vezes um é mais alto que o outro, e é por isso que P(dados) não faz falta. Depois de muitos passos, o tempo que o andarilho passa em cada região é proporcional à altura dela. Anotando onde ele está de tempos em tempos, temos uma amostra do posterior. A mesma imagem explica os três detalhes que vão aparecer na prática:

- ***Burn-in*.** O andarilho começa num lugar sorteado, em geral baixo, e os primeiros passos são a subida até os morros. Essas primeiras anotações não representam o posterior, e são descartadas (aqui, os primeiros 25%).

- **Cadeias aquecidas.** O risco é o andarilho ficar preso num morro e nunca cruzar o vale até outro morro, talvez mais alto. Por isso, cada corrida tem quatro andarilhos: um “frio”, que anda na paisagem verdadeira e é o único cujas posições anotamos, e três “aquecidos”, que andam numa versão achatada da mesma paisagem, com vales mais rasos, e cruzam de um morro para o outro com facilidade. De tempos em tempos, dois andarilhos propõem trocar de lugar; assim, o frio pode pular para um morro que um aquecido encontrou. Cada andarilho é uma *cadeia*; o conjunto é o MC<sup>3</sup> (*Metropolis-coupled* MCMC).

- **Duas corridas independentes.** Uma *corrida* é um grupo de quatro andarilhos. Rodamos duas, soltas em pontos diferentes da paisagem e sem contato uma com a outra. Se os andarilhos frios das duas acabam passando o mesmo tempo nos mesmos morros, temos um bom sinal de que os dois exploraram a paisagem toda: as corridas *convergiram*. Se cada um passou a maior parte do tempo num morro diferente, não convergiram. É isso que vamos medir (seção [A corrida convergiu?](#sec-convergencia)): o ASDSF compara as árvores das duas corridas, e o PSRF compara os parâmetros.

No fim, a probabilidade posterior de um clado é simplesmente a fração das árvores amostradas (depois do *burn-in*) em que aquele clado aparece.

### Delimitando linhagens com uma árvore: o mPTP

O modelo PTP (*Poisson Tree Processes*) (Zhang et al. 2013) parte de uma ideia simples: numa árvore de um único gene, os ramos *entre* espécies e os ramos *dentro* das espécies crescem de jeitos diferentes. Os ramos entre espécies refletem a especiação, um processo lento, com ramos longos; os ramos dentro das espécies refletem a coalescência dentro de cada população, com ramos curtos e muitos nós. O PTP modela o número de substituições nos ramos como dois processos exponenciais e procura o ponto da árvore em que se passa de um para o outro. Cada clado abaixo desse ponto é uma espécie (para nós, uma *linhagem*, veja abaixo).

O mPTP (*multi-rate* PTP) (Kapli et al. 2017) permite que cada espécie tenha a sua própria taxa do processo de dentro (a opção `--multi`), já que espécies com populações de tamanhos diferentes ou amostradas de jeitos diferentes têm ramos internos de tamanhos diferentes. O mPTP encontra a delimitação de máxima verossimilhança (`--ml`) e estima o suporte de cada nó com um MCMC (`--mcmc`).

Duas exigências do mPTP vão ditar a preparação da árvore:

- a árvore precisa ter comprimentos de ramo em substituições por sítio (uma árvore de máxima verossimilhança ou Bayesiana, e não um cronograma);

- a árvore precisa ser *binária* (sem politomias) e enraizada.

**Linhagem não é espécie**

Com um único gene mitocondrial, o mPTP delimita *linhagens mitocondriais*. Isso é uma hipótese a ser testada com outros dados (genes nucleares, morfologia, distribuição, canto, ecologia), e não uma lista de espécies. Métodos de um único loco tendem a dividir demais quando há estrutura populacional forte, e o DNA mitocondrial, herdado só pela mãe e sem recombinação, conta a história de uma única molécula.

## Parte 1: a filogenia Bayesiana com o MrBayes

### Instalação

Vamos usar a versão 3.2.9 do MrBayes. Na página da versão, <https://github.com/NBISweden/MrBayes/releases/tag/v3.2.9>, baixe o executável do seu sistema.

- **Windows.** O arquivo é o `mb_3.2.9_win_x64.exe`. **Copie esse arquivo para a pasta do tutorial** (a pasta onde está o `anolis_ND2.nex`) e abra-o de lá, com dois cliques: o MrBayes procura os arquivos na pasta em que ele está, e não vai achar o `anolis_ND2.nex` se ficar na pasta `Downloads`.

- **Mac e Linux.** Ponha o executável na pasta do tutorial e, no Terminal, dentro dela, deixe-o executável com `chmod +x` seguido do nome do arquivo. No Mac, rode antes `xattr -d com.apple.quarantine` seguido do nome do arquivo, para tirar o bloqueio que o Mac põe em programas baixados. Depois, rode-o com `./` seguido do nome (por exemplo, `./mb`, se você renomear o arquivo para `mb`).

- **WSL.** Quem quiser usar o WSL (o Linux dentro do Windows; instalação na seção [Instalando o mPTP](#sec-instalarmptp)) é só seguir as instruções do Linux, dentro do Ubuntu.

Vamos usar também o Tracer (<https://github.com/beast-dev/tracer/releases>) e o FigTree (<https://github.com/rambaut/figtree/releases>), que precisam do Java.

### O arquivo NEXUS

Abra o `anolis_ND2.nex` num editor de texto simples (Bloco de Notas, TextEdit em modo texto puro, ou o próprio RStudio). Um arquivo NEXUS é feito de blocos, cada um começando com `begin` e terminando com `end;`. O primeiro é o bloco de dados:

    #NEXUS
    begin data;
      dimensions ntax=110 nchar=924;
      format datatype=dna missing=? gap=-;
      matrix
    CR_SP.MRT6281    AACACAATAGCAATCATCCC...
    ...
    ;
    end;

O número de terminais (`ntax`) e de sítios (`nchar`) precisa bater com a matriz. Depois dela vem o bloco do MrBayes, com os comandos que o programa vai executar em sequência:

    begin mrbayes;
      set autoclose=yes nowarn=yes;
      log start filename=anolis_ND2_log.txt replace;
      outgroup BR_ACB.CHUNB14548;

      charset ND2_pos1 = 1-924\3;
      charset ND2_pos2 = 2-924\3;
      charset ND2_pos3 = 3-924\3;
      partition codons = 3: ND2_pos1, ND2_pos2, ND2_pos3;
      set partition = codons;

      lset applyto=(all) nst=mixed rates=invgamma;
      unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all);
      prset applyto=(all) ratepr=variable;

      mcmc ngen=100000 nruns=2 nchains=4 samplefreq=100 printfreq=1000
           diagnfreq=1000 relburnin=yes burninfrac=0.25;
      sump;
      sumt;
      sumt contype=allcompat conformat=simple outputname=anolis_ND2_allcompat;
      log stop;
    end;

Linha por linha:

- `set autoclose=yes nowarn=yes`: no fim da corrida, não pergunte se queremos continuar, e não peça confirmação para sobrescrever arquivos.

- `log start ...`: guarde tudo o que aparecer na tela no arquivo `anolis_ND2_log.txt`.

- `outgroup`: a sequência de *A. brasiliensis* enraíza as árvores que o MrBayes escreve (o MrBayes não estima a raiz; isso é só para desenhar).

- `charset` e `partition`: o ND2 é um gene codificador de proteína, e as três posições do códon evoluem de jeitos muito diferentes (a terceira, em geral, muito mais rápido). `1-924\3` quer dizer “do sítio 1 ao 924, de três em três”. As três posições viram três partições, cada uma com o seu modelo.

- `lset nst=mixed rates=invgamma`: o modelo de substituição de cada partição. Com `nst=mixed`, o MrBayes não usa um modelo fixo: ele *amostra* entre todos os modelos da família do GTR (do JC ao GTR) durante a própria corrida, em proporção à probabilidade posterior de cada um (Huelsenbeck et al. 2004). Assim, não precisamos escolher o modelo antes. `rates=invgamma` inclui a variação de taxas entre sítios (gama) e uma proporção de sítios invariáveis.

- `unlink`: por padrão, as partições compartilham os parâmetros; aqui, cada partição tem as suas próprias frequências de bases, taxas de substituição, forma da gama e proporção de invariáveis.

- `prset ratepr=variable`: cada partição tem a sua taxa média de evolução (a terceira posição pode evoluir mais rápido que as outras).

- `mcmc`: a corrida. `ngen`: número de gerações (os passos do andarilho), aqui 100 mil, para caber na aula; `nruns`: corridas independentes; `nchains`: cadeias (andarilhos) por corrida; `samplefreq`: amostra a cada 100 gerações (1.000 árvores por corrida, no total); `printfreq`: mostra o andamento na tela a cada 1.000 gerações; `diagnfreq`: calcula o diagnóstico de convergência a cada 1.000; `relburnin` e `burninfrac`: descarta os primeiros 25% das amostras.

- `sump`: resume os parâmetros (arquivos `.p`).

- `sumt`: resume as árvores (arquivos `.t`) numa árvore de consenso (seção [A árvore de consenso](#sec-consenso)).

- o segundo `sumt`: um consenso diferente, que só vai servir para o mPTP (seção [A árvore de consenso](#sec-consenso)).

**Questão 1.** Abra o MrBayes e digite `help`. Depois, `help lset` e `help mcmc`. Qual é o valor padrão de `nchains`? E o que faz a opção `temp`?

### Rodando

Abra o MrBayes (no Windows, com dois cliques no executável que você copiou para a pasta do tutorial; no Mac, no Linux e no WSL, com `./` e o nome do executável no Terminal, dentro da pasta) e digite:

    execute anolis_ND2.nex

O MrBayes lê os dados, monta o modelo e começa a corrida. A tela mostra o andamento a cada 1.000 gerações:

```
50000 -- [-7747.075] (-7790.338) (-7770.700) (-7820.862) * (-7773.493) ... -- 0:01:33
```

O primeiro número é a geração; depois vêm o logaritmo da verossimilhança das quatro cadeias da primeira corrida, um asterisco e as quatro da segunda. Os colchetes marcam a cadeia fria (o andarilho cujas posições anotamos); quando ela muda de lugar na linha, é porque dois andarilhos trocaram de posição. No fim da linha, o tempo que falta. A cada 1.000 gerações, aparece também

     Average standard deviation of split frequencies: 0.037924

que é o nosso primeiro diagnóstico de convergência (próxima seção).

A corrida de 100 mil gerações leva uns 3 minutos. Enquanto ela roda, leia a próxima seção.

### A corrida convergiu? {#sec-convergencia}

Três diagnósticos, do mais rápido ao mais completo:

1.  **O desvio padrão médio das frequências dos clados (ASDSF).** Para cada clado, o MrBayes compara a sua frequência nas duas corridas. Se as duas amostram as mesmas árvores, as frequências são parecidas, e a média dos desvios padrão tende a zero. O critério mais usado é ASDSF < 0,01; entre 0,01 e 0,05, a corrida pode servir para uma ideia geral, mas não para publicar.

2.  **O PSRF** (*potential scale reduction factor*) (Gelman & Rubin 1992). Na tabela do `sump`, compara a variação de cada parâmetro dentro e entre as corridas. Deve ficar perto de 1,00 (na prática, até 1,01 ou 1,02).

3.  **O ESS** (*effective sample size*). As amostras de uma cadeia são autocorrelacionadas: uma amostra se parece com a anterior. O ESS diz a quantas amostras independentes as nossas equivalem. Deve ser maior que 200 (alguns aceitam 100). O `sump` mostra o ESS mínimo (`minESS`) e o médio (`avgESS`) das duas corridas.

**Questão 2.** Quando a corrida terminar, procure no `anolis_ND2_log.txt` o último valor do ASDSF, e na tabela do `sump` o PSRF e o ESS de cada parâmetro. A corrida convergiu? Algum parâmetro ficou mal amostrado?

**O Tracer**

Abra o Tracer e carregue os dois arquivos de parâmetros, `anolis_ND2.nex.run1.p` e `anolis_ND2.nex.run2.p` (*File \> Import Trace File*). O Tracer descarta 10% como *burn-in* por padrão; para ficar igual ao MrBayes, troque para 25% clicando duas vezes no valor da coluna *Burn-In* de cada arquivo e digitando 25000 (o Tracer conta o *burn-in* em gerações, e não em amostras). Para cada parâmetro, veja:

- na aba *Trace*, o traço da cadeia ao longo das gerações: deve parecer uma “taturana peluda”, sem tendência de subir ou descer;

- na aba *Estimates*, o ESS (em vermelho quando é menor que 100, em amarelo entre 100 e 200);

- selecionando as duas corridas juntas na aba *Marginal Density*: as duas distribuições devem se sobrepor.

**Questão 3.** No Tracer, olhe o traço do `LnL` e de dois parâmetros com ESS baixo. Como eles se parecem? O que você faria para melhorar a amostragem?

> **Na vida real.** As 100 mil gerações daqui são só para a corrida caber na aula. Numa análise de verdade, rodaríamos pelo menos 10 milhões de gerações, e até todos os critérios serem atendidos, num computador mais rápido ou num servidor como o CIPRES (<https://www.phylo.org>). Também escolheríamos o esquema de partições e os modelos com um programa como o PartitionFinder ou o ModelFinder, em vez de partir direto para as posições do códon.

### A árvore de consenso {#sec-consenso}

O MrBayes amostrou cerca de 2.000 árvores (1.000 por corrida), das quais cerca de 1.500 ficam depois do *burn-in*. Para resumir essas árvores numa só, o `sumt` faz uma **árvore de consenso de maioria de 50%**: entram na árvore os clados que aparecem em mais da metade das árvores amostradas, e o número em cada nó é a probabilidade posterior do clado. Onde nenhum arranjo passa de 50%, a árvore fica com uma politomia. **É essa árvore, o consenso de 50%, que normalmente interpretamos e publicamos numa análise com o MrBayes**, porque ela só mostra os clados com mais apoio do que o contrário. Ela está no arquivo `anolis_ND2.nex.con.tre`.

O segundo `sumt`, com `contype=allcompat`, faz um consenso diferente: depois dos clados com mais de 50%, ele acrescenta, em ordem de frequência, os clados com menos de 50% que ainda forem compatíveis com os que já estão na árvore. O resultado é uma árvore quase totalmente resolvida. Ela *não* é para interpretar, porque mostra clados com pouco apoio como se fossem tão confiáveis quanto os outros; só a fazemos porque o mPTP exige uma árvore sem politomias (seção [Preparando a árvore](#sec-preparar)). Ela está em `anolis_ND2_allcompat.con.tre`.

Os comprimentos dos ramos nas duas árvores são as medianas do posterior, em substituições por sítio.

**O FigTree**

Abra o `anolis_ND2.nex.con.tre` no FigTree. Quando ele perguntar como chamar os rótulos dos nós, aceite o padrão. Depois:

- selecione o ramo de `BR_ACB.CHUNB14548` e clique em *Reroot* para enraizar no grupo externo;

- em *Tree \> Increasing Node Order*, ordene os nós para a árvore ficar mais fácil de ler;

- marque *Node Labels* no painel da esquerda e escolha *prob* em *Display*: é a probabilidade posterior de cada clado (de 0 a 1).

**Questão 4.** As sequências de *A. chrysolepis* (`CR_`) e de *A. tandai* (`TA_`) formam dois clados? Quantos clados bem apoiados (probabilidade posterior ≥ 0,95) você vê dentro de cada espécie? Há politomias? Onde?

## Parte 2: delimitando linhagens com o mPTP

### Preparando a árvore {#sec-preparar}

O mPTP precisa de uma árvore binária, enraizada e em formato Newick. O script `preparar_mptp.R` faz isso a partir do consenso `allcompat`:

1.  lê a árvore (`read.nexus`, do pacote `ape`);

2.  enraíza em *A. brasiliensis* (`root`);

3.  tira os dois grupos externos (`drop.tip`): eles são espécies diferentes, distantes das outras, e não queremos delimitá-los; o mPTP até tem opções para isso (`--outgroup` e `--outgroup_crop`), mas tirá-los antes é mais simples;

4.  se ainda houver alguma politomia, resolve-a ao acaso, com ramos novos de comprimento zero (`multi2di`). Na árvore de referência não há nenhuma, e o passo não muda nada. Mas o `allcompat` de uma corrida curta, como a de vocês, pode sair com politomias (quando os clados que as resolveriam nunca foram amostrados), e aí o mPTP se recusa a ler a árvore. Os ramos novos têm comprimento zero, menor que o ramo mínimo (próxima seção), e o mPTP não os conta; mas eles decidem quais sequências podem formar uma linhagem juntas, e por isso o script fixa a semente do sorteio;

5.  grava a árvore em Newick (`anolis_mptp.tre`) e um FASTA sem os grupos externos (`anolis_mptp.fas`), que o mPTP vai usar para calcular o ramo mínimo.

Para todos chegarem ao mesmo resultado, já deixamos prontos esses dois arquivos a partir de uma corrida de 3 milhões de gerações, e é com eles que vamos trabalhar. Se você quiser preparar a sua própria árvore, ponha o `preparar_mptp.R` na pasta da corrida, abra no RStudio, use *Session \> Set Working Directory \> To Source File Location* e rode. O script vai gravar por cima dos arquivos do professor; guarde uma cópia deles antes.

### Instalando o mPTP {#sec-instalarmptp}

O mPTP é um programa de linha de comando, sem janela: tudo é feito digitando comandos.

- **Windows.** Baixe o arquivo para Windows em <https://github.com/Pas-Kapli/mptp/releases>, descompacte, e copie o `mptp.exe` para a pasta da prática. O mPTP do Windows roda pelo **PowerShell**: abra a pasta da prática no Explorador de Arquivos, clique com o botão direito num espaço vazio e escolha *Abrir no Terminal* (ou, no Windows 10, segure Shift e escolha *Abrir janela do PowerShell aqui*). Nos comandos abaixo, troque `mptp` por `.\mptp.exe`; por exemplo, `.\mptp.exe --help`.

- **Linux.** Baixe o arquivo `.tar.gz` para Linux na mesma página, descompacte, e o programa é o `mptp` da pasta `bin`.

- **WSL.** Quem fizer a prática no WSL segue as instruções do Linux. Para instalar o WSL: no PowerShell como administrador, `wsl --install`, reinicie e abra o Ubuntu. As pastas do Windows ficam em `/mnt/c/` (por exemplo, `/mnt/c/Users/seu_nome/Downloads`).

- **Mac.** Baixe o `mptp-0.2.5-macos-aarch64.tar.gz` na mesma página (é a versão para os Macs com processador Apple, M1 ou mais novo), descompacte com dois cliques, e copie o `mptp` da pasta `bin` para a pasta da prática. No Terminal, dentro da pasta, rode uma vez:

      xattr -d com.apple.quarantine mptp
      chmod +x mptp

  O primeiro comando tira o bloqueio que o Mac põe em programas baixados; o segundo deixa o arquivo executável. Depois, chame o programa como `./mptp`.

Teste com `mptp --help` (`./mptp --help` no Mac e no Linux, se o programa estiver na pasta da prática; `.\mptp.exe --help` no Windows). Os comandos abaixo foram quebrados em mais de uma linha só para caber na página: no Terminal ou no PowerShell, digite cada um numa linha só.

### O comprimento mínimo de ramo

Sequências idênticas ou quase idênticas ficam ligadas por ramos muito curtos, que não refletem eventos de coalescência e atrapalham o modelo. O mPTP ignora os ramos mais curtos que um valor mínimo, e esse valor pode ser calculado a partir do alinhamento:

    mptp --minbr_auto anolis_mptp.fas --tree_file anolis_mptp.tre
         --output_file anolis_minbr

(tudo numa linha só). O valor aparece na tela, numa linha que começa com `Minimum branch length (--minbr) should be set to`, e no arquivo `anolis_minbr.txt`. Use o valor do arquivo, com seis casas decimais, para todos chegarem aos mesmos números.

**Questão 5.** Qual é o comprimento mínimo de ramo sugerido? A quantas substituições no ND2 inteiro (924 sítios) ele corresponde, aproximadamente?

### A delimitação de máxima verossimilhança

Com o valor do ramo mínimo (troque `X` pelo valor da seção anterior):

    mptp --ml --multi --tree_file anolis_mptp.tre --minbr X
         --output_file anolis_ml

Na tela, e no arquivo `anolis_ml.txt`:

- `Null-model score` e `Best score for multi coalescent rate`: os logaritmos da verossimilhança do modelo nulo (todas as sequências são uma espécie só) e da melhor delimitação;

- `LRT computed p-value`: o teste da razão de verossimilhanças entre os dois; se `LRT: passed`, a delimitação é significativamente melhor que o modelo nulo;

- `Number of delimited species`: o número de linhagens;

- depois, os blocos `Species 1:`, `Species 2:`, …, com as sequências de cada linhagem.

O arquivo `anolis_ml.svg` mostra a árvore com as linhagens: abra no navegador (Chrome, Firefox, Safari). Nele, os ramos do processo entre espécies aparecem em verde, e os ramos de dentro de cada linhagem, em vermelho: cada grupo vermelho é uma linhagem.

**Questão 6.** Quantas linhagens o mPTP encontrou? Alguma linhagem mistura sequências de *A. chrysolepis* e de *A. tandai*? Alguma tem uma sequência só?

### O suporte da delimitação

A delimitação de máxima verossimilhança é uma só, mas outras delimitações podem ser quase tão boas. Para saber quanto confiar em cada linhagem, o mPTP faz um MCMC sobre as delimitações possíveis, nesta mesma árvore. Vamos rodar duas cadeias, de 10 milhões de gerações cada, amostrando a cada 10.000 e descartando o primeiro milhão:

    mptp --mcmc 10000000 --multi --mcmc_runs 2 --mcmc_sample 10000
         --mcmc_burnin 1000000 --tree_file anolis_mptp.tre --minbr X
         --output_file anolis_mcmc --seed 2026

Leva poucos segundos. O MCMC do mPTP é bem mais simples que o do MrBayes, porque a árvore está fixa: só a delimitação muda. Na tela, no fim:

- `CCI` e `HPD`: o intervalo de credibilidade do número de linhagens;

- `ML average support based on run with seed ...`: o suporte médio das linhagens da delimitação de máxima verossimilhança, em cada corrida;

- `Average standard deviation of support values among runs`: a diferença entre as duas corridas (deve ser perto de zero, como o ASDSF do MrBayes).

O arquivo `anolis_mcmc.2026.combined.svg` mostra a árvore com o suporte das duas corridas juntas: a cor vai do preto (suporte 1, nó entre linhagens) ao vermelho (suporte 0, nó dentro de uma linhagem), e o valor aparece escrito nos nós com suporte acima de 0,5. Os valores de todos os nós estão no `anolis_mcmc.2026.combined.tree`, que abre no FigTree (escolha *label* em *Node Labels*). O suporte de cada nó é a fração das amostras do MCMC em que aquele nó estava no processo *entre* espécies, ou seja, a fração em que ele separava linhagens diferentes. Para cada linhagem, olhe dois nós:

- o nó logo *acima* da linhagem (o ancestral dela com a linhagem vizinha): perto de 1 quer dizer que ela foi separada das vizinhas em quase todas as amostras;

- o nó da *base* da linhagem (o ancestral comum das sequências dela): perto de 0 quer dizer que ela ficou inteira em quase todas as amostras; um valor alto quer dizer que, em muitas amostras, o mPTP a dividiu em mais de uma.

A mesma informação, resumida para o número de linhagens, está nos arquivos `.stats`: a porcentagem das amostras com 14, 15, 16… linhagens.

**Questão 7.** Qual é o suporte médio da delimitação de máxima verossimilhança? Pelo arquivo `.stats`, o número de linhagens da delimitação de máxima verossimilhança é também o mais frequente no MCMC? Quais linhagens o MCMC costuma dividir?

> **Na vida real.** O resultado depende da árvore: valeria repetir o mPTP com a árvore de máxima verossimilhança (de um IQ-TREE ou RAxML) e com árvores amostradas do posterior, e comparar. Também valeria comparar com outros métodos de delimitação (ABGD, ASAP, GMYC), e, sobretudo, testar as linhagens com genes nucleares e outros dados antes de chamar qualquer uma de espécie.

## Referências

- Gelman, A. & Rubin, D. B. 1992. Inference from iterative simulation using multiple sequences. *Statistical Science* 7: 457–472.
- Huelsenbeck, J. P., Larget, B. & Alfaro, M. E. 2004. Bayesian phylogenetic model selection using reversible jump Markov chain Monte Carlo. *Molecular Biology and Evolution* 21: 1123–1133.
- Kapli, P., Lutteropp, S., Zhang, J., Kobert, K., Pavlidis, P., Stamatakis, A. & Flouri, T. 2017. Multi-rate Poisson tree processes for single-locus species delimitation under maximum likelihood and Markov chain Monte Carlo. *Bioinformatics* 33: 1630–1638.
- Paradis, E. & Schliep, K. 2019. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics* 35: 526–528.
- Rambaut, A., Drummond, A. J., Xie, D., Baele, G. & Suchard, M. A. 2018. Posterior summarization in Bayesian phylogenetics using Tracer 1.7. *Systematic Biology* 67: 901–904.
- Ronquist, F., Teslenko, M., van der Mark, P., Ayres, D. L., Darling, A., Höhna, S., Larget, B., Liu, L., Suchard, M. A. & Huelsenbeck, J. P. 2012. MrBayes 3.2: efficient Bayesian phylogenetic inference and model choice across a large model space. *Systematic Biology* 61: 539–542.
- Vieira de Oliveira, A. P., D’Angiolella, A., Azevedo, J. A. R., Vitt, L. J., Colli, G. R., Choueri, E. L., Gamble, T., Taucce, P. P. G., Rodrigues, M. T., Noonan, B. P., Avila-Pires, T. C. & Werneck, F. P. 2025. The influence of historical climate cycles on the diversification and demography of species from the Neotropical lizard group *Anolis chrysolepis*. *Biological Journal of the Linnean Society* 144: blaf006.
- Zhang, J., Kapli, P., Pavlidis, P. & Stamatakis, A. 2013. A general species delimitation method with applications to phylogenetic placements. *Bioinformatics* 29: 2869–2876.
