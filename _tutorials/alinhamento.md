---
layout: page
title: "Alinhamento múltiplo de sequências"
description: "MUSCLE no MEGA, MAFFT no terminal e concatenação no R"
description_en: "Multiple sequence alignment with MUSCLE (MEGA), MAFFT (command line) and concatenation in R."
order: 2
toc:
  sidebar: left
---

Tutorial da disciplina de pós-graduação [Métodos Computacionais em Sistemática e Filogeografia]({{ '/teaching/metodos_computacionais/' | relative_url }}).

<div class="alert alert-secondary" role="alert">
<strong>Arquivos da prática:</strong> <a href="{{ '/assets/tutorials/alinhamento/pratica_alinhamento.zip' | relative_url }}">pratica_alinhamento.zip</a>. Individualmente: <a href="{{ '/assets/tutorials/alinhamento/04_concatenar.R' | relative_url }}">04_concatenar.R</a>, <a href="{{ '/assets/tutorials/alinhamento/16S_ids.fasta' | relative_url }}">16S_ids.fasta</a>, <a href="{{ '/assets/tutorials/alinhamento/COI_ids.fasta' | relative_url }}">COI_ids.fasta</a>.
</div>

## Alinhamento e homologia primária

Apesar de muitos sistematas darem pouca importância ao alinhamento e o tratarem apenas como parte do processo para a obtenção de uma árvore filogenética, esta etapa é extremamente importante para a análise filogenética como um todo, pois é nela que definimos a **hipótese de homologia primária** entre as nossas sequências. Em outras palavras, nesta etapa assumimos que as bases alinhadas em uma mesma coluna **divergiram de um mesmo estado ancestral**.

Sendo esta uma etapa tão importante, existem vários programas e algoritmos comumente utilizados para alinhar sequências de DNA (e também de aminoácidos, que não abordaremos no curso), como o ClustalW (Larkin et al. 2007), o T-Coffee (Notredame et al. 2000), o MUSCLE (Edgar 2004) e o MAFFT (Katoh & Standley 2013). Nesta prática usaremos os dois últimos.

Continuamos o trabalho com a filogenia de Hylodidae do tutorial [Aquisição de sequências no GenBank]({{ '/tutorials/genbank/' | relative_url }}). Agora que você adquiriu as sequências de 16S rRNA e COI no GenBank, precisa alinhá-las. Usaremos as matrizes baixadas por número de acesso, menores e mais controladas, com um exemplar por espécie. Elas incluem *Limnomedusa macroglossa* (Alsodidae), que servirá como raiz da árvore nas próximas análises.

## MUSCLE no MEGA 12 {#sec-mega}

O MEGA (Kumar et al. 2024) traz o MUSCLE embutido. Vamos usá-lo para alinhar o COI:

1.  No MEGA 12, vá em *Align \> Edit/Build Alignment*, escolha *Retrieve sequences from a file* e abra `sequencias/COI_ids.fasta`. Está tudo desalinhado.

2.  Selecione todas as sequências (*Edit \> Select All*, ou **Ctrl+A**) e vá em *Alignment \> Align by MUSCLE*. Se o MEGA oferecer as opções *DNA* e *Codons*, escolha *DNA*.

3.  Na janela de parâmetros do MUSCLE, mantenha os valores padrão e clique em *Compute*. Aguarde um pouquinho e pronto: as sequências estão alinhadas.

O COI é um gene codificante, então é interessante verificar se as sequências estão corretas: se começam no códon certo e se os possíveis *gaps* são múltiplos de três, por exemplo. Isso influenciará etapas posteriores, por isso é interessante que o alinhamento comece sempre na posição 1 do códon. A tradução também serve para checar se há *stop codons* no meio das sequências, o que poderia indicar que você sequenciou um pseudogene ou o gene errado.

1.  Clique na aba *Translated Protein Sequences*. Quando aparecer a pergunta *The Current Genetic Code is: Standard. Is this correct?*, clique em *No*.

2.  Na janela seguinte, como se trata de um fragmento mitocondrial de anfíbios, marque *Vertebrate Mitochondrial* e confirme.

3.  Verifique se há *stop codons* (indicados pelo símbolo `*`) em alguma posição. Se houver, provavelmente o alinhamento não começa na posição 1 do códon, ou o gene sequenciado não é o que você queria.

## MAFFT no terminal {#sec-mafft}

O MAFFT é um programa um pouco mais completo que o MUSCLE, por ter vários algoritmos diferentes (*estratégias*) para alinhar as sequências. Alguns são mais rápidos e menos precisos, ideais para um grande número de sequências ou para fragmentos codificantes, sem estrutura secundária; outros são mais precisos e mais lentos, ideais, por exemplo, para fragmentos com estrutura secundária, como o 16S. Por serem mais lentos, estes últimos não podem ser utilizados em grandes conjuntos de sequências. As duas estratégias que usaremos são:

| Estratégia                | Opções no comando               | Uso nesta prática |
|:--------------------------|:--------------------------------|:------------------|
| FFT-NS-2 (rápida, padrão) | `--retree 2`                    | COI               |
| E-INS-i (precisa, lenta)  | `--genafpair --maxiterate 1000` | 16S               |

A FFT-NS-2 parte de uma matriz de distâncias confiável, mas não é tão lenta quanto outras estratégias do MAFFT. A E-INS-i lida bem com fragmentos como o 16S, em que regiões conservadas se alternam com regiões cheias de inserções e deleções; ela faz iterações, refinando o alinhamento até chegar a um resultado satisfatório, e `--maxiterate 1000` aumenta o número de iterações possíveis para 1000. Em todos os comandos usaremos também `--preservecase`, para que as bases sejam mantidas em letra maiúscula. Para saber mais sobre as estratégias, veja <https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html>.

### Instalação

**Windows (Ubuntu no WSL)**

No terminal do Ubuntu, digite os comandos abaixo; o `sudo` pede a senha do Ubuntu (que não aparece enquanto é digitada) e o `apt` pergunta se deve continuar (responda `Y`):

    $ sudo apt update
    $ sudo apt install mafft

**macOS**

Na página do MAFFT (<https://mafft.cbrc.jp/alignment/software/>), entre em *Mac OS X*, baixe o pacote padrão (`.pkg`), abra-o e siga as instruções do instalador. Quem usa o Homebrew pode, em vez disso, digitar `brew install mafft`.

**Conferindo**

Feche e abra o terminal e digite:

    $ mafft --version

Deve aparecer o número da versão (por exemplo, `v7.505`, seguido da data). Lembre-se de que o símbolo `$` representa o *prompt* e não deve ser digitado; no seu terminal ele pode ser outro, como `%`.

### Indo até a pasta da prática

O MAFFT precisa encontrar os arquivos, então primeiro vamos até a pasta `pratica_alinhamento` com o `cd` (prática 2 da disciplina).

- **Windows:** os arquivos do Windows ficam, no Ubuntu, dentro de `/mnt/c`. Se a pasta estiver em *Downloads*, por exemplo, o comando é o abaixo, trocando `Nome` pelo seu usuário do Windows (a tecla **Tab** completa os nomes):

      $ cd /mnt/c/Users/Nome/Downloads/pratica_alinhamento

  Outra opção: abra a pasta no Explorador de Arquivos, clique na barra de endereço, digite `wsl` e tecle **Enter**; o Ubuntu abre já dentro da pasta.

- **macOS:** digite `cd` seguido de um espaço, arraste a pasta do Finder para a janela do Terminal e tecle **Enter**.

Confira se você está no lugar certo: `pwd` mostra o caminho da pasta atual, e a lista mostrada por `ls` deve incluir as pastas `R` e `sequencias`.

    $ pwd
    $ ls

Agora crie a pasta onde ficarão os alinhamentos:

    $ mkdir alinhamentos

### Alinhando um marcador

Vamos começar pelo COI, com a estratégia FFT-NS-2:

    $ mafft --preservecase --retree 2 \
          sequencias/COI_ids.fasta > alinhamentos/COI_alinhado.fasta

A barra invertida `\` no fim da primeira linha indica que o comando continua na linha de baixo: tecle **Enter** logo depois dela e continue digitando (você também pode escrever tudo em uma única linha, sem a barra). Por partes:

- `mafft --preservecase --retree 2`: o programa e as opções (seção [MAFFT no terminal](#sec-mafft));

- `sequencias/COI_ids.fasta`: o arquivo de entrada, com as sequências desalinhadas;

- `> alinhamentos/COI_alinhado.fasta`: o sinal `>` redireciona a saída do MAFFT, isto é, o alinhamento, para um arquivo, como no prática 2 da disciplina. **Não use o mesmo nome do arquivo de entrada**, ou ele será apagado antes de ser lido!

Enquanto trabalha, o MAFFT mostra mensagens de progresso na tela; elas não vão para o arquivo. Quando o *prompt* voltar, confira com `ls alinhamentos` se o arquivo foi criado.

### Alinhando os dois marcadores de uma vez {#sec-laco}

Para o 16S, bastaria repetir o comando trocando o nome dos arquivos e a estratégia. Com muitos marcadores, porém, repetir comandos dá trabalho e abre espaço para erros. Um **laço** (*loop*) executa os mesmos comandos para cada item de uma lista. Digite as linhas abaixo, teclando **Enter** ao fim de cada uma:

    for marcador in 16S COI
    do
      if [ "$marcador" = "16S" ]
      then
        mafft --preservecase --genafpair --maxiterate 1000 \
          sequencias/${marcador}_ids.fasta > alinhamentos/${marcador}_alinhado.fasta
      else
        mafft --preservecase --retree 2 \
          sequencias/${marcador}_ids.fasta > alinhamentos/${marcador}_alinhado.fasta
      fi
    done

Depois da primeira linha, o terminal troca o *prompt* por um sinal de continuação (`>` no bash do Ubuntu, `for>` no zsh do macOS), indicando que espera o resto do laço. Os comandos só são executados depois que você digita `done`. Os espaços no início das linhas são opcionais: servem apenas para facilitar a leitura.

Como o laço funciona:

O alinhamento do 16S leva mais tempo; espere o *prompt* voltar e confira:

    $ ls alinhamentos
    16S_alinhado.fasta  COI_alinhado.fasta

O `COI_alinhado.fasta` foi refeito pelo laço, com o mesmo resultado do comando anterior.

Abra os dois alinhamentos no MEGA (*Align \> Edit/Build Alignment \> Retrieve sequences from a file*) ou em um editor de texto e compare-os: **qual a principal diferença entre o alinhamento do 16S e o do COI? Qual dos dois demorou mais para ser concluído? Por quê?** Para medir o tempo, repita um dos comandos com `time` na frente (por exemplo, `time mafft --preservecase --retree 2 ...`).

## Concatenando os alinhamentos no R {#sec-concatenar}

Nas próximas práticas, os dois marcadores serão analisados juntos, em uma única matriz: o 16S nas primeiras colunas e o COI em seguida. Essa união se chama **concatenação**, e deve ser feita **pelo nome da espécie**: a linha de *Hylodes nasus* da matriz final deve juntar o 16S e o COI de *H. nasus*, embora os dois fragmentos tenham números de acesso diferentes. Quando uma espécie tem apenas um dos marcadores, o outro é preenchido com `?`, o símbolo de **dado faltante**.

Abra o RStudio, crie um projeto na pasta `pratica_alinhamento` (*File \> New Project \> Existing Directory*) e abra o script `R/04_concatenar.R`. Execute-o linha a linha (**Ctrl+Enter**, ou **Cmd+Enter** no macOS), acompanhando os resultados no console. O script usa o pacote ape (Paradis & Schliep 2019), instalado no tutorial [Aquisição de sequências no GenBank]({{ '/tutorials/genbank/' | relative_url }}). Os principais passos são descritos abaixo.

**1. Ler os alinhamentos**

`read.FASTA()` lê o arquivo e `as.matrix()` transforma o alinhamento em uma matriz, com uma linha por sequência e uma coluna por posição; `dim()` mostra o número de linhas e de colunas:

    s16 <- as.matrix(read.FASTA("alinhamentos/16S_alinhado.fasta"))
    coi <- as.matrix(read.FASTA("alinhamentos/COI_alinhado.fasta"))
    dim(s16)
    dim(coi)

**2. Separar a espécie do número de acesso**

Duas funções com expressões regulares (prática 1 da disciplina) separam as partes do nome:

    especie <- function(x) sub("_[A-Z]+_?[0-9]+$", "", x)
    acesso  <- function(x) sub("^.*_([A-Z]+_?[0-9]+)$", "\\1", x)

Assim, `Crossodactylus_schmidti_AY843579` vira `Crossodactylus_schmidti` e `AY843579`. O script monta uma tabela com o acesso usado para cada espécie em cada marcador e, em seguida, renomeia as linhas das duas matrizes apenas com o nome da espécie.

**3. Completar com dados faltantes**

`union()` reúne os nomes das duas matrizes, e `setdiff()` mostra quais espécies faltam em cada uma:

    nomes <- sort(union(rownames(s16), rownames(coi)))
    setdiff(nomes, rownames(s16))   # espécies sem 16S
    setdiff(nomes, rownames(coi))   # espécies sem COI

A função `completar()` acrescenta a cada matriz uma linha de `?` para cada espécie que falta e põe as linhas na mesma ordem (a de `nomes`).

**4. Concatenar**

Com as duas matrizes com as mesmas linhas, na mesma ordem, `cbind()` junta as colunas:

    concatenada <- cbind(s16_completo, coi_completo)
    dim(concatenada)

O número de colunas da matriz concatenada deve ser a soma das colunas do 16S e do COI.

**5. Gravar**

O script cria a pasta `matriz` com:

- `hylodidae_concatenado.fasta`, `.phy` e `.nex`: a mesma matriz nos formatos FASTA, PHYLIP e NEXUS, pedidos por diferentes programas de análise filogenética;

- `particoes.txt`: as colunas ocupadas por cada fragmento, que serão usadas para analisar cada marcador com seu próprio modelo de evolução:

      DNA, 16S = 1-2475
      DNA, COI = 2476-3187

- `tabela_acessos.csv`: os números de acesso usados para cada espécie, a tabela que acompanha a publicação de um trabalho.

Os números de colunas citados acima são os obtidos com o MAFFT v7.471; com outra versão, eles podem variar um pouco.

## Para refletir

1.  No MEGA, por que foi preciso trocar o código genético para *Vertebrate Mitochondrial*? O que aconteceria com o código *Standard*?

2.  O que seria preciso mudar no laço da seção [Alinhando os dois marcadores de uma vez](#sec-laco) para alinhar também um terceiro marcador, por exemplo o citocromo *b*, em um arquivo `cytb_ids.fasta`, com a estratégia FFT-NS-2?

3.  Quantas espécies há na matriz concatenada? Quais delas não têm o 16S e quais não têm o COI? Quantas colunas tem a matriz?

4.  Por que a concatenação é feita pelo nome da espécie, e não pela posição das sequências nos arquivos ou pelo número de acesso?

## Problemas comuns {#ap-problemas}

- **`mafft: command not found`.** O MAFFT não foi instalado, ou o terminal foi aberto antes da instalação: feche-o e abra-o de novo. No macOS, com o pacote `.pkg`, o programa fica em `/usr/local/bin/mafft`.

- **`No such file or directory`.** Você não está na pasta `pratica_alinhamento`, ou o nome do arquivo foi digitado errado. Confira com `pwd` e `ls`.

- **Arquivo alinhado vazio (0 bytes).** O MAFFT não conseguiu ler a entrada (por exemplo, porque o nome do arquivo de entrada está errado); leia a mensagem de erro na tela. Se o erro for sobre a pasta `alinhamentos`, ela ainda não existe: rode `mkdir alinhamentos`.

- **O terminal fica mostrando `>` e não executa nada.** O laço ou o `if` não foi fechado: falta digitar `fi` ou `done`. Se errar no meio, tecle **Ctrl+C** e comece de novo.

- **`cannot open file` no R.** A pasta de trabalho não é `pratica_alinhamento`, ou os alinhamentos ainda não foram feitos. Confira com `getwd()` e `list.files()`.

- **`there is no package called ape`.** Instale o pacote como no tutorial [Aquisição de sequências no GenBank]({{ '/tutorials/genbank/' | relative_url }}) (`install.packages("ape")`).

## Referências

- Edgar, R. C. 2004. MUSCLE: multiple sequence alignment with high accuracy and high throughput. *Nucleic Acids Research* 32(5): 1792–1797.
- Katoh, K. & Standley, D. M. 2013. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution* 30(4): 772–780.
- Kumar, S., Stecher, G., Suleski, M., Sanderford, M., Sharma, S. & Tamura, K. 2024. MEGA12: Molecular Evolutionary Genetic Analysis version 12 for adaptive and green computing. *Molecular Biology and Evolution* 41(12): msae263.
- Larkin, M. A., Blackshields, G., Brown, N. P., Chenna, R., McGettigan, P. A., McWilliam, H., Valentin, F., Wallace, I. M., Wilm, A., Lopez, R., Thompson, J. D., Gibson, T. J. & Higgins, D. G. 2007. Clustal W and Clustal X version 2.0. *Bioinformatics* 23(21): 2947–2948.
- Notredame, C., Higgins, D. G. & Heringa, J. 2000. T-Coffee: a novel method for fast and accurate multiple sequence alignment. *Journal of Molecular Biology* 302(1): 205–217.
- Paradis, E. & Schliep, K. 2019. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics* 35: 526–528.
