---
layout: post
title: Construindo uma árvore de máxima verossimilhança com o software IQ-TREE3
date: 2025-05-28 17:00:00
description: A step-by-step beginner's guide 
tags: IQ-TREE2 MaximumLikelihood SubstitutionModel 
categories: Phylogenetics
thumbnail: assets/img/iqtree.png
---

Aproveitando a aula da disciplina de **Biogeografia**, hoje faremos um tutorial para construção de uma árvore de máxima verossimilhança e, para isso, vamos utilizar o programa <a href="https://https://iqtree.github.io/">IQ-TREE3</a>. Apesar de haver alguns software que são normalmente utilizados para inferência baseada em máxima verossimilhança, este é um dos que eu mais gosto pois é um programa rápido e bastante completo. Antes de começar, faça o download do programa no <a href="https://https://iqtree.github.io/">website</a>, de acordo com seu sistema operacional. Neste tutorial, usaremos Windows 11 e, se você usa o mesmo sistema e quer ir direto ao ponto, pode baixá-lo <a href="https://github.com/iqtree/iqtree3/releases/download/v3.0.1/iqtree-3.0.1-Windows.zip">aqui</a>.

Primeiramente, precisamos de um alinhamento. Vamos utilizar a matriz de um [fragmento de mtDNA](/assets/files/p_rupicola.fasta) (ou DNA mitocondrial) do gene 16S, do trabalho em que <a href="https://bioone.org/journals/journal-of-herpetology/volume-54/issue-2/19-114/A-New-Rupicolous-Species-of-the-Pristimantis-conspicillatus-Group-Anura/10.1670/19-114.short">Taucce et al. (2020)</a> descrevem o anuro *Pristimantis rupicola*. A matriz já está alinhada e pronta para ser utilizada. Agora, descompacte o arquivo que você baixou do site do **IQ-TREE2** e coloque a matriz *p_rupícola.fasta* dentro da pasta com os arquivos descompactados, em uma pasta com o programa propriamente dito, chamada *bin*.

Depois de tudo pronto, a primeira coisa que vamos fazer é clicar no programa, que está na pasta *bin*. Clique na versão do programa chamada **iqtree3-click.exe**. Você verá uma janela assim:

<div class="image-container mt-3">                                                                                          <figure id="fig1">                                                                                                          {% include figure.liquid loading="eager" path="assets/img/posts/2025-05-28/fig1.png" class="img-fluid rounded z-depth-1" zoomable=true %}
         <figcaption style="font-size: 0.9em; text-align: center;"><strong>Figura 1:</strong> Primeira janela do programa IQ-TREE3 na versão de clicar</figcaption>
    </figure>                                                                                                           </div>  

Digite *y* e aperte enter. Você verá uma mensagem pedindo para que possamos inserir os parâmetros. Para realizar uma análise de máxima verossimilhança, primeiramente precisaremos inferir o melhor modelo de substituição nucleotídica. No entanto, o IQ-TREE pode fazer isso de forma automática, logo antes da busca pela árvore mais verossímil. Isso salva um monte de tempo, então vamos usar! Faremos então a busca pelo modelo de substituição nucleotídica que mais se encaixa nos dados, depois a busca pela árvore mais verossímil seguida de 1000 réplicas de bootstrap, para que possamos saber o nível de suporte dos nossos clados. Para isso, vamos digitar um comando de cada vez. Primeiro, *-s*, e então, enter. Esse primeiro comando significa que vamos fazer a busca pelo melhor modelo e então pela melhor árvore. Aperte *e* de extend, enter e então escreva o nome do arquivo, que está na pasta, no nosso caso *p_rupicola.fasta*, e então enter de novo. Depois digite *e*, aperte enter e escreva *-B*, comando para as réplicas de bootstrap ultrarápido, implementado no IQ-TREE, e aperte enter. Depois de apertar *e* novamente, digite o número de réplicas de bootstrap: 1000 está de bom tamanho. Aperte enter e então sua tela deve estar assim:


<div class="image-container mt-3">
    <figure id="fig2">
        {% include figure.liquid loading="eager" path="assets/img/posts/2025-05-28/fig2.png" class="img-fluid rounded z-depth-1" zoomable=true %}
         <figcaption style="font-size: 0.9em; text-align: center;"><strong>Figura 2:</strong> Janela do IQ-TREE3 depois de todos os comandos digitados</figcaption>
    </figure>
</div>

Se tudo deu certo, confirme apertando *y* e então enter. Sua análise vai começar a rodar! Depois de tudo acabado, é só fechar a janela e dar uma olhada na pasta. Dentre os vários arquivos, está a nossa árvore: *p_rupicola.fasta.treefile*. Pronto! Agora é só usar um visualizador de árvore como por exemplo o <a href="https://github.com/rambaut/figtree/releases">FigTree</a> e comparar os resultados com a minha árvore!

<div class="image-container mt-3">
    <figure id="fig2">                                                                                                          {% include figure.liquid loading="eager" path="assets/img/posts/2025-05-28/fig3.png" class="img-fluid rounded z-depth-1" zoomable=true %}                                                                                                        <figcaption style="font-size: 0.9em; text-align: center;"><strong>Figura 3:</strong> Árvore de máxima verossimilhança seguida de 1000 réplicas de bootstrap inferida no programa IQ-TREE3</figcaption>
    </figure>
</div>
