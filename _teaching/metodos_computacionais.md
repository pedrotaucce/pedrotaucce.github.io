---
layout: page
title: Métodos Computacionais em Sistemática e Filogeografia
description: "Disciplina de pós-graduação (Unesp), ofertada bienalmente"
img: /assets/img/teaching/coalescencia.png
importance: 1
category: Graduate
---

*Computational Methods in Systematics and Phylogeography: a biennial graduate course, taught in a hybrid format (in-person part at Unesp Bauru).*

Disciplina **bienal**, ofertada nos Programas de Pós-Graduação em **Biociências Interunidades** (campi de Bauru e Assis) e em **Biodiversidade** (campus de São José do Rio Preto) da Unesp. Tem **formato híbrido**, com a parte presencial na Unesp Bauru.

- **Professores:** Pedro Paulo Goulart Taucci e Marcus Thadeu Teixeira Santos
- **Carga horária:** 60 horas
- **Vagas:** 10

{% include figure.liquid path="assets/img/teaching/metodos_turma_2026.jpg" class="img-fluid rounded z-depth-1" alt="First class of the course / Primeira turma da disciplina em 2026" %}

*The first class of the course (2026), during the in-person classes at the School of Sciences, Unesp Bauru. / Primeira turma da disciplina (2026), durante as aulas presenciais na Faculdade de Ciências da Unesp Bauru.*

## Tutoriais

Os [tutoriais práticos da disciplina]({{ '/tutorials/' | relative_url }}) estão disponíveis no site, em português:

{% assign tutoriais = site.tutorials | sort: "order" %}
<ol>
{% for t in tutoriais %}
  <li><a href="{{ t.url | relative_url }}">{{ t.title }}</a></li>
{% endfor %}
</ol>

## Objetivo

O estudo da biodiversidade avança no mesmo passo da tecnologia, beneficiando-se de ferramentas computacionais em diversas áreas, como sistemática, filogeografia, biogeografia e genética da conservação. O aluno deverá concluir a disciplina compreendendo como aplicar o rico arcabouço de softwares disponíveis para o estudo da sistemática molecular em sua dissertação ou tese, trabalhando, quando possível, com seus próprios dados.

## Ementa

Noções de bioinformática: manipulação de arquivos, expressões regulares e ambiente UNIX; obtenção e alinhamento de sequências; seleção de modelos de substituição nucleotídica; inferência filogenética: máxima verossimilhança e inferência bayesiana; métodos de delimitação molecular de espécies; teoria coalescente e inferência de árvores de espécies; R aplicado à sistemática e filogeografia; filogeografia baseada em modelos; filogeografia comparada; introdução à genômica.

## Metodologia de ensino

Aulas expositivas ministradas com auxílio de software de apresentação de slides e projetor; discussão de artigos científicos; aulas práticas com o uso de computadores e tutoriais, que ficarão disponíveis online. Aulas teórico-práticas online e apresentação de trabalho presencial nos dois últimos dias da disciplina. É necessário que os alunos possuam computadores em que consigam fazer os tutoriais da disciplina.

## Forma de avaliação

Além da participação nas aulas e discussões, os alunos serão divididos em duplas e desenvolverão um trabalho ao longo da disciplina, utilizando algumas das técnicas e métodos aprendidos. Os dados podem ser próprios ou de conjuntos já publicados, aplicando, necessariamente, métodos diferentes. Ao final da disciplina, os alunos compartilharão os resultados na forma de uma pequena apresentação oral e um relatório em formato de artigo.

## Bibliografia básica

- Carstens, B.C., Pelletier, T.A., Reid, N.M., & Satler, J.D. (2013). How to fail at species delimitation. Molecular Ecology, 22, 4369-4383.
- Lemey, P., Salemi, M., & Vandamme, A.-M. (2009). The phylogenetic handbook: a practical approach to phylogenetic analysis and hypothesis testing. 2ª edição. University Press, New York, NY, USA.
- Paradis, E. (2020). Population Genomics with R (1st ed.). Chapman and Hall/CRC.
- Revell, L. J., & Harmon, L. J. (2022). Phylogenetic Comparative Methods in R: An authoritative introduction to the latest comparative methods in evolutionary biology. Princeton University Press, Princeton & Oxford, 430 pp.
- Scornavacca, C., Delsuc, F., & Galtier, N. (2020). Phylogenetics in the Genomic Era. No commercial publisher, authors open access book, 568 pp.
