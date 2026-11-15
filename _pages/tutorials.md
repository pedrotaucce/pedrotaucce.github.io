---
layout: page
permalink: /tutorials/
title: Tutoriais da Pós
description: Tutoriais da disciplina de pós-graduação Métodos Computacionais em Sistemática e Filogeografia (Unesp Bauru).
nav: true
nav_order: 7
---

Material prático da disciplina [Métodos Computacionais em Sistemática e Filogeografia]({{ '/teaching/metodos_computacionais/' | relative_url }}), na ordem em que costumamos trabalhar: da obtenção das sequências à inferência filogenética.

{% assign tutoriais = site.tutorials | sort: "order" %}
<ol>
{% for t in tutoriais %}
  <li>
    <a href="{{ t.url | relative_url }}"><strong>{{ t.title }}</strong></a><br>
    <span class="text-muted">{{ t.description }}</span>
  </li>
{% endfor %}
</ol>
