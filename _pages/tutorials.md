---
layout: page
permalink: /tutorials/
title: Tutorials
description: Hands-on tutorials for the graduate course Computational Methods in Systematics and Phylogeography.
nav: true
nav_order: 7
---

These tutorials accompany the graduate course [Métodos Computacionais em Sistemática e Filogeografia]({{ '/teaching/metodos_computacionais/' | relative_url }}) (Computational Methods in Systematics and Phylogeography), taught at Unesp. They follow the usual workflow, from retrieving sequences to phylogenetic inference.

<div class="alert alert-info" role="alert">
<strong>Note:</strong> the tutorials are written in <strong>Brazilian Portuguese</strong>. Each entry below has a short English summary.
</div>

{% assign tutoriais = site.tutorials | sort: "order" %}
<ol>
{% for t in tutoriais %}
  <li>
    <a href="{{ t.url | relative_url }}"><strong>{{ t.title }}</strong></a><br>
    <span class="text-muted">{{ t.description_en }}</span>
  </li>
{% endfor %}
</ol>
