---
layout: page
permalink: /publications/
title: Publications
description: Peer-reviewed papers by LaBESTA members, in reverse chronological order.
nav: true
nav_order: 2
---

<!-- _pages/publications.md -->

<!-- Bibsearch Feature -->

{% include bib_search.liquid %}

<div class="publications">

<h3>Peer Reviewed papers</h3>
{% bibliography --query @article %}

<h3>Book Chapters</h3>
{% bibliography --query @incollection %}


</div>
