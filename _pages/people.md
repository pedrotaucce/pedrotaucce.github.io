---
layout: page
permalink: /people/
title: People
description: Members of LaBESTA
nav: true
nav_order: 8
---

<div class="people-layout">
  <nav class="people-nav" aria-label="People sections">
    {% for group in site.data.people.groups %}
    <a href="#{{ group.id }}">{{ group.title }}</a>
    {% endfor %}
    <a href="#join">Join the lab</a>
  </nav>

  <div class="people-content">
    {% for group in site.data.people.groups %}
    <section id="{{ group.id }}">
      <h2 class="category">{{ group.title }}</h2>
      <div class="people-grid{% if group.wide %} people-grid-wide{% endif %}">
        {% for p in group.people %}
        <div class="card person{% if group.wide %} person-wide{% endif %}">
          <div class="person-photo">
            {% if p.photo != blank %}
            <img src="{{ p.photo | relative_url }}" alt="Photo of {{ p.name }}" loading="lazy"{% if p.photo_position %} style="object-position: {{ p.photo_position }}"{% endif %}>
            {% else %}
            <img class="person-placeholder" src="{{ '/assets/img/people/placeholder.svg' | relative_url }}" alt="" loading="lazy">
            {% endif %}
          </div>
          <div class="card-body">
            <div class="person-name">{{ p.name }}</div>
            <div class="person-role">{{ p.role }}</div>
            {% if p.affiliation %}<div class="person-affiliation">{{ p.affiliation }}</div>{% endif %}
            {% if p.bio != blank %}<div class="person-bio">{{ p.bio | markdownify }}</div>{% endif %}
            {% if p.lattes != blank %}<a class="person-link" href="{{ p.lattes }}" target="_blank" rel="noopener">Lattes CV</a>{% endif %}
          </div>
        </div>
        {% endfor %}
      </div>
    </section>
    {% endfor %}

    <section id="join">
      <h2 class="category">Join the lab</h2>
      <div class="join-lab">
        <p><strong>Postgraduate students.</strong> Prof. Pedro Taucci advises in two graduate programs at Unesp: the <strong>PPG Biociências Interunidades</strong> (Bauru and Assis campuses) and the <strong>PPG Biodiversidade</strong> (São José do Rio Preto campus). Selection processes usually take place around <strong>October/November</strong> (Biociências Interunidades) and <strong>May/June</strong> (Biodiversidade) of each year. Check each program's website for the calls and requirements.</p>
        <p><strong>Prospective postdoctoral fellows.</strong> Get in touch to look for funding together, for example through a <a href="https://fapesp.br" target="_blank" rel="noopener">FAPESP postdoctoral fellowship</a>.</p>
        <p>Contact details are at the bottom of the <a href="{{ '/' | relative_url }}">Home</a> page.</p>
      </div>
    </section>
  </div>
</div>
