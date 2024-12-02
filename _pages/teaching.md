---
layout: page
permalink: /teaching/
title: Teaching
description: Materials for undergraduate and graduate courses taught.
nav: true
nav_order: 6
display_categories: [Undergraduate, Graduate]
horizontal: false
---

<div class="teaching">
{% if site.enable_course_categories and page.display_categories %}
  <!-- Display categorized courses -->
  {% for category in page.display_categories %}
  <a id="{{ category }}" href=".#{{ category }}">
    <h2 class="category">{{ category }}</h2>
  </a>
  {% assign categorized_courses = site.teaching | where: "category", category %}
  {% assign sorted_courses = categorized_courses | sort: "importance" %}
  <!-- Generate cards for each course -->
  {% if page.horizontal %}
  <div class="container">
    <div class="row row-cols-1 row-cols-md-2">
    {% for course in sorted_courses %}
      {% include teaching/courses_horizontal.liquid %}
    {% endfor %}
    </div>
  </div>
  {% else %}
  <div class="row row-cols-1 row-cols-md-3">
    {% for course in sorted_courses %}
      {% include teaching/courses.liquid %}
    {% endfor %}
    </div>
  </div>
  {% endif %}
  {% endfor %}
{% else %}
  <!-- Display courses without categories -->
  {% assign sorted_courses = site.teaching | sort: "importance" %}
  <!-- Generate cards for each course -->
  {% if page.horizontal %}
  <div class="container">
    <div class="row row-cols-1 row-cols-md-2">
    {% for course in sorted_courses %}
      {% include teaching/courses_horizontal.liquid %}
    {% endfor %}
    </div>
  </div>
  {% else %}
  <div class="row row-cols-1 row-cols-md-3">
    {% for course in sorted_courses %}
      {% include teaching/courses.liquid %}
    {% endfor %}
    </div>
  </div>
  {% endif %}
{% endif %}
</div>
