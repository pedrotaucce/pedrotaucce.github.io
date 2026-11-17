# LaBESTA website

Website of the **Laboratório de Biologia Evolutiva e Sistemática de Tetrápodes Atuais** (Laboratory of
Evolutionary Biology and Systematics of Current Tetrapods) at São Paulo State University (Unesp), Bauru, Brazil.

Live site: <https://labesta.org>

Built with [Jekyll](https://jekyllrb.com/) and the [al-folio](https://github.com/alshedivat/al-folio) theme
(MIT license, see `LICENSE`). Every push to `master` is built and published by GitHub Actions
(`.github/workflows/deploy.yml`).

## Where things live

| What | Where |
|---|---|
| People (names, roles, photos, bios, Lattes) | `_data/people.yml`; photos in `assets/img/people/` |
| News | one short file per item in `_news/` |
| Publications | `_bibliography/papers.bib`; PDFs in `assets/pdf/`, thumbnails in `assets/img/publication_preview/` |
| Courses (Teaching) | `_teaching/`; syllabi sources in `ementas/`; card images in `assets/img/teaching/` |
| Tutorials (Portuguese) | `_tutorials/`; downloadable files in `assets/tutorials/<tutorial>/` |
| Projects | `_pages/projects.md` (under construction) |
| CV | `assets/json/resume.json` and the PDF in `assets/pdf/CV/` |
| Site settings (title, social links, navigation) | `_config.yml` |

## Adding a news item

Create `_news/YYYY-MM-DD_short-name.md`:

```markdown
---
layout: post
date: 2026-12-01 08:00:00
inline: true
related_posts: false
---

Text of the news item, with [links](https://example.com) and *italics*.
```

## Adding a paper

Add an entry to `_bibliography/papers.bib` with the fields `website` (DOI link), `abstract`, `pdf` (file name in
`assets/pdf/`) and `preview` (file name in `assets/img/publication_preview/`) to get the Abs / PDF / Website buttons.

## Images

Please resize photos before uploading (about 1600 px on the long side is plenty): large originals slow down the
build and the site.

## Running the site locally

You need Ruby (3.3 or newer) and ImageMagick.

```bash
bundle install
bundle exec jekyll serve
```

Then open <http://localhost:4000>.
