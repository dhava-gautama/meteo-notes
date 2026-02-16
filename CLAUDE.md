# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is a **documentation-only** repository — a personal knowledge base of Markdown notes on meteorological modeling, observations, and forecast verification. There is no source code, build system, or test suite.

## Structure

- `models/` — Numerical model guides (WRF, MPAS, ROMS, COAWST, SWAN, WW3)
- `verification/` — Forecast verification methods: theory (`forecast-verification-notes.md`), MET/METplus tooling (`met-verification-guide.md`), Python-only approach (`python-verification-guide.md`)
- `observations/` — Surface (`surface-obs-guide.md`) and ocean (`ocean-obs-guide.md`) observation guides
- `verification/images/` — Figures referenced by the verification notes

## Writing Conventions

- Each guide is a single, comprehensive Markdown file covering end-to-end workflows
- Guides include inline code blocks (bash, Python, Fortran, namelist snippets) as practical examples — these are reference snippets, not runnable project code
- Topics have a regional focus on **Indonesian waters and BMKG** (Indonesian meteorological agency) where applicable
- Model guides reference specific software versions (e.g., WRF v4.7, MPAS v8.3, ROMS, SWAN v41.51, WW3 v6.07, MET v12.0/METplus v6.0)

## When Editing or Adding Notes

- Keep each topic as a self-contained Markdown file — avoid splitting into many small files
- Update `README.md` when adding or renaming guide files (it maintains a table of contents with one-line descriptions)
- Preserve the existing depth and style: comprehensive single-file guides with code examples, configuration references, and workflow steps
