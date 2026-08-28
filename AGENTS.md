# Agent Instructions

Project guidance for AI coding agents working in `shigella` --- Codex, Cursor,
Gemini, opencode, and anything else that reads `AGENTS.md`.
Claude Code reads [`CLAUDE.md`](CLAUDE.md), which says the same things.

[`.github/copilot-instructions.md`](.github/copilot-instructions.md) is the
**source of truth** for repository-specific structure, style, and workflow.
This file is a short orientation and a list of the conventions that most often
catch people out.
Defer to `copilot-instructions.md` for anything not stated here, and edit that
file rather than this one when the two would disagree.

Authoritative style guide: the
[UCD-SERG Lab Manual](https://ucd-serg.github.io/lab-manual/)
(source: <https://github.com/UCD-SERG/lab-manual>).

## Project context

`shigella` is an R package for multivariate Bayesian hierarchical modeling of
antibody response trajectories after confirmed *Shigella* infection.
It is maintained by UCD-SERG (Aiemjoy and Morrison labs).

## Repository layout

- `R/`, `man/`, `NAMESPACE`, `DESCRIPTION` --- package source and generated docs
- `inst/extdata/` --- Stan and JAGS model files
- `inst/scripts/` --- chapter analysis scripts, not package code
- `tests/` --- test suite
- `vignettes/articles/` --- Quarto manuscripts and the dissertation source
- `inst/WORDLIST` --- accepted spellings for {spelling}
- `_typos.toml` --- where accepted terms for crate-ci/typos go, at the repo root
- `.github/workflows/` --- CI, mostly thin callers of `Morrison-Lab/gha`

## Conventions that trip people up

Each of these fails CI rather than merely being frowned upon.

### Never touch `Version:` in a PR

`DESCRIPTION`'s `Version:` must match `main`'s exactly.
The `version-check` workflow fails a PR that changes it.
`bump-dev-version` bumps the dev counter after every merge to `main`, so a PR
never needs to.

This inverted on 2026-08-27.
The older convention required every PR to increment the version, so older
branches, older instructions, and habit all point the wrong way.
If a deliberate version change is genuinely needed --- a release --- apply the
`no version increment` label.

### Two accept-lists, with different owners

Adding a word to the wrong one silently fails.

- `inst/WORDLIST` --- read by `spelling::spell_check_package()`.
  Covers `DESCRIPTION`, `man/*.Rd`, vignette sources, and root Markdown.
  Sorted by codepoint, so uppercase entries precede lowercase ones.

- `_typos.toml`, at the repo root --- read by crate-ci/typos.
  Covers everything the R spellchecker cannot see, diff-scoped to what a PR
  adds.
  Put a domain abbreviation under `[default.extend-words]`, and prefer fixing
  a real typo in the source over listing it here.
  Create the file if it is not there yet; typos discovers it automatically.

### Documentation is generated

Run `devtools::document()` after editing any roxygen2 comment.
`R-check-docs.yml` re-runs `roxygen2::roxygenise()` and fails when the result
differs from what is committed.
Never hand-edit `man/` or `NAMESPACE`.

### CI is mostly `Morrison-Lab/gha`

Most workflows under `.github/workflows/` are thin callers of reusable
workflows, pinned to `@v2`.
Two things follow.

Read the callee at its pinned tag before changing a caller's `with:` or
`secrets:` block.
A caller that passes a secret the pinned tag does not declare is rejected
before any job starts, and that failure has no logs, no annotations, and no
check run --- so it is invisible in `gh pr checks` and on the PR.
`gh run list --workflow <file>` is what surfaces it.

Prefer fixing a shared problem in `Morrison-Lab/gha` over reintroducing a
hand-maintained workflow here.

## Before you open a PR

Run these locally; each has a CI counterpart that will fail otherwise.

```r
devtools::document()
lintr::lint_package()
spelling::spell_check_package()
devtools::check()
```

Add a `NEWS.md` entry under `# shigella (development version)`, or apply the
`no changelog` label when the change genuinely does not warrant one.

Write Markdown with semantic line breaks --- one sentence per line.
`check-new-line-breaks` scans the lines a PR adds and fails on lines packing
more than one sentence together.

## Known-failing checks

`lint-changed-lines` currently reports a backlog in `inst/scripts/`, carried in
from the dissertation branches.
It is pre-existing rather than something a given PR introduced.
Do not treat clearing it as in scope for unrelated work, and do not suppress it
per line either; the disposition is tracked separately.
