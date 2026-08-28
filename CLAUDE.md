# Claude Code Instructions

Project guidance for Claude Code --- the CLI, the IDE extension, and the
`@claude` GitHub Action.

[`AGENTS.md`](AGENTS.md) is the cross-agent contract and carries the same
conventions; read it rather than duplicating its content here.
[`.github/copilot-instructions.md`](.github/copilot-instructions.md) is the
**source of truth** for repository-specific structure, style, and workflow.
This file adds only what is specific to Claude Code.

Authoritative style guide: the
[UCD-SERG Lab Manual](https://ucd-serg.github.io/lab-manual/)
(source: <https://github.com/UCD-SERG/lab-manual>).

## Start here

Read [`AGENTS.md`](AGENTS.md) first.
Its "Conventions that trip people up" section covers the four things that fail
CI here: the inverted `Version:` rule, the two separate accept-lists,
generated documentation, and the `Morrison-Lab/gha` caller pins.

## Review scope

When reviewing a PR in this repository, judge it against
`.github/copilot-instructions.md` and the lab manual, in that order.
`copilot-instructions.md` carries a "Things Not to Flag" section; honour it.

Two review notes specific to this repository:

- `inst/scripts/` holds chapter analysis code rather than package code.
  It is held to a looser standard than `R/`, and it carries a known lint
  backlog.

- Stan and JAGS files under `inst/extdata/` are not lintable as R, and their
  numerical content --- priors, parameter names, simulation semantics --- is
  not to be changed for style.

## Editing workflows

`.github/workflows/` is mostly thin callers of `Morrison-Lab/gha` reusable
workflows.
The integrated `GITHUB_TOKEN` cannot push changes under `.github/workflows/`,
so an `@claude` run that edits one needs a `WORKFLOW_TOKEN` secret; without it
the push is rejected and the commits are posted to the thread as a patch.

The two Claude workflows are a pair.
`claude.yml` dispatches `claude-code-review.yml`, so keep that filename, and
keep both pinned to the same major tag.

## Data handling

Never commit real *Shigella* participant data.
`check-phi` scans added lines on every PR, and `check-secrets` scans the whole
git history --- but neither is a substitute for not committing it.
Synthetic values that trip `check-phi` can be exempted with a `phi-allow`
comment on the line.
