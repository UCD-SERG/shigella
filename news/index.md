# Changelog

## shigella (development version)

- Initial CRAN submission.

### New API features

### Bug fixes

### Internal changes

- added `.lintr` configuration file

- Moved helper functions into `R/` subdirectory (#1)

- simplified spell-check GitHub Action (#3)

- Migrated the Claude Code and NEWS-changelog workflows to
  `d-morrison/gha` reusable workflows, and removed the obsolete
  upstream-fix watcher (#29)

- Finished the `Morrison-Lab/gha` migration: replaced the
  hand-maintained spellcheck, changed-file lint, and version-check
  workflows with their gha counterparts, added the `bump-dev-version`
  companion, and added checks for workflow, YAML, and Markdown lint,
  junk files, typos, committed secrets, PHI, broken links, and semantic
  line breaks (#32)

- Moved the Claude, Claude-review, and NEWS callers from `@v1` to `@v2`,
  fixing the startup failure every Claude Code run hit because `@v1`
  does not declare the `ANTHROPIC_API_KEY` secret the callers pass (#32)

- Recorded the roxygen2 version in `Config/roxygen2/version`, which
  roxygen2 8.1.0 uses in place of `RoxygenNote` (#32)

- Pinned `use-ai-config: true` explicitly on the Claude and
  Claude-review `Morrison-Lab/gha@v2` callers (#36). The reusable
  workflows already default that input to true, so the bot and the
  reviewer were already installing the shared `ai-config` plugin on
  every run. Setting it in the caller records that intent and keeps it
  from changing if the `@v2` default ever moves.

- Added `AGENTS.md` and `CLAUDE.md`, short orientation files for AI
  coding agents that defer to `.github/copilot-instructions.md` as the
  source of truth, plus a `.claude/settings.json` permissions allowlist
  (#35)

- Corrected `.github/copilot-instructions.md`, whose CI section was
  stale, and recorded the inverted version convention and the two
  separate spelling accept-lists there (#35)

## shigella 0.0.0

Started development.
