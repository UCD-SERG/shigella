---
editor: 
  markdown: 
    wrap: 72
---

# Copilot Instructions for shigella

## Repository Overview

**shigella** is an R package for multivariate Bayesian hierarchical
modeling of antibody response trajectories following confirmed Shigella
infection. Chapter 2 extends the univariate Chapter 1 approach (handled
upstream in the `serodynamics` package) to a Kronecker-structured
covariance model ($\Sigma_B \otimes \Sigma_P$) implemented in Stan.

-   **Type**: R package (statistical modeling)
-   **Language**: R (\>= 4.1.0), Stan (cmdstanr backend)
-   **Key Dependencies**: cmdstanr (\>= 0.9.0), serodynamics, posterior,
    cli, dplyr, MASS, tibble
-   **Lifecycle**: Experimental — Chapter 2 simulation is the active
    branch (`chapter2-stan-simulation`)
-   **Repository owner**: UCD-SERG (Aiemjoy / Morrison labs)

## Lab-Wide Guidance

Follow the [UCD-SeRG Lab Manual](https://ucd-serg.github.io/lab-manual/)
for culture, reproducibility, GitHub workflows, coding practices, and AI
tools usage. If the web version is inaccessible, refer to the [source
files on GitHub](https://github.com/UCD-SERG/lab-manual).

## Repository Structure

This package follows strict standard R package layout. **Do not create
subdirectories under `R/`, and do not place `R/`, `inst/`, or
`vignettes/` folders inside any other directory** (e.g., a former
`chapter2/R/` subfolder was a structural error and was removed; do not
recreate that pattern).

### Layout

-   **`R/`**: All exported and internal R functions. One function per
    file. File name must match function name (dotted internal helpers
    like `.helper_fn()` go in `R/helper_fn.R` — drop the leading dot in
    the filename only).
-   **`inst/stan/`**: Stan model files (`model_1.stan`, `model_2.stan`).
-   **`inst/examples/`**: One example script per exported function,
    named `<function_name>-examples.R`.
-   **`tests/testthat/`**: Unit tests, named `test-<function_name>.R`.
-   **`vignettes/`**: Quarto vignettes (`.qmd`). Chapter 2
    manuscript-style documentation goes here.
-   **`scripts/`**: Shiva HPC orchestration scripts (Phase 0/1
    diagnostic runs). Listed in `.Rbuildignore` — **not** part of the
    package build.
-   **`slurm/`**: SLURM sbatch files for Shiva. Also in `.Rbuildignore`.
-   **`man/`**: Auto-generated documentation from roxygen2 — **do not
    edit directly**.

### Configuration Files

-   **`DESCRIPTION`**: Package metadata. Keep `Imports` minimal
    (currently cli, dplyr, MASS, tibble). Optional functionality belongs
    in `Suggests`.
-   **`NAMESPACE`**: Auto-generated — **do not edit**.
-   **`.lintr`** and **`.lintr.R`**: Custom lintr configuration.
-   **`.Rbuildignore`**: Must exclude `^scripts$`, `^slurm$`, and
    `^chapter2$` (last for historical safety).

-   **`_typos.toml`**: Where accepted terms for crate-ci/typos go, at the repo root, created on demand.
-   **`AGENTS.md`** and **`CLAUDE.md`**: Short orientation files for AI coding agents, both deferring to this file as the source of truth.
    Record a repository convention here and let them point at it.

## Review Priorities — What Copilot Should Flag

### 1. R Package Structure

-   All R functions live directly in `R/`. **No subdirectories**.
-   Every exported function has:
    -   `@title`, `@description`, `@param` for every argument,
        `@return`, `@export`, and either `@examples` inline or
        `@example inst/examples/<fn>-examples.R`.
-   One function per file. File name matches function name. Internal
    helpers are dotted (`.helper_fn()`) and live in `R/helper_fn.R`.
-   Stan model files belong in `inst/stan/`, never in a nested folder.

### 2. Examples That Actually Run

-   Every exported function has a matching example file in
    `inst/examples/<function_name>-examples.R`.
-   Examples must NOT depend on confidential Shigella data files. Use
    `sim_correlated_case_data()` to generate synthetic data instead.
-   For Stan-fitting examples, use minimal settings (n = 5 subjects, 100
    warmup + 100 sampling, 1 chain) so the example completes within \~5
    minutes on CI. Wrap heavy fits in
    `if (interactive() && requireNamespace("cmdstanr", quietly = TRUE))`.

### 3. Unit Tests

-   Every function in `R/` should have a test in
    `tests/testthat/test-<function_name>.R`.
-   For Stan-fitting tests, use minimal MCMC settings and consider
    `testthat::skip_on_ci()` if the test exceeds \~5 minutes.
-   Tests for `sim_correlated_case_data()` must verify the generated
    data has the intended structure: correct dimensions, correct
    ground-truth correlation sign, the `truth` attribute is present.

### 4. Real-Data Handling

-   Real Shigella data files (anything matching `dL_clean_*.rda`,
    `data/*.rda` containing patient data, or `*.xlsx` with patient
    records) MUST be in `.gitignore` and MUST NOT be committed.
-   Flag any PR that adds files with `clean_`, `dL_`, or
    patient-identifying content.

### 5. Stan Model Files

-   Stan files (`.stan`) belong in `inst/stan/`.
-   Compiled Stan binaries (no extension, or `.exe`/`.hpp`) belong in
    `.gitignore`.
-   When a Stan model file is modified, verify:
    -   Likelihood computed in log space (no `exp` followed by `log`).
    -   Priors are weakly informative (sigma typically ≤ 5 on log-scale
        parameters).
    -   LKJ priors use `lkj_corr_cholesky`, not `lkj_corr` directly.
-   Both the full Kronecker model (`model_2.stan`) and any future
    block-covariance simplification (`model_2_block.stan`) should be
    preserved. The convergence proposal identifies all four
    cross-biomarker / cross-parameter covariance structures as
    scientifically relevant; block simplification is an operational
    trade-off, not a project-scope redefinition.

### 6. Shiva HPC Scripts

-   `scripts/` and `slurm/` are Shiva-only directories. They are not
    built into the package and not expected to run in GitHub Actions. Do
    NOT flag them for missing tests.
-   However: the R functions they call (e.g., `run_phase0_diagnostic()`,
    `run_phase1_diagnostic()`) live in `R/` and DO need to pass GitHub
    Actions tests.
-   Scripts must NOT contain function definitions or `source()` calls.
    They contain only `library(shigella)`, parameter setup, and function
    invocations. Flag any function definition or `source()` inside
    `scripts/`.
-   Phase 0 scripts are designed to run under `salloc` (interactive
    SLURM allocation), NOT on the login node directly.

### 7. Vignettes and Documentation

-   Vignettes are Quarto `.qmd` files rendered with the `quarto`
    package.
-   Chapter 2 manuscript-style documentation belongs in a single
    `vignettes/chapter2.qmd` with the standard 4-section structure
    (Introduction / Methods / Results / Discussion), all math in one
    file. Do not fragment into separate methods/demo files.

## Style Preferences

-   **Follow the tidyverse style guide**: https://style.tidyverse.org
-   **Native pipe**: `|>` not `%>%`
-   **Naming**: snake_case for functions, arguments, and most objects.
    Uppercase acronyms allowed (e.g., `Omega_B`, `tau_P`). Constants use
    UPPER_SNAKE_CASE.
-   **Maximum line length**: 120 characters (per `.lintr`).
-   **No `library()` in package code**: Use `::` or declare in
    `DESCRIPTION` Imports.
-   **`source()` is banned in the package**. Use `library(shigella)` for
    end-user contexts (vignettes, examples, scripts) or
    `devtools::load_all()` for developer workflows.
-   **No `<<-`** (global assignment).
-   **User-facing messages** use `cli::cli_inform()`, `cli::cli_warn()`,
    `cli::cli_abort()` — not `message()`, `warning()`, `stop()`. See
    `.lintr.R` for the enforced list.
-   **No developer names in user-facing messages.** Strip any reference
    to specific people ("Kwan Ho", "Ezra", etc.) from `cli::cli_*()`,
    `cat()`, `message()` output.
-   **Decompose long functions and long loops**: any R function \> \~40
    lines or any loop body \> \~5 lines should be split into named
    helpers (one function per file).

## Build and Development Workflow

### Setup

``` r
# Install development dependencies
install.packages(c("devtools", "remotes"))
remotes::install_github("UCD-SERG/serodynamics", upgrade = "never",
                        dependencies = NA)
install.packages("cmdstanr",
                 repos = c("https://stan-dev.r-universe.dev",
                           "https://cloud.r-project.org"))
cmdstanr::install_cmdstan()

# Install the package itself
devtools::install(".", upgrade = "never")
```

For minimal install (skipping `Suggests`), use base R:

``` bash
cd ~/shigella
R CMD INSTALL .
```

### Documentation Generation

**Always regenerate documentation after modifying roxygen2 comments.**

``` r
devtools::document()
```

`man/` and `NAMESPACE` are auto-generated. Do not edit directly.

### Package Checking

``` r
devtools::check()
```

Allow \~5 minutes. Acceptable post-cleanup state: 0 errors, 0 warnings,
≤1 NOTE (the `serodynamics:::calc_fit_mod` triple-colon use is a known
tolerated exception, marked with `# nolint: namespace_linter`).

### Testing

``` r
devtools::test()
```

Tests live in `tests/testthat/` and use testthat 3.0+ with edition 3
config (see `DESCRIPTION`).

### Linting

``` r
lintr::lint_package()
```

**Lint must be clean on every commit.** Treat lint failures with the
same urgency as test failures. See `.lintr.R` for the full custom
config. Custom rules include:

-   `cli::cli_*()` enforced over `message()`/`warning()`/`stop()`.
-   Native pipe `|>` enforced over `%>%`.
-   snake_case with uppercase acronyms allowed.
-   `source()` flagged.

If a lint must be suppressed, use `# nolint: <linter_name>` on the
specific line with a justifying comment, and mention it in your PR
summary.

### Spelling

``` r
spelling::spell_check_package()
```

Custom words live in `inst/WORDLIST`.
A second, separate accept-list lives in `_typos.toml` — see "Two Accept-Lists" below.

## Continuous Integration

Most workflows under `.github/workflows/` are thin callers of reusable
workflows in
[`Morrison-Lab/gha`](https://github.com/Morrison-Lab/gha), pinned to `@v2`.
Read the callee at its pinned tag before changing a caller's `with:` or
`secrets:` block: a caller passing a secret the pinned tag does not declare is
rejected before any job starts, and that failure produces no logs, no
annotations, and no check run.

Hand-maintained here, because `gha` models no equivalent:

1.  **R-CMD-check.yaml** — Runs `R CMD check` on multiple platforms.
2.  **R-check-docs.yml** — Verifies `roxygen2::roxygenise()` output matches committed `man/` and `NAMESPACE`.
3.  **check-readme.yaml** — Renders `README.Rmd`.
4.  **pkgdown.yaml** — Builds the pkgdown site.
5.  **test-coverage.yaml** — Coverage, uploaded to Codecov.
6.  **pr-commands.yaml** — `/document` and `/style` PR commands.
7.  **copilot-setup-steps.yml**, **phase0-debug.yaml** — Copilot agent environment setup, and a manual Phase 0 diagnostic run.

Called from `gha`:

-   **spellcheck.yml** — {spelling} over package prose, accepting `inst/WORDLIST`.
-   **check-typos.yml** — crate-ci/typos over the lines a PR adds, accepting `_typos.toml`.
-   **lint-changed-lines.yml** — lintr over the lines a PR adds or modifies, rather than whole changed files.
-   **version-check.yml** and **bump-dev-version.yml** — the version pair described under "Versioning" below.
-   **news.yaml** — requires a `NEWS.md` entry, bypassed with the `no changelog` label.
-   **check-junk-files.yml**, **check-secrets.yml**, **check-phi.yml**, **check-links.yml**, **check-new-line-breaks.yml** — hygiene and safety checks.
-   **lint-workflows.yml**, **lint-yaml.yml**, **lint-markdown.yml** — currently warn-only, while a pre-existing backlog is worked down.
-   **claude.yml** and **claude-code-review.yml** — the `@claude` agent and its reviewer, which are a pair.
    `claude.yml` dispatches the reviewer by filename, so keep both names and keep their pins in step.

**`lint-changed-lines` currently fails** on a pre-existing backlog in `inst/scripts/`, carried in from the dissertation branches.
It is not something a given PR introduced.
Do not treat clearing it as in scope for unrelated work.

## Versioning

**Do not change `DESCRIPTION`'s `Version:` in a pull request.**
It must match `main`'s exactly, and `version-check` fails a PR that changes it.
`bump-dev-version` bumps the dev counter after every merge to `main`, so a PR never needs to.

This inverted on 2026-08-27.
The older convention required every PR to increment the version, so older branches and older habits point the wrong way.
Apply the `no version increment` label when a deliberate version change is genuinely needed, such as a release.

## Two Accept-Lists

Adding a word to the wrong one silently fails.

-   **`inst/WORDLIST`** — read by `spelling::spell_check_package()`, covering `DESCRIPTION`, `man/*.Rd`, vignette sources, and root Markdown.
    Sorted by codepoint, so uppercase entries precede lowercase ones.

-   **`_typos.toml`**, at the repo root — read by crate-ci/typos, covering everything the R spellchecker cannot see.
    Put a domain abbreviation in `[default.extend-words]`, and prefer fixing a real typo in the source over listing it here.

## Things Not to Flag

-   `scripts/` and `slurm/` directories at the repo root — intentional
    Shiva-only directories listed in `.Rbuildignore`.
-   The `serodynamics:::calc_fit_mod` triple-colon call in
    `R/run_mod_stan.R` — known tolerated until upstream exports it.
-   Stan and JAGS files — not lintable as R code.
-   Tests appropriately skipped on CI with explanations.

## Contact

-   **Repository owner**: UCD-SERG (Aiemjoy / Morrison labs)
-   **Primary maintainer (Chapter 2)**: Kwan Ho Lee
    (ksjlee\@ucdavis.edu)
-   **Co-advisors**: Kristen Aiemjoy, Douglas Ezra Morrison

## Trust These Instructions

When making changes:

1.  **ALWAYS** keep all function definitions inside `R/`, one per file.
2.  **ALWAYS** preserve numerical/statistical logic (priors, parameter
    names, simulation semantics) unless there is a clear bug.
3.  **ALWAYS** run `lintr::lint_package()` and fix every issue in the
    same commit. Never push a commit that introduces or leaves a lint.
4.  **ALWAYS** run `devtools::document()` after modifying roxygen2.
5.  **ALWAYS** run `devtools::check()` and `devtools::test()` locally
    before requesting review.
6.  **NEVER** add a `source()` call or define a function inside
    `scripts/` or `vignettes/`.
7.  **NEVER** commit real Shigella data or files matching
    `dL_clean_*.rda` or patient data spreadsheets.
8.  **NEVER** modify `inst/stan/*.stan` files for stylistic reasons —
    only for clearly documented bug fixes.

9.  **NEVER** change `DESCRIPTION`'s `Version:` in a pull request, because `bump-dev-version` owns it — see "Versioning" above.

Only search for additional information if these instructions are
incomplete or incorrect for your specific task.
