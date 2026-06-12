#!/bin/bash
export GIT_COMMITTER_NAME="Claude"
export GIT_COMMITTER_EMAIL="claude@anthropic.com"
export GIT_AUTHOR_NAME="Claude"
export GIT_AUTHOR_EMAIL="claude@anthropic.com"
git commit -m "fix: clear 2 WARNINGs + 1 NOTE from R CMD check (Batch 1)

A. Remove utils from Imports: in DESCRIPTION (utils::globalVariables
   works without it; utils is always available as a base package).
   data.table stays -- 8 :: calls remain in prep_heatmap_data.R.

B. serocalculator::as_case_data -> serocalculator:::as_case_data
   with TODO comment; R CMD check reports it unexported, demoting
   WARNING to the existing ::: NOTE.

D. Replace [serodynamics::run_mod_pop()] and [serocalculator::as_case_data()]
   auto-links with plain backtick code in both .R source and .Rd
   files -- fixes Missing link WARNINGs.

E. Add ^\.lintr\$ to .Rbuildignore -- removes Found hidden files NOTE.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
