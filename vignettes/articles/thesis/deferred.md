# Chapters 2 and 3 — deferred items

Raised in review, judged real, and deliberately not acted on before filing.
None affects an argument, a number, or a conclusion.
A reviewer who flags any of these is right —
each entry names the rule it falls under and says why it was left.

Both chapters go into one dissertation, so they share one list.
Each entry names the chapter it belongs to.

Last reviewed: 2026-09-09.
Chapter 3 at `cde69d6`, chapter 2 at `a925bdf`, dissertation at `6d31d7c`.

---

- [ ] **ch3 L1147** — 0.99 re-quoted from Results L1111,
      with "three of six" restated as "half the sites".
      Name the result, do not re-quote it.
      **Rule: a discussion may name a result rather than re-quoting its value.**

- [ ] **ch3 L1192** — 29% re-quoted from Results L1070,
      together with the claim that the simulation does not reproduce it.
      Near-verbatim;
      cut one.
      **Rule: say it once, and refer to the section holding the detail.**

- [ ] **ch3 L1038 / L1178 / L104** — the same 76% is "of the total" in Results
      and "of the correction" in the abstract and Discussion.
      Pick one noun;
      the abstract and Discussion already agree, so Results is the outlier.

- [ ] **ch3 Abstract** — 357 words.
      Under the dissertation ceiling;
      a fifth over the 300-word journal limit.
      Trim at submission, not at filing.
      **Rule: an abstract should meet its target journal's word limit.**

- [ ] **ch3 Anchor ids** — `tbl-banks`, `tbl-two-banks`, `tbl-cells`, `sec-ch3-banks`
      still carry vocabulary the prose retired.
      Invisible to readers;
      costs every future reviewer five grep hits to re-adjudicate.
      **Rule: a term appearing in no paper you cite should be dropped, or defined where it first appears.**

- [ ] **ch3 Anchor id** — `sec-ch3-honest` → `sec-ch3-prespec`, matching its heading.
      Produces two false positives on check 6 as it stands.

- [ ] **ch3 L1001 / L1066** — two Results sentences of near-identical shape
      describing different quantities (interval width; point estimate).
      Name the quantity at the head of each —
      a careful reader misread the second as the first on two separate passes.

- [ ] **ch3 in the dissertation** — the chapter refers to its source as the kinetics analysis sixteen times
      and cites `@lee2026correlated` five times, and never says chapter 2.
      Naming the analysis at first mention is what is asked for,
      and it is Kristen's own comment.
      Correct in the standalone, where the companion paper is a paper;
      in the dissertation it is chapter 2, and the citation stands where a cross-reference would serve.
      Fixing it in the standalone would be wrong,
      and the text sits inside a transplanted region,
      so a thesis-only fix is reverted by the next re-transplant.
      **Rule: name the analysis you mean rather than referring to it obliquely.**

- [ ] **ch2 L812 / L815** — the coverage figure's internal labels read `null (c = 0)` and `correlation present`
      while every caption and every prose use now says null and alternative.
      The figure has to be regenerated on Mercury to change them,
      and a legend may stay descriptive because it has no preceding sentence to carry the definition,
      so this is a decision about the figure rather than about the prose.
      **Rule: one term per object, used consistently in prose.**

- [ ] **ch2 L1012** — "computed without a model" and "use no model at all"
      say the same thing eleven words apart.
      The bold marks the second as emphasis rather than repetition,
      and the sentence carries the chapter's strongest claim about what is evidence and what is commentary,
      so tightening it means choosing which half to cut.
      **Rule: say it once.**

- [ ] **ch1 L420** — the convergence-diagnostics paragraph narrows the pooled model's failure to *S. sonnei* IgG,
      while Methods L351, Results L471 and Discussion L537 all carry it unqualified as *S. sonnei*.
      L471 extends it explicitly ("*S. sonnei* IgA showed a similar pattern")
      and the abstract and the closing chapter both use the unqualified form.
      One paragraph reporting less than the other three, not a contradiction;
      a two-word addition to L420 would settle it.
      Chapter 1 is submitted, so this is a note for the next revision rather than a change to make now.

- [ ] **ch1 — *S. flexneri* 6, a question to be ready for** — the abstract and L2822 say
      that pooling antigenically distinct serotypes produced implausible decay.
      L307 calls Sf6's O-antigen structurally distinct, Sf6 was pooled,
      and no convergence failure is reported for it;
      L491 gives its pooled MAE as 0.19--0.22, which looks good.
      The answer is L477:
      under the pooled model Sf6 shows the flattest trajectories, near baseline for most participants,
      and a flat trajectory is easy to predict,
      so a low MAE there is not evidence the model is right.
      L543 adds that no alternative was fitted (n = 5)
      and L541 groups Sf6 with *S. sonnei* as requiring dedicated assays.
      Not a defect in the text — an answer that lives in the body and not in the abstract.

- [ ] **Appendix A — where a line-by-line comparison should start** — the appendix
      reproduces the supplementary methods submitted with chapter 1,
      and carries three things the supplement does not:
      a section on priors,
      a closing paragraph on why chapter 2 cannot use the Wishart prior,
      and — inside the shared hierarchical-parameterization section — the clause
      "which is the assumption chapter 2 relaxes".
      The first two are additions after the shared text;
      the third is an insertion within it,
      so the appendix cannot be treated as a copy with separable additions.
      Any line-by-line comparison against the submitted supplement should start at that clause.
      Chapter 1 is submitted, so nothing here is to be changed now.

Considered and declined:

- [ ] **ch3 L306 — Kristen's `the same` → `kinetic` (tracked changes, runs 38-39)** —
      declined, and the request it carries is met elsewhere.
      Her paired edit "in the companion analysis" → "described previously" (runs 40-41)
      has no target: the sentence now reads "as in the longitudinal cohorts".
      Applying 38-39 alone breaks the clause they lead into —
      "so the parameter draws and the specimens they are applied to are on one scale"
      follows from the assay being *the same* as the one behind the draws,
      and "kinetic ELISA as described previously" drops both the sameness and the referent.
      It would also orphan the (ELISA) definition, which L318 uses and nothing else supplies.
      Her "as described previously" is the same request as her Methods note,
      and that is now carried by L304, "described in detail previously".
      **Rule: an edit is applied only where its target survives and its logic still closes.**

- [ ] **ch3 L156 — Kristen's goal-statement placeholder (tracked changes, run 25)** —
      declined; the slot is empty but its job is done.
      She deleted the signpost at the head of that paragraph
      ("The last clause is the one this chapter examines", run 20)
      and left her placeholder for a goal statement at the tail.
      Those are one move: relocate the signpost.
      `f3938c7` independently rewrote the head into today's L154,
      "and it is that consequence this chapter examines" — a sentence she never saw.
      Filling the tail as well would give one paragraph a signpost at both ends,
      and make three such statements between L154 and Objectives at L290.
      **Rule: a placeholder is closed by the request being met, not by the slot being filled.**

From Kristen's tracked changes — notes not closed:

- [ ] **ch3 Introduction — note A, the base method** — she asked for three things to be
      introduced in the Introduction: the implementation, the within-host decay model,
      and how the model pairs with cross-sectional data.
      The seam sentence at L143 names the kinetic model and the combination with the
      survey data, and the implementation is named one line into the next section, at L147.
      The third piece is done, and the first two belong to that section rather than to
      the paragraph.
      Closing them fully would mean writing the paragraph she deleted.

- [ ] **ch3 L158–159 — note C, the contrast Kristen was drawing** — she asked for recency
      to be presented as information in addition to correlation.
      The chapter establishes the separation and quantifies it, at L105 and L1025–1027:
      the shared infection time accounts for 76% of the correction
      and requires no result about correlation.
      But the Introduction never draws that contrast.
      L159 contrasts the pair with a single marker, not with correlation,
      and its colon clause states the general point itself rather than only illustrating it,
      so a sentence adding the correlation contrast above it pre-empts its own conclusion.
      Her run 29, a lowercasing, is the tail of this note and depends on it.
      **The gap is that the Introduction contrasts the pair with one marker
      where she wanted it contrasted with correlation.**

- [ ] **ch3 Methods — note E, the longitudinal case data** — she asked for the case data
      and the model to be described where the cross-sectional data are.
      Between L304 and the posterior parameter draws at L427
      the chapter never says what the longitudinal cohorts were: who, how many, what design.
      In the dissertation this is covered by chapters 1 and 2;
      in the standalone manuscript it is not.

- [ ] **ch3 — note F, a Table 1** — she asked for a study-population table
      with columns by site and rows for age, sex and collection dates,
      modelled on the SEES paper.
      The rows exist but are split three ways:
      `tbl-cells` carries sample size, age, censoring and the observed correlation by site,
      `tbl-collect` carries the collection dates,
      and sex appears only as a clause at L324.
      Building the single table needs participant-level data that is not in this repository.

Process:

- [ ] **ch3 — Kristen's tracked changes were not applied in the first pass** — her review
      of 2026-08-26 left ten comments on the abstract and seven tracked edits to it.
      The comments were all addressed; none of the edits were.
      Across the chapter, fourteen of her eighteen prose deletions were still in the file
      until `bc8fd64`.
      Her comments were tracked as a numbered list and her edits were not,
      and the pass that applied the comments went in inside a 750-line rewrite
      whose message names neither her nor the abstract,
      so the history shows when the text changed but not whether she was answered.
      This is why the four entries above exist.
      **Rule: a review is closed against the reviewer's own file, not against the commit log.**

Settled — do not undo:

- [ ] **reference.docx — the dissertation's copy has diverged deliberately** — the custom
      `Table` style carries `tblCellMar` left and right at 0
      in `vignettes/articles/thesis/reference.docx`,
      and at 80 in the two chapter copies.
      The zero pulls the caption rule inside the one-inch margin,
      which UC Davis requires of figures and tables;
      the journals have no such rule and their copies are correct as they are.
      **Anyone synchronising the three files should not copy 80 back over the zero.**

Not mine to decide:

- [ ] **ch3 Abstract** — the first sentence names only "enteric fever",
      and participant selection is not stated.
      Both rules trace to Kristen's comments 7 and 10 and the wording is hers.
      **Rules: an abstract should name both diseases and say how participants were selected.
      OWNER: Kristen — ask, do not decide.**

- [ ] **ch3 Abstract — she wants the level reported** — she inserted
      "and the seroincidence estimate was higher across all sites"
      and deleted "The point estimates are not in question, and no published incidence
      estimate is shown to be wrong;
      what is too confident is the uncertainty around them."
      Those are one move: report that the level went up,
      and drop the sentence that contradicts it.
      The chapter argues the opposite at L110, L935 and L1129,
      while L1053 does report that the joint estimate exceeded the product estimate
      in six of six sites — so her substantive point is in the body and not in the abstract.
      This is a disagreement with a co-author, not an oversight.
      **OWNER: Kristen — ask, do not decide.**
