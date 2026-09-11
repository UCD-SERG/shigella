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

- [ ] **ch3 — the cluster-robust result is no longer promised in advance** —
      Kristen's second round removed the sentence in the Introduction
      ("We report it whether or not it does...") and rewrote the Methods sentence
      that carried the same commitment ("we report it either way").
      Both are her rule 8 working as designed: commentary about the reporting process,
      stated once in Methods or not at all.
      The chapter still does the thing — Results says
      "The cluster-robust check did not support the argument, and we report it as such" —
      but it no longer says in advance that it will.
      **Rule: a prespecification is only a prespecification if it is recorded before the result.
      This one is now a description of what happened.**

- [ ] **ch3 — Results claims a prespecification that Methods no longer records** —
      Results says "Two conclusions were pre-specified in @sec-ch3-app",
      but the paragraph in that section which named them and called them pre-specified
      was removed when her new opening named them instead.
      Her sentence names the two conclusions; it does not say they were fixed in advance.
      **The claim in Results is still true of what was done and is no longer supported by what §2.9 says.**

- [ ] **ch3 — two anchors lost their only in-text reference** —
      `@eq-contrasts` and `@sec-ch3-banks`.
      Both fell out of her replacement sentences:
      the contrasts equation was cited in the sentence her §2.9 opening replaced
      and in the "established internally" paragraph she deleted,
      and the parameter-draw section was cited only from the old Objectives.
      Neither is an unresolved reference — both are still defined, and the render is clean.
      `@tbl-design` was orphaned the same way and was given a reference,
      because a table that appears in the List of Tables with no path to it is a reader's problem;
      an equation and a section are not.

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

From Kristen's second round — content she asked for and we do not have:

Grouped by what each one needs, because the groups close in different ways.

NEEDS THE PARTICIPANT-FLOW FIGURE REGENERATED (Mercury)

- [ ] **ch3 @fig-ch3-flow — three changes to the figure, one of which changes a claim** —
      she asked for "rows" to be replaced by "observations";
      for the flow to start at "population-based serosurveys" rather than at the
      top-level exclusion, which she says only existed because all the data arrived in
      one spreadsheet; and for "community participants" to become
      "single index participant per household".
      **The third is not a relabelling. It changes what the figure asserts about who was
      sampled — from a community sample to one index participant per household —
      and the surrounding text and @tbl-cells describe the first.
      Anyone regenerating this figure should settle which is true before drawing it,
      because the label and the text cannot both be right.**
      TO CLOSE: the plotting script for `ch3_fig_flow.png` has to be re-run on Mercury,
      where the participant-level data lives. The first two changes are label edits in
      that script. The third needs the sampling question below answered first.

NEEDS FACTS FROM KRISTEN OR FROM THE SEES PAPER

- [ ] **ch3 §2.1 — the sampling design is not described** — she asked for
      "two phase geographically random sample, census first and then households randomly
      selected". The chapter says only that households were sampled in high- and
      low-endemicity sites.
      TO CLOSE: one or two sentences, from the SEES design paper [@aiemjoy2022] or from her.
      Nothing here can be inferred from the repository.

- [ ] **ch3 §1.1 — background on HlyE** — she asked for how the antigen was identified and
      where it has been evaluated, and separately for "+references to charles work"
      to follow the sentence describing it.
      Her tracked edit to that sentence was applied; the note attached to it was not,
      because it names a person rather than a citation.
      TO CLOSE: the references themselves. Which of Charles's papers, and whether the
      background is a clause or a sentence, are hers to say.

NEEDS NUMBERS WE CAN COMPUTE BUT DID NOT REPORT

- [ ] **ch3 §2.1 — counts behind the sex percentages** — she marked `48%` with "(n/N)" and
      `35%` and `53%` with "(n/n)(n/n)", asking for the counts as well as the proportions.
      TO CLOSE: a calculation, not a fact — the numerator and denominator overall and for
      the two extreme sites. The participant-level data is on Mercury; the percentages in
      the chapter were computed there and only the percentages were carried back.

- [ ] **ch3 §2.1 — median age and IQR** — she asked for them overall and by site.
      @tbl-cells already carries an age distribution by site, so the by-site half may be a
      restatement rather than a gap; the overall figures are not in the chapter at all.
      TO CLOSE: the same calculation on Mercury, and a decision about whether the overall
      pair goes in the text or a row goes in the table.

NEEDS A DECISION RATHER THAN A FACT

- [ ] **ch3 Abstract — "the true seroincidence"** — two comments on the same words,
      "simulated true?" and "simulated?". In the coverage study the truth is the simulated
      λ, so "true" is correct as written; she is asking whether the abstract should say so.
      Nothing is missing and nothing needs computing.
      TO CLOSE: a wording decision — leave "true", or write "the simulated seroincidence"
      and accept that the abstract then reads as though the chapter only simulated.
      **Neither of the two comments proposes wording, so this is not a tracked edit to apply.**

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

- [ ] **ch3 §2.5 — her second rewrite was taken and her first was not** — two of
      her comments supply a replacement for the same idea:
      `401576599` rewrites "The same function evaluates all three...",
      and `1880105609` rewrites "Placing all three forms in one function...",
      twenty-five lines later.
      Both open on the single-implementation point and both contain the sentence
      "any differences between the likelihoods can be attributed directly to...".
      `1880105609` is introduced as a "potential re-write",
      which means she wrote it knowing the first existed, so it was taken and the first was not.
      A reader meeting the same conclusion twice in twenty-five lines
      is worse served than one meeting it once.
      **This was the author's decision, not hers.
      Reverting it means `401576599` goes in as well, not instead —
      and with it "The same function evaluates all three." and "Two contrasts follow,",
      both of which her own writing guide quotes as examples of what to avoid.**

For the journal revision:

- [ ] **chapter 1 — "Section B5 in S1 Text" was a broken reference in the submitted
      version** — the supplement has A1 and B1 through B4; there is no B5.
      The text it pointed at, in a sentence about Gamma and Wishart priors,
      is what Appendix A's Section 8.6 now holds —
      a section the supplement does not have.
      The dissertation resolves it to the appendix.
      **The PLOS version still carries the dangling pointer
      and should be corrected at revision.**

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
      Her second round deletes the same sentence again, on a file that still carried it,
      so she has now asked for it twice in two independent passes.
      That is evidence about what she wants, not a change of position —
      but it means the question to put to her is which of the two she meant, not whether she meant it.
      **OWNER: Kristen — ask, do not decide.**
