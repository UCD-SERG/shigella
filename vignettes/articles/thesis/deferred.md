# Chapters 2 and 3 — deferred items

Raised in review against `KR_writing_standard.md`, judged real, and deliberately not
acted on before filing. None affects an argument, a number, or a conclusion. A reviewer
who flags any of these is right — that is why they are here rather than in Section 11 of
the standard, which lists items that must not be flagged.

Both chapters go into one dissertation, so they share one list. Each entry names the
chapter it belongs to.

Last reviewed: 2026-09-05. Chapter 3 at 628f67d, chapter 2 at a925bdf, dissertation at
599b546.
Section numbers refer to version 3 of the standard.

---

- [ ] **ch3 L1147** — 0.99 re-quoted from Results L1111, with "three of six" restated as
      "half the sites". Name the result, do not re-quote it. **§6.2**

- [ ] **ch3 L1192** — 29% re-quoted from Results L1070, together with the claim that the
      simulation does not reproduce it. Near-verbatim; cut one. **§6.2 §6.1**

- [ ] **ch3 L1038 / L1178 / L104** — the same 76% is "of the total" in Results and "of the
      correction" in the abstract and Discussion. Pick one noun; the abstract and
      Discussion already agree, so Results is the outlier.

- [ ] **ch3 Abstract** — 357 words. Under the dissertation ceiling; a fifth over the
      300-word journal limit. Trim at submission, not at filing. **§9.6**

- [ ] **ch3 Anchor ids** — `tbl-banks`, `tbl-two-banks`, `tbl-cells`, `sec-ch3-banks` still
      carry the vocabulary §2.2 retired. Invisible to readers; costs every future
      reviewer five grep hits to re-adjudicate. **§2.2**

- [ ] **ch3 Anchor id** — `sec-ch3-honest` → `sec-ch3-prespec`, matching its heading.
      Produces two false positives on check 6 as it stands.

- [ ] **ch3 L1001 / L1066** — two Results sentences of near-identical shape describing
      different quantities (interval width; point estimate). Name the quantity at the
      head of each — a careful reader misread the second as the first on two separate
      passes.

- [ ] **ch3 in the dissertation** — the chapter refers to its source as the kinetics
      analysis sixteen times and cites `@lee2026correlated` five times, and never says
      chapter 2. §2.5 asks for "the kinetics analysis of chapter 2" at first mention, and
      it is Kristen's own comment. Correct in the standalone, where the companion paper is
      a paper; in the dissertation it is chapter 2, and the citation stands where a
      cross-reference would serve. Fixing it in the standalone would be wrong, and the text
      sits inside a transplanted region, so a thesis-only fix is reverted by the next
      re-transplant. **§2.5**

- [ ] **ch2 L812 / L815** — the coverage figure's internal labels read `null (c = 0)` and
      `correlation present` while every caption and every prose use now says null and
      alternative. The figure has to be regenerated on Mercury to change them, and §2.7
      permits a legend to stay descriptive, so this is a decision about the figure
      rather than about the prose. **§2.7**

- [ ] **ch2 L1012** — "computed without a model" and "use no model at all" say the same
      thing eleven words apart. The bold marks the second as emphasis rather than
      repetition, and the sentence carries the chapter's strongest claim about what is
      evidence and what is commentary, so tightening it means choosing which half to
      cut. **§6.1**

Not mine to decide:

- [ ] **ch3 Abstract, §9.1 and §9.4** — the first sentence names only "enteric fever", and
      participant selection is not stated. Both rules trace to Kristen's comments 7 and
      10 and the wording is hers. **OWNER: Kristen — ask, do not decide.**
