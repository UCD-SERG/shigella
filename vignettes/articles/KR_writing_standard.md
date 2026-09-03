# Manuscript Writing Standard — SERG dissertation chapters

**Sources**
1. Kristen Aiemjoy's 34 review comments on `chapter3_v1_Kristen_reviewed.docx` (2026-08-26)
2. Luby S, Southern D. *A Guide to Quantitative Writing in the Health Sciences* (2022) —
   the reference Kristen named in comments #17 and #32

**Scope** — chapter1.qmd, chapter2.qmd, chapter3.qmd, and the thesis-only sections of
thesis_ch1_ch2_ch3.qmd.

⚠ This document is written to be handed to a reviewer, human or automated, as the
standard against which a draft is checked. Every rule below traces to a source. Where
the two sources conflict, §0.2 says which wins.

---

# 0. How to use this

## 0.1 Rules are numbered by origin

```
[K##]  Kristen's comment number
[G#.#] Guideline section number
[D]    Derived — follows from both, or from a decision we made
```

★ When reporting a violation, cite the tag. A finding with no tag is an opinion, and
opinions go in a separate list.

## 0.2 ⚠ Precedence when the two sources disagree

```
★★★ Kristen's explicit instruction wins over a general guideline principle.
```

**Worked example.** [G3.2] says avoid lengthy background in the introduction. [K25] says
"previous method goes in background." Following [K25] put four equations in the
introduction of chapter 3. ★ That is correct as it stands. The guideline yields.

⚠ A reviewer must not flag the equations-in-introduction as a violation.

## 0.3 What a review should and should not do

```
★ Report        rule violated · location · the text · a concrete replacement
⚠ Do not        rewrite silently · improve prose that violates nothing ·
                flag style preferences as violations
★★ Separate     "violates a rule" from "I would have written it differently"
```

---
---

# 1. Voice and person

## 1.1 [K16][K32][G1.4-Strong] Active voice

⚠ Kristen's overall verdict on v1: *"reads like a stream of consciousness. Rewrite in
clear active scientific voice."*

```
★ Target      passive constructions under 25% of sentences
★ Measure     count "was/were/is/are + past participle"
```

**Before / after from v2:**
```
before  "The likelihood was evaluated across a grid of λ."
after   "We evaluated the likelihood across a grid of λ."
```

★ v1 → v2 moved from 37.2% passive to 23.4%. ⚠ That is the standard, not an aspiration.

## 1.2 [G5.4] ★★★ "We" must be attributable

⚠ The guideline forbids attributing work to "we" when it was done by people not on the
author line.

**Violation found in chapter 3 v2:**
```
"★ We measured HlyE IgG and IgA by the same ELISA…"
```
⚠ The next sentence says the laboratories did it.

```
★ Fix     "HlyE IgG and IgA were measured by the same ELISA at the study laboratories"
★★ Rule   Field and laboratory work performed by study staff takes passive voice or
          names the actor. Analysis performed by the authors takes "we".
```

## 1.3 [D] ★★ No first person singular

⚠ These are multi-author manuscripts.

```
Violation found:  "a second version of my own code"
Fix:              "a second implementation written for this analysis"
```

★ Scan for `\bmy\b`, `\bI\b`, `\bmine\b` outside quotations.

---
---

# 2. Terminology

## 2.1 [K2] ★★ Name the quantity you estimate

```
❌ transmission · transmission intensity · force of infection
✅ seroincidence · infection rate
```

⚠ Kristen: *"we are not measuring transmission, we are measuring the infection rate."*

★ One survivor is acceptable — "transmission year" in a seasonality context — because it
names a period, not the estimand.

## 2.2 [K30] ★★★ Do not invent terms

⚠ Kristen: *"why 'bank'? Use a standard term, or define it clearly on first use."*

```
❌ bank · parameter bank · arm S · cells
✅ posterior parameter draws · paired draws · independent draws · shuffled draws
✅ study sites
```

★ **Test**: if a term appears in no paper you cite, either drop it or define it in the
sentence where it first appears.

## 2.3 [K18][K4] ★★ Expand every abbreviation on first use — **in each document part**

⚠⚠ The abstract and the body are read independently. A journal indexes the abstract on
its own.

```
Violation found:  SEES expanded in the abstract, undefined at its first body use
Fix:              expand in both
```

★ Same for HlyE, MFI, ELISA, STROBE.

## 2.4 [K12] ★★ Introduce the antigen before you use it

```
Violation found:  HlyE first used in §1.6; defined in §2.1, ~4,500 characters later
```

★ Definition goes where the reader first meets the word.

## 2.5 [K31] ★ Say which analysis you mean

```
❌ the companion analysis · the other paper
✅ the kinetics analysis of chapter 2  (first mention)
✅ the kinetics analysis                (thereafter)
```

⚠ In the integrated dissertation, "chapter 2" is the right anchor. ★ When the chapter
separates for journal submission, replace it with the citation.

## 2.6 [K5] ★ Define statistical shorthand before using it

```
"nominal"  →  define once: "an interval constructed to have 95% coverage"
```

---
---

# 3. Structure

## 3.1 [K15][K32] ★★ Standard informative headings

⚠ Kristen: *"use manuscript-style informative section headings."*

```
❌ "What remains" · "The seam this opens" · conversational headings
✅ noun phrases naming the content
✅ Kristen accepted: "Study population and sampling" · "Measurement and detection limits"
```

## 3.2 [G5.6] ★ No sub-headings in the Discussion

★ chapter 3 v2 complies. ⚠ Keep it that way.

## 3.3 [K25][G3.2] ★★★ IMRAD role separation

```
Introduction  ★ the existing method, what it assumes, why the assumption matters
Methods       ★ the new method only
Results       numbers, no interpretation
Discussion    interpretation, no new numbers
```

⚠ **See §0.2** — [K25] permits equations in the introduction when they describe the
existing method. That is not a violation.

## 3.4 [K19][G4.13] ★★★ No numbered or bulleted lists

⚠ Absolute. v1 had 15+; v2 has zero.

```
★ Convert every list to prose. Three items become one sentence with semicolons or
  a short paragraph.
```

## 3.5 [K20][K21] ★ Participant flow as a figure, per STROBE

```
❌ a table of exclusions
✅ a flow diagram
```

⚠ Kristen asked for it as Figure 1. ★ In chapter 3 it is Figure 2, behind the assumption
schematic. **That is an open question with her, not a settled violation.**

## 3.6 [D] ★★ Appendix numbering must not collide

```
Violation found:  "# Supporting material {.unnumbered}" with level-2 children
                  → children inherit the previous numbered section and render 4.1–4.4,
                    which reads as a subsection of "4 Discussion"
Fix:              number the parent
```

⚠ In the thesis this does not arise — the appendices are numbered top-level sections
there. ★ A reviewer must check which document it is looking at.

---
---

# 4. Sentences

## 4.1 [G1.4-Short] ★★ One idea per sentence

```
★ Target   under 40 words
⚠ chapter 3 v2: 74 sentences (17%) exceed it
```

**Longest offender and its fix:**
```
before  "The estimator has been validated against blood-culture surveillance in these
         same populations, where the rank ordering of seroincidence across countries
         agreed with the rank ordering of culture-confirmed incidence, and it is that
         agreement — not the level of any single estimate — that the field relies on
         when it uses seroincidence to compare places." (52 words)

after   "The estimator has been validated against blood-culture surveillance in these
         same populations. The rank ordering of seroincidence across countries agreed
         with the rank ordering of culture-confirmed incidence. It is that agreement,
         not the level of any single estimate, that the field relies on when it uses
         seroincidence to compare places."
```

## 4.2 [K22] ★ Complete sentences only

```
❌ "What remains."  ·  "Nothing else departs."
✅ "The final analysis sample comprised 2,613 participants."
```

## 4.3 [K29][K32] ★★★ No conversational asides

⚠ Kristen flagged four by name in v1: *"Nothing else departs"* · *"the change is not
cosmetic"* · *"the whole content of the chapter"* · *"that statement as a picture"*.

★★ **What they have in common**: each is vague or refers to the writing rather than to
the science. That is the test — **not** whether a sentence sounds informal.

```
❌ "the check costs seconds"                    → vague about cost
❌ "That comparison is not made anywhere"        → evasive passive
❌ "This matters more than a generic reassurance would."  → asserts its own importance
```

⚠⚠ **A concrete, checkable statement is not a violation even if it is vivid.**

```
✅ KEEP  "An interval can be widened by anyone. Multiplying it by two makes it wider
          without making it right."
```
★ This is specific, verifiable, and carries the chapter's argument. Flattening it loses
content.

## 4.4 [G4.5][G4.6] ★ Numerals

```
★ Numbers under ten in words, except with units
★ Never begin a sentence with a numeral
```

---
---

# 5. Objectivity

## 5.1 [G5.5] ★★★ No psychological or subjective framing

⚠ The guideline forbids "interesting", "surprising", "striking", and any statement of
the authors' reaction.

**Violations found in chapter 3 v2:**

| before | after |
|---|---|
| the **most interesting thing** this chapter found | the finding this chapter is least able to explain |
| the **honesty** of the uncertainty | the calibration of the interval |
| the **honest** reporting choice | the appropriate reporting choice |
| the kind that **does not announce itself** | the error is not self-evident in the output |

★★★ The first of these matters most. ⚠ Kristen said the same thing out loud about the
same result: *"I wouldn't say it's the most interesting finding, because it's not quite
a finding. We don't really know what it is yet."*

## 5.2 [G1.4-Specific] ★ No subjective intensifiers

```
❌ striking · remarkable · dramatic · surprisingly · notably · clearly
★ If the number is large, the number says so.
```

## 5.3 [D] ★★ State what you did not do, plainly

```
❌ "That comparison is not made anywhere in this chapter."
✅ "We did not make that comparison."
```

⚠ Passive avoidance reads as concealment. ★ Owning it reads as care.

---
---

# 6. Repetition

## 6.1 [G2.6] ★★★ Say it once

**Violation found:** detection-limit handling appears in §1.5, §2.1, and §4.2.

⚠ Kristen's comment [K24] was *"move it up"*, which is not *"put it in three places"*.

```
★ Keep    the general statement where the method is first described
★ Keep    the specific statement where equations need it
⚠ Delete  the third
```

## 6.2 [G2.6] ★★ The Discussion does not restate Results numbers

```
Violation found:  89.7% and 98.0% appear in Results and again in Discussion
                  46% appears three times
```

★ **Rule**: the Discussion may name a result. It may not re-quote its value unless the
comparison is new.

## 6.3 [D] ⚠ Cross-references replace repetition

```
❌ "for the reason given there"          ← ambiguous antecedent
✅ "@sec-ch3-methods"                     ← a reference resolves
```

★ **Test**: if "there", "above", or "as noted" cannot be replaced by a section anchor,
the sentence is unclear.

---
---

# 7. Equations and notation

## 7.1 [K14] ★★ Do not quote source text; summarize and show the equation

⚠ Kristen: *"summarize rather than quote directly, and show it as an equation."*

## 7.2 [K28] ★★★ Concept before notation

⚠ Kristen on the τ paragraph: *"hard to follow. Put the concept before the notation,
delete 'atom', and cut or move the t/τ and p/f symbol discussion."*

```
★ Order   what the quantity means  →  the symbol for it  →  the equation
❌        symbol  →  definition  →  meaning
```

## 7.3 [D] ★★ Notation must match across documents

```
★ chapter 3 uses τ, and explains why it departs from the software's T
⚠ Slides and the thesis must use τ too
```

★ **Test**: pick three symbols and grep every document that discusses them.

## 7.4 [D] ⚠ LaTeX renders

```
Violation found:  "$ ho_c$"  ← backslash lost from \rho
```

★ **Automated check**: search for `$ ` followed by a Greek letter name — `ho`, `ambda`,
`heta`, `au`, `igma`, `lpha`, `eta`, `elta` — with no preceding backslash.

---
---

# 8. Claims and attribution

## 8.1 [K1] ★★★ Do not claim validation against a method whose limits you are showing

⚠ Kristen: *"how can you validate against the existing method when we know its
limitations?"*

```
❌ "we validated the joint likelihood against the established estimator"
✅ "we verified that it reproduces the established estimator when only one biomarker
    is measured"
```

★★ **Same rule elsewhere**: a coverage figure computed with your own product form is not
a statement about the published package.

```
❌ "89.7% under the established estimator"
✅ "89.7% when the two isotype log-likelihoods are summed"
✅ "89.7% under the product form"
```

## 8.2 [K3] ★★ Assumptions before you relax them

```
★ Order   state the assumption  →  say you relax it  →  show what changes
```

## 8.3 [K6] ★ Define every performance term

```
"accuracy"  →  "relative root mean squared error"
"coverage"  →  "the proportion of intervals containing the true value"
```

## 8.4 [D] ⚠ Every number traces to a source in the same document

★ **Test**: every figure in an abstract, introduction, bridge, or discussion must appear
in the results or an appendix of the same document. A summary may not introduce a number.

## 8.5 [D] ★★★ Every external figure traces to a citation that actually contains it

⚠⚠ This one has bitten us.

```
Violation found:  "an estimated 200 million episodes of shigellosis annually
                   [@kotloff2013gems; @troeger2018]"

★ Kotloff 2013 is a seven-site case-control study — it cannot produce a global count
★ Khalil & Troeger 2018 reports deaths and per-1000 rates, not a global episode total
⚠ No cited source contains 200 million
```

★ **Rule**: before citing a figure, confirm the cited paper states that figure. A
plausible number in a review that cites the paper is not the same thing.

---
---

# 9. Abstract

## 9.1 [K7] ★ Name both diseases in the first sentence
## 9.2 [K9] ★★ The longitudinal model comes before the estimator that uses it
## 9.3 [K10] ★ Say how participants were selected
## 9.4 [K4] ★★ Expand abbreviations — the abstract stands alone
## 9.5 [D] ⚠ Length

```
Journal        under 300 words, or the journal's limit
Dissertation   ★ no limit, but two pages double-spaced is the practical ceiling
```

---
---

# 10. Automated checks a reviewer can run

★ These need no judgment. Run them first.

```
1  passive voice %          (was|were|is|are|been) + past participle ÷ sentences
2  "we" count               and whether any is unattributable  [G5.4]
3  first person singular    \bmy\b  \bI\b  \bmine\b            [D]
4  lists                    ^\s*[-*+]  or  ^\s*\d+[.)]          [G4.13]
5  sentence length          count over 40 words                 [G1.4]
6  subjective words         interesting|striking|surprising|remarkable|honest|
                            notably|clearly|obviously            [G5.5]
7  banned terms             transmission intensity|bank|arm S|cells|companion  [K2][K30][K31]
8  broken LaTeX             \$\s+(ho|ambda|heta|au|igma|lpha|eta|elta)\b        [D]
9  undefined abbreviations  every all-caps token; is it expanded before first use? [K18]
10 ambiguous reference      \b(there|above|as noted)\b without a nearby anchor  [D]
11 numeral at sentence start                                     [G4.6]
12 repeated result values   figures appearing in both Results and Discussion   [G2.6]
13 cross-reference integrity  every @ref resolves; every anchor is used
14 number provenance       every figure in a summary section appears in Results [D]
```

---
---

# 11. ⚠ What a reviewer must NOT flag

★ These have been decided. Flagging them wastes a round.

| Item | Why it stands |
|---|---|
| Equations in the introduction of chapter 3 | ★ [K25] required it. Overrides [G3.2] |
| "Multiplying it by two makes it wider without making it right" | ★★ Concrete and checkable. Carries the argument |
| Figure 2 as the flow diagram | ⚠ Open with Kristen, not a settled violation |
| "transmission year" in the seasonality sentence | ★ Names a period, not the estimand |
| τ rather than T | ★ Deliberate departure from the software, explained in the text |
| Appendix numbering in the thesis | ✅ Correct there; the issue is the standalone chapter only |

---
---

# 12. ★ Review output format

```
For each finding:

  [tag]  section · line
  rule   one sentence
  found  the exact text
  fix    the exact replacement
  grade  MUST / SHOULD / CONSIDER

Then, separately:

  OPINIONS  — things you would write differently that violate no rule
  QUESTIONS — things you could not resolve from the document
```

⚠ **Do not mix the three.** A MUST buried among opinions gets missed.

★★ And state what you could not check:

> *"I could not verify the burden figures because I did not read the cited papers."*

⚠ Silence about an unchecked item reads as a pass.
