# Manuscript Writing Standard — SERG dissertation chapters

**Version 3.** Rebuilt from primary sources on 2026-09-04, then corrected on the same
day by a review that used it. ⚠ Version 1 carried two numeric thresholds that appear in
neither source (§0.4); ⚠ version 2 promoted three targets that cannot hold at once
(§0.5).

---

## Sources, and what I actually read

```
[K##]  Kristen Aiemjoy's review comments on chapter3_v1_Kristen_reviewed.docx
       ★ 34 comments, read verbatim from word/comments.xml

[G#.#] Luby S, Southern D. The Pathway to Publishing: A Guide to Quantitative Writing
       in the Health Sciences. Springer, 2022. 188 pp.
       ★ Section numbers verified against the table of contents.
       ⚠ Sections quoted below were read in full. Sections cited without a quotation
         were located but not read end to end, and are marked ⚠ unread.

[D]    Derived by us. ★ Every one says what it was derived from and why.
```

⚠ **A rule with no tag is an opinion.** Opinions go in a separate list.

## ⚠ Scope — Kristen asked for both chapters

> *"please do the same for chapter 2 as well before i review"* [K33]

★★ Her comments were written on chapter 3, but she asked for chapter 2 to be revised
the same way **before she reviews it**. ⚠ So this standard governs both, and a chapter 2
review is not an extrapolation — it is what she asked for.

★ Where a rule does not transfer, §2.6 and §11 say so.

---

# 0. How to use this

## 0.1 Precedence

```
★★★ Kristen's explicit instruction wins over a general guideline principle.
```

**Worked example.** [G3.2] separates the roles of Introduction and Methods. [K25] says
*"previous method goes in background, new methods go in methods section."* Following
[K25] put four equations in chapter 3's introduction. ★ That is correct as it stands.

⚠ A reviewer must not flag it.

## 0.2 What a review does

```
★ Report        rule violated · location · the text · a concrete replacement
⚠ Do not        rewrite silently · improve prose that violates nothing ·
                report a measurement as a violation when the source sets no threshold
★★ Separate     "violates a rule" from "I would have written it differently"
```

## 0.3 ★★★ Measurements are not thresholds

⚠ The guideline names three numeric targets, all in [G8.18]. ★★ Only the first can be
treated as a target — see §0.5.

```
✅ [G8.18]  average words per sentence  < 25       ★ enforce
⚠ [G8.18]  Flesch Reading Ease         > 50       ★ report only
⚠ [G8.18]  Flesch-Kincaid grade level  16–18      ★ report only
```

★★ For anything else — passive voice percentage, longest sentence, count of a phrase —
**measure it and report the number, then point at the sentences.** ⚠ Do not call a
percentage a violation.

## 0.4 ⚠⚠ Two errors in version 1, corrected here

### ★ "passive voice under 25% of sentences"

⚠ **Not in either source.** [G5.3] gives no percentage at all — see §1.1, which now
quotes it. The number 25 was taken from [G8.18], where it is the target for **average
words per sentence**, a different quantity. Version 1 also stated that chapter 3's
measured 23.4% *"is the standard, not an aspiration"*, ★ which converted an observation
into a rule.

### ★ "sentences under 40 words"

⚠ **Not in either source.** [G1.4-Short] says *"Use short sentences containing only one
idea in each. Split complex sentences."* — no count. [G8.18] targets an **average** below
25, not a per-sentence cap.

★★ Both are removed. What replaces them is in §1.1 and §4.1.

## 0.4a ⚠ Three errors in version 2, corrected here

★ All three were found by a review that used version 2 and checked it against the
sources rather than against my summary of what had changed.

```
① I said version 2 added four Kristen comments. ★ It added eight —
   K8, K11, K13, K23, K26, K27, K33, K34. My check matched verbatim strings and
   missed the two long replacement passages and two that version 1 had paraphrased.

② I said it added twelve guideline sections. ★ Thirteen — I omitted
   [G1.4-Sequential], which §6.3 rests on.

③ I said §11 lists eighteen items. ★ It lists fourteen. Eighteen is the count of
   checks in §10. ⚠ The wrong number reached a review instruction.
```

★★ **A count in this document is a fact about this document.** ⚠ Verify it against the
file rather than against a change note.

## 0.5 ⚠⚠ The three [G8.18] targets cannot hold at once

★★★ This was found by a review that applied version 2, not by reading the guideline.

```
To reach Flesch-Kincaid grade 16 at 23 words per sentence, prose needs about
1.85 syllables per word. ⚠ At that density Flesch Reading Ease falls below 25 —
half of target 2.

Conversely, prose that holds Reading Ease above 50 at 23 words per sentence tops
out near grade 12.
```

★★ So targets 1 and 2 are jointly reachable; ⚠ target 3 is not reachable with them.

★ Chapter 3 measures 22.9 words per sentence, Reading Ease 38–49, and grade 12–14. It
sits at the joint optimum of targets 1 and 2. **Its grade level is a consequence of
meeting the other two, not a defect.**

### ★ How to use them

```
1  average words per sentence < 25   ★ a target. Report it and treat a miss as a finding
2  Flesch Reading Ease                ⚠ report the number. Do not grade
3  Flesch-Kincaid grade level         ⚠ report the number. Do not grade
```

⚠ And say which syllable counter produced the reading scores. ★ A regex heuristic and a
hyphenation library can land on opposite sides of 50 for the same text.

---
---

# 1. Voice and person

## 1.1 [K16][K32][G5.3][G1.4-Strong] ★★ Active voice — with the exceptions the source names

⚠ Kristen's verdict on v1, verbatim:

> *"this reads like train of thought. try to use clear active voice scientific writing"* [K16]

> *"I would like you to use ACTIVE scientific voice and a general-to-specific structure,
> consistent with The Pathway to Publishing"* [K32]

★ And [G1.4-Strong]:

> *"Use the verb as the center of gravity of your sentence. If the verb is weak, the
> sentence is weak. For example, instead of 'we did an interview,' write 'we
> interviewed.' Use active voice instead of passive."*

### ⚠⚠ But [G5.3] explicitly permits passive, and names when

> *"Although passive voice is used in many scientific articles, **especially in the
> methods section**, active voice is increasingly common."*

> *"**When to Use Passive Voice** — It is not always an error to use passive voice.
> Passive voice is particularly useful, even recommended, in two situations:*
> *1. When it is more important to draw our attention to the person or thing acted upon*
> *2. When the actor in the situation is not important: Passive voice is especially
> helpful in scientific or technical writing or laboratory reports where the process or
> principle being described is of ultimate importance."*

### ★★★ So the rule is

```
★ Where you can say the same thing either way, choose active.       [G5.3]
★ Methods sections legitimately carry more passive.                 [G5.3]
⚠ Report the percentage. ★ Do not call it a violation.              [D, from §0.3]
★ Point at specific sentences where active would have served.
```

★ The test [G5.3] gives is not a ratio:

> *"its imprecisions risks conveying that the authors are **unwilling to specify who took
> the action**."*

⚠ **That** is what to look for — a passive that hides an actor who matters.

### ★★ Kristen wrote a replacement sentence, verbatim [K26]

> *"again use active voice throughout the manuscript: **We used the estimator implemented
> in serocalculator and followed the notation in its methodology vignette so readers can
> compare each step directly.**"*

⚠ Note what happened to that sentence afterwards. [K1] forbids claiming validation
against a method whose limits you are showing, and the chapter's own body says the
comparator is internal. So the verb had to change again — *"We describe the estimator
implemented in serocalculator…"* — because §1.5 exists to describe the prior method
rather than to report running it.

★★★ **Two of Kristen's comments can pull in different directions.** When they do, the
one about the claim wins over the one about the voice.

### ★ Reference measurements

```
chapter 3 v1  37.2%   ⚠ Kristen: "reads like train of thought"
chapter 3 v2  24.0%   ★ Kristen has not seen it
chapter 2     26.8%
```

⚠ ★ These are context, not targets.

## 1.2 [G5.4] ★★★ "We" must be attributable

★ Verbatim:

> *"Work that is conducted by fieldworkers or other members of the team who are not on
> the author line should not be attributed to the authors."*
>
> *✗ We interviewed households at baseline* → *✓ Trained enumerators interviewed
> households at baseline*

★★ In these chapters that means:

```
✅ "We evaluated the likelihood across a grid"        analysis, authors did it
✅ "Antibody responses were measured by ELISA"        laboratory, staff did it
⚠ "We measured HlyE IgG and IgA"                      ★ violation
```

## 1.3 [D] ★★ No first person singular

⚠ Derived from [G5.4]: these are multi-author manuscripts, and "my" names one.

```
Found in chapter 3 v2:  "a second version of my own code"
Fixed to:               "a second implementation written for this analysis"
```

★ Scan for `\bmy\b`, `\bI\b`, `\bmine\b` outside quotations.

---
---

# 2. Terminology

## 2.1 [K2] ★★ Name the quantity you estimate

⚠ Verbatim: *"we measure infection rate not transmission"*

```
❌ transmission · transmission intensity · force of infection
✅ seroincidence · infection rate
```

★ "transmission year" naming a period is acceptable — it is not the estimand.

## 2.2 [K30][K34][G6.5] ★★★ Do not invent labels

⚠ Kristen, twice:

> *"I still don't understand why we are using the new term 'bank' can we use standard
> terminology or else the concept/terminology needs to be clearly introduced at first
> use"* [K30]

> *"what is arm S? rather than infecting new nomenclature can we describe what each
> thing is"* [K34]

> *"use epidemiology terms 'study sites' not 'cells'"* [K22]

★ And [G6.5] names the mechanism:

> *"Study teams commonly develop some study-specific vocabulary, for example, Group A
> and Group B and Phase 1 and Phase 2. The study team becomes so familiar with these
> labels and their underlying characteristics that they use these labels in everyday
> conversation within the study team."*

```
❌ bank · parameter bank · arm S · cells
✅ posterior parameter draws · paired draws · independent draws · study sites
```

★ **Test**: if a term appears in no paper you cite, either drop it or define it in the
sentence where it first appears.

## 2.3 [K4][K18][G4.1] ★★ Expand every abbreviation on first use — in each document part

⚠ Kristen: *"this acronym will be for abstract readers"* [K4] and *"sell out accronym at
fist use"* [K18]

★★ The abstract and the body are read independently; a journal indexes the abstract on
its own. **Expand in both.**

★ Applies to: SEES, SEAP, HlyE, MFI, ELISA, MCSE, STROBE, LPS.

## 2.4 [K12] ★★ Introduce the antigen before you use it

⚠ Verbatim: *"introduce the antigens / hlye"*

★ The definition goes where the reader first meets the word, not four sections later.

## 2.5 [K31] ★ Say which analysis you mean

⚠ Verbatim: *"'companion anlalysis' is vague do you mean the longitudinal model you
estimate in chapter 2?"*

```
❌ the companion analysis · the other paper
✅ the kinetics analysis of chapter 2   (first mention)
✅ the kinetics analysis                 (thereafter)
```

## 2.6 [K6] ★ Define statistical shorthand before using it

⚠ Verbatim: *"how is accuracy defined here?"*

⚠⚠ **This one is journal-dependent.** Kristen's comment is on a manuscript for a clinical
readership. A term that needs defining for *Lancet Microbe* may be standard in
*Statistics in Medicine*.

```
★ chapter 3 (Lancet Microbe)      ★★ define: nominal · accuracy · coverage
                                  ⚠ Kristen's comment was written ON this manuscript
★ chapter 2 (Statistics in Med)   ⚠ these are standard. Defining them is not required
```

⚠⚠ **Read the direction carefully.** [K6] was written on chapter 3, for a clinical
readership. In chapter 3 the definitions **are required**. In chapter 2 they are not.

★ A reviewer working on chapter 3 who finds `nominal` undefined has found a violation.
A reviewer working on chapter 2 who finds the same thing has not.

## 2.7 [G5.9] ★★★ One term per object, throughout

★ Verbatim:

> *"To avoid mind-numbing repetition, authors commonly vary word choice and style
> throughout the manuscript. Although such variation can engage readers, **if it is
> applied to scientific terms, it risks confusing readers.**"*
>
> *✓ Define the term injury … in the methods section. Use the term injury consistently
> throughout the manuscript.*

★★ This is the rule behind chapter 2's null/alternative problem: nine phrasings for one
condition, none of them the standard pair.

```
★ Define the pair once, where the design is stated.
★ Use it in prose thereafter.
⚠ Figure and table labels may stay descriptive — a legend has no preceding sentence
  to carry the definition.
```

---
---

# 3. Structure

## 3.1 [K15][K32] ★★ Informative section headings

⚠ Verbatim: *"use manuscript/informative section headings"* [K15]

★ And [K32] gives examples of both kinds:

> *"Several current headings read like a train of thought rather than section headings,
> for example 'Two more arrangements of the same integral,' 'Two assumptions, not one,'
> and 'The shortfall is width, not location.' Please use headings that tell the reader
> what the section is about in a formal, scientific way, such as **'Alternative
> likelihood formulations,' 'Decomposition of independence assumptions,' or 'Coverage
> and interval calibration.'**"*

## 3.2 [G5.6] ★ No sub-headings in the Discussion

⚠ unread. Section located at p. 71.

## 3.3 [K25][G3.2] ★★★ IMRAD role separation

⚠ Kristen: *"previous method goes in background, new methods go in methods section"*

```
Introduction  ★ the existing method, what it assumes, why the assumption matters
Methods       ★ the new method only
Results       numbers, no interpretation
Discussion    interpretation, no new numbers
```

★ **See §0.1** — [K25] permits equations in the introduction when they describe the
existing method.

## 3.4 [K19][G4.13] ★★★ No numbered or bulleted lists

⚠ Verbatim: *"dont use numbered/bulleted list"*

★ Absolute. Convert every list to prose.

## 3.5 [K20][K21] ★ Participant flow as a figure, per STROBE

⚠ Verbatim: *"use a flow diagram instead and this can be included as figure 1 of the
manuscript (as per STROBE guidelines)"* [K20] and *"replace with flow diagram / figure
1"* [K21]

```
❌ a table of exclusions
✅ a flow diagram
```

⚠ Chapter 2 has it as Figure 1. Chapter 3 has it as Figure 2, behind an assumption
schematic. ★ That difference is open with Kristen, not a settled violation.

## 3.6 [G3.8] ★★★ No new data in the Discussion

★ Verbatim:

> *"The role of the discussion is to tell the reader what the authors believe the results
> mean. It is a violation of the standard IMRAD format …"*

⚠ unread beyond this sentence.

★★ This is the rule behind chapter 3's age-stratified correlation: the Discussion
reported four values that appeared nowhere in Results.

## 3.7 [D] ★★ Appendix numbering must not collide

⚠ Derived from a rendering fault, not from either source.

```
"# Supporting material {.unnumbered}" with level-2 children
  → children inherit the previous numbered section and render 4.1–4.4,
    reading as subsections of "4 Discussion"
★ Fix: number the parent
```

⚠ In the thesis this does not arise — the appendices are numbered top-level sections
there. ★ Check which document you are looking at.

---
---

# 4. Sentences

## 4.1 [G1.4-Short][G8.18] ★★ One idea per sentence

★ Verbatim [G1.4-Short]:

> *"Use short sentences containing only one idea in each. Split complex sentences. Cut
> unnecessary information elements and only include those data that relate to the point
> of your paper. … Remember, 'if it's only nice to know, it ought to go.'"*

★★★ And the one measurable target, verbatim [G8.18]:

> *"1. Average words per sentence should be <25. Strive to be concise."*

```
✅ Measure    average words per sentence — ★ target < 25
⚠ Do not     apply a per-sentence word cap. Neither source sets one.
★ Report     the longest sentences and whether each carries more than one idea
```

⚠ **A 45-word sentence carrying one idea is not a violation.** A 30-word sentence
carrying three is.

## 4.2 [G8.18] ★ Readability, measurable

> *"2. Readability: (a) Flesch Reading Ease on a scale of 1–100. … Strive for >50.
> (a) Flesch-Kincaid Grade level … Target a grade level of 16–18."*

⚠ These are the only other numbers either source gives. ★ Report them.

## 4.3 [K23] ★ Complete sentences only

⚠ Verbatim: *"awkward transition/incomplete sentence"*

```
❌ "What remains."  ·  "Nothing else departs."
```

## 4.4 [K32] ★★★ No conversational asides

⚠ Verbatim, and this is the fullest statement of the problem Kristen had:

> *"In many sections, the chapter reads more like an internal technical memo than a
> dissertation chapter written for a scientific audience. Phrases such as **'nothing else
> departs,' 'the change is not cosmetic,' 'the whole content of the chapter,' and 'that
> statement as a picture'** feel conversational or defensive, and they do not help orient
> the reader."*

★★★ **What those four share**: each is vague, or refers to the writing rather than to
the science. ⚠ That is the test — **not** whether a sentence sounds informal.

```
❌ "the check costs seconds"                             vague about cost
❌ "That comparison is not made anywhere in this chapter"  evasive passive
❌ "This matters more than a generic reassurance would"    asserts its own importance
✅ "Multiplying it by two makes it wider without making it right"
   ★ specific, checkable, carries the argument. Keep it.
```

## 4.5 [K32] ★★ Paragraphs open with the message

★ Verbatim:

> *"Each paragraph should generally begin with the main message or interpretation,
> followed by the supporting technical detail. Right now, many paragraphs begin with
> notation, implementation details, or algebraic distinctions, and the main point only
> becomes clear later."*

## 4.6 [G4.5][G4.6] ★ Numerals

⚠ unread. [G4.5] "Failure to Spell Out an Isolated Numeral < 10", [G4.6] "Starting a
Sentence with a Numeral".

---
---

# 5. Objectivity

## 5.1 [G5.5] ★★★ No psychological or subjective framing

★ Verbatim:

> *"What interests or surprises people varies and often depends upon their personal
> experiences … when you are writing a scientific manuscript, you should focus on the
> ideas relevant to the issues examined in your study and the consistency of ideas and
> theories with available evidence."*
>
> *✗ We were surprised to find that people admitted to using alcohol …*

★★ Found and fixed in chapter 3 v2:

| before | after |
|---|---|
| the **most interesting thing** this chapter found | the finding this chapter is least able to explain |
| the **honesty** of the uncertainty | the calibration of the interval |
| the **honest** reporting choice | the appropriate reporting choice |
| the kind that **does not announce itself** | the error is not self-evident in the output |

⚠⚠ The first mattered most. Kristen said the same thing aloud about the same result:
*"I wouldn't say it's the most interesting finding, because it's not quite a finding."*

## 5.2 [G1.4-Specific][G6.3] ★ No subjective intensifiers

★ Verbatim [G1.4-Specific]:

> *"Don't use qualifiers, which are imprecise and judgmental. Avoid words such as 'very,'
> 'rather,' or 'much.' Choose adjectives carefully. Don't use adjectives that imply
> subjectivity and/or emotion, for example, 'It was a very large outbreak.' **What does
> very mean? How big is large?**"*

```
❌ striking · remarkable · dramatic · surprisingly · notably · clearly · greatly · markedly
★ If the number is large, the number says so.
```

## 5.3 [D] ★★ State what you did not do, plainly

⚠ Derived from [G5.3]'s test — passive that hides an actor.

```
❌ "That comparison is not made anywhere in this chapter."
✅ "We did not make that comparison."
```

---
---

# 6. Repetition

## 6.1 [G2.6] ★★★ Say it once

⚠ [G2.6] "Repeating Information", p. 39. ⚠ unread in full, but [G2.3.1] gives the
adjacent principle verbatim:

> *"Best practice is to refer to a prior article that provided details and then offer a
> succinct summary."*

★★ Applied within a document: **refer to the section that has the detail, then summarize.**

⚠ Found in chapter 3 v2: detection-limit handling appeared in §1.5, §2.1 and §4.2.
[K24] said *"this is part of the general method, can be moved above"* — ★ which is not
"put it in three places."

## 6.2 [G3.8][D] ★★ The Discussion may name a result, not re-quote its value

⚠ Derived from [G3.8]. ★ The exception: an opening restatement at the head of a
Discussion is conventional and helps a reader who has not read the Results.

## 6.3 [D] ★★ Cross-references replace backward pointers

⚠ Derived from [G1.4-Sequential]:

> *"Take the reader by the hand through the sequence of thoughts, step-by-step, without
> any leaps or missing links."*

```
❌ "for the reason given there"          ambiguous antecedent
✅ "@sec-ch3-methods"                     a reference resolves
```

★ **Test**: if "there", "above", or "as noted" cannot be replaced by a section anchor,
the sentence is unclear. ⚠ A pointer to a few lines up inside the same subsection is
ordinary signposting and is fine.

---
---

# 7. Equations and notation

## 7.1 [K14] ★★ Summarize the source; show the equation

⚠ Verbatim: *"i dont think you need to quote the description verbatim - instead you can
summarize and possibly show the math"*

## 7.2 [K28] ★★★ Concept before notation

⚠ Verbatim, and this is the most detailed instruction Kristen gave:

> *"I found this paragraph hard to understand, and I think the issue is **not just wording
> but structure.** The reader needs a clear explanation of the modeling idea **before the
> notation.** Please revise to explain first that, in the serosurvey, the infection date
> is unknown and some participants may never have been infected, so time since infection
> is treated as a latent variable and integrated over. Avoid technical shorthand like
> 'an atom at never infected' unless you define it in plain language. The discussion of
> symbol choices (t vs. τ, p vs. f) feels distracting here and **should be shortened or
> moved out of the main text.**"*

```
★ Order   what the quantity means  →  the symbol for it  →  the equation
❌        symbol  →  definition  →  meaning
```

⚠ Note *"shortened **or** moved"* — shortening satisfies it.

### ★★ And she wrote the replacement herself [K27]

> *"The model assumes that infections occur randomly over a participant's lifetime
> according to a Poisson process with rate λ, the seroconversion rate in infections per
> person-year. For a participant of age a, the model first distinguishes between
> participants who have never been infected and those who have been infected at least
> once. Among participants who have been infected, τ denotes the time since the most
> recent infection and follows a truncated exponential distribution."*

★★★ Read the order she used: **the assumption in words → the rate and its units → the
two cases → then τ.** The symbol arrives last and only after the reader knows what it
measures.

⚠ ★ This paragraph is the clearest single statement of what [K28] and [K32] are asking
for. When a passage is hard to follow, compare it to this one.

## 7.3 [D] ★★ Notation must match across documents

⚠ Derived. ★ chapter 3 uses τ and explains why it departs from the software's T; the
slides and the thesis must use τ too.

## 7.4 [D] ⚠ LaTeX renders

⚠ Derived from a defect found in chapter 3 v2: `$ ho_c$`, a lost backslash.

★ Search for `$ ` followed by a Greek letter name with no preceding backslash — `ho`,
`ambda`, `heta`, `au`, `igma`, `lpha`, `eta`, `elta`.

---
---

# 8. Claims and attribution

## 8.1 [K1] ★★★ Do not claim validation against a method whose limits you are showing

⚠ Verbatim: *"how can it be validated against current use when we know current use has
limitations?"*

```
❌ "we validated the joint likelihood against the established estimator"
✅ "we verified that it reproduces the established estimator when only one biomarker
    is measured"
```

★★ Same rule elsewhere: a coverage figure computed with your own product form is not a
statement about the published package.

## 8.2 [K3] ★★ Assumptions before you relax them

⚠ Verbatim: *"the assumptions need to be described first before sayign that they are
relaxed"*

## 8.3 [K9] ★ The longitudinal model comes first

⚠ Verbatim: *"the base longitudinal model needs to be describe first"*

## 8.3a [K13] ★★ Describe the prior method's own combination rule

⚠ Verbatim: *"now describe how serocalculator combines information from multiple
biomarkers"*

★★ This is the comment that put four equations in chapter 3's introduction. Kristen did
not ask for a citation to the software; she asked for the combination rule to be
**written out**, so the reader can see what the new method changes.

★ Paired with [K3] and [K25], the order is: state the existing rule → state what it
assumes → then relax it.

## 8.4 [G2.2.1][G2.3.1] ★★★ Cite the primary source, and check it says what you say

★ Verbatim [G2.3.1]:

> *"In this form of the error, the author cites an article that cites the original
> observation. Standard scientific practice is to cite the primary observation. It is a
> flagrant error if you cite an article that makes a similar point to the argument you
> want to make in your article, and the article that you are citing perhaps, in its
> introduction, cites the primary articles."*

⚠⚠ **This has bitten us.**

```
Found:  "an estimated 200 million episodes of shigellosis annually
         [@kotloff2013gems; @troeger2018]"

★ Kotloff 2013 is a seven-site case-control study — it cannot produce a global count
★ Khalil & Troeger 2018 reports deaths and per-1000 rates, not a global episode total
⚠ No cited source contains 200 million
```

★ **Rule**: before citing a figure, confirm the cited paper states that figure. A
plausible number in a review that cites the paper is a secondary source.

## 8.5 [D] ⚠ Every number traces to a source in the same document

⚠ Derived from [G3.8].

★ **Test**: every figure in an abstract, introduction, bridge, or discussion must appear
in the results or an appendix of the same document.

### ⚠⚠ Exception — where the generating code lives

★★★ **A missing writer script in this repository is not a finding.**

Tables and figures are produced by scripts that run on the **Mercury compute server**,
which are not all mirrored to GitHub. The committed CSV under `tab/` is the citable
artifact and the analysis is reproducible.

```
★ In scope     a number in the prose with no CSV and no table behind it
⚠ Out of scope a CSV that exists and is committed, but whose writer script is not here
```

## 8.6 [G7.3][G7.4] ★★ Decimal places, and where rounding changes a claim

⚠ unread. ★ But the distinction we needed is [D]:

```
✅ Descriptive rounding   median 0.639 → "0.64"        ★ fine
⚠ Boundary claims        lower bound 0.525 → "0.53"   ★★ changes the fact
                         "no smaller than 0.53" is false when 0.525 is in the table
```

---
---

# 9. Abstract

⚠ Kristen wrote replacement text for two abstract passages. ★ Both are quoted here
because the wording is hers, and changing it is a question for her rather than a
decision for us.

## 9.1 [K7] ★ Name both diseases

⚠ Verbatim: *"enteric fever should be mentioned somewhere"*

★ Note what this does **not** say: it does not ask for typhoid and paratyphoid to be
listed separately.

## 9.2 [K8] ★★ Kristen's replacement opening, verbatim

> *"Cross-sectional serosurveys can estimate the seroincidence rate of enteric fever by
> leveraging antibody decay dynamics estimated from confirmed cases. These methods
> recover the same rank-order burden as blood-culture-based clinical incidence
> estimates, while requiring a fraction of the time, cost, and infrastructure. The
> established estimator combines biomarkers by adding their log-likelihoods, which
> treats one participant's two measurements as though they came from two independently
> infected people."*

## 9.3 [K11] ★★ Kristen's replacement closing, verbatim

> *"Seroincidence remains a useful, lower-cost approach for ranking enteric fever burden
> and supporting vaccine-prioritization decisions, but when multiple biomarkers are
> combined, the current method likely reports uncertainty that is too narrow."*

## 9.4 [K10] ★ Say how participants were selected

⚠ Verbatim: *"randomly selected population"*

## 9.5 [K4] ★★ The abstract stands alone — expand abbreviations there too

## 9.6 [G8.8][D] ⚠ Length

```
Journal        the journal's limit — Lancet Microbe and Stat Med differ
Dissertation   ★ no limit. Two pages double-spaced is a practical ceiling
```

---
---

# 10. Automated checks

★ These need no judgment. ⚠ Numbers 1, 2 and 3 have sourced targets; **the rest report
counts and nothing more.**

```
★ SOURCED TARGET  [G8.18] — one only, see §0.5
 1  average words per sentence        ★ target < 25

★ REPORTED, NOT GRADED  [G8.18]
 2  Flesch Reading Ease               ⚠ report the number and the syllable counter used
 3  Flesch-Kincaid grade level        ⚠ report the number. Unreachable alongside 1 and 2

★ COUNTS — report, do not grade
 4  passive voice %                   (was|were|is|are|been) + past participle
                                      ⚠ no threshold exists. Point at sentences where
                                      an actor is hidden [G5.3]
 5  longest sentences                 report the top ten and whether each has one idea
 6  "we" instances                    and whether any is unattributable  [G5.4]
 7  first person singular             \bmy\b  \bI\b  \bmine\b
 8  lists                             ^\s*[-*+]  or  ^\s*\d+[.)]   [G4.13] — zero
 9  subjective words                  interesting|striking|surprising|remarkable|honest|
                                      notably|clearly|greatly|markedly|very|rather|much
10  retired terms                     transmission intensity|bank|arm S|cells|companion
11  broken LaTeX                      \$\s+(ho|ambda|heta|au|igma|lpha|eta|elta)\b
12  undefined abbreviations           every all-caps token; expanded before first use?
                                      ⚠ check the abstract separately  [K4]
13  ambiguous back-references         \b(there|above|as noted)\b without a nearby anchor
14  numeral at sentence start         [G4.6]
15  repeated result values            figures in both Results and Discussion  [G3.8]
                                      ⚠ strip image attributes first — {width=90%} is
                                      markup, not prose
16  cross-reference integrity         every @ref resolves; every equation label is cited
17  number provenance                 every figure in a summary section appears in
                                      Results  ⚠ stop at the CSV (§8.5)
18  one term per object               list every phrasing used for each key condition
                                      [G5.9]
```

## ⚠ On sentence segmentation

★ A `.qmd` wraps sentences across lines. **A line-based splitter overcounts sentences and
splits long ones**, which understates both the average length and the passive fraction.

```
★ Correct rule  a prose line ending in terminal punctuation closes a sentence;
                otherwise it continues to the next line
⚠ Skip          code fences, YAML, #| options, table rows, equation displays
```

⚠ Chapter 3 measured 667 sentences under the naive rule and 532 under the correct one.

---
---

# 11. ⚠ What a reviewer must NOT flag

★ These have been decided. Flagging any of them wastes a round.

| Item | Why it stands |
|---|---|
| Equations in the introduction of chapter 3 | ★ [K25] required it. Overrides [G3.2] |
| "Multiplying it by two makes it wider without making it right" | ★★ Concrete and checkable. Carries the argument |
| Figure 2 as chapter 3's flow diagram | ⚠ Open with Kristen, not a settled violation |
| "transmission year" in the seasonality sentence | ★ Names a period, not the estimand |
| τ rather than T | ★ Deliberate departure from the software, explained in the text |
| Appendix numbering in the thesis | ✅ Correct there; the issue is the standalone chapter |
| A committed CSV whose writer script is not in this repository | ★★ Expected. See §8.5 |
| "A script held under version control produces every quantity we report" | ★ Accurate. The scripts exist and are versioned, just not all on GitHub |
| `nominal`, `accuracy`, `coverage` undefined in chapter 2 | ★★ Standard in Statistics in Medicine. [K6] was written for a clinical readership |
| The Discussion's opening restatement of a headline result | ★ Conventional. §6.2 |
| "95%" as a nominal coverage level | ★ Not a result value |
| Section anchors defined but never referenced | ★ Deliberate. Equation labels are different — those consume a number |
| Forward signposting inside a subsection | ★ "The forms below", "Everything below" — ordinary, not §6.3 |
| chapter 2 names SEAP, chapter 3 does not | ★ Follows the two designs. To be checked once both sit in the dissertation |

---
---

# 12. ★ Review output format

```
For each finding:

  [tag]  section · line
  rule   one sentence, ★ with the source's own words if short enough
  found  the exact text
  fix    the exact replacement
  grade  MUST / SHOULD / CONSIDER

Then, separately:

  MEASUREMENTS — the numbers from §10, with the three sourced targets marked
  OPINIONS     — things you would write differently that violate no rule
  QUESTIONS    — things you could not resolve from the document
```

⚠ **Do not mix them.** A MUST buried among opinions gets missed.

★★ And state what you could not check:

> *"I could not verify the burden figures because I did not read the cited papers."*

⚠ Silence about an unchecked item reads as a pass.

---
---

# 13. ⚠ Rules whose source I have not read in full

★ Recorded so nobody mistakes a located section for a read one.

```
[G3.2]  Confusing the Role of Introduction, Methods, Results and Discussion   p.42
[G3.8]  Presenting New Data in the Discussion — ★ first sentence read only    p.55
[G4.1]  Using Nonstandard Acronyms                                            p.65
[G4.5]  Failure to Spell Out an Isolated Numeral < 10                         p.68
[G4.6]  Starting a Sentence with a Numeral                                    p.68
[G4.13] Using Bulleted Lists Rather Than Sentences                            p.72
[G5.6]  Using Excessive Subheadings in the Discussion                         p.71
[G6.3]  Using Adjectives and Qualifiers                                       p.74
[G6.5]  Using Nondescriptive Numeric or Alphabetical Labels — ★ opening read  p.83
[G7.3]  Using Too Many Decimal Places                                         p.98
[G7.4]  Using Too Few Decimal Places                                          p.98
[G8.8]  Exceeding the Journal Word Limit                                      p.112
```

⚠ If a review turns on one of these, ★ read the section before deciding.
