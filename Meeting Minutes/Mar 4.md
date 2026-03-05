# Agenda - March 4

## Meeting Objectives
* Go over changes to our rarefaction methodology
* Review the rubric and comments for our proposal.
* GitHub reorganization

### Questions to ask:
*  Wilcox vs T-test (can we assume parametric distribution for Aim 1? In other words is shannon classification here. We saw differences in significance w/ T-test being significnat but Wilcox not)
*  Rarefaction depth and changing it for aim 3? (we dropped only 2 samples, but both from female intestinal group which already only had 3)
*  Issues with statistical power of only 3 samples?
*  Linear model vs 2 way ANOVA for parametric multiple predictors?
*  Visuals?



### 1. Rarefaction Decision
* **Rarefaction Impact:** Our original plan involved a rarefaction depth of 10,222 reads/sample.
* **Sample Loss:** This depth previously excluded 2 intestinal-type samples that were below the read count.
* **Context:** Those 2 samples represent 2 out of our 3 total female intestinal samples, and losing them would compromise our ability to model sex differences significantly.
* **Our Proposed Decision:** We will continue with the rarefaction depth of **10,222 reads/sample** for our original aim 1 (now aim 1a), but for the gender-based aim, which was originally aim 3 (now aim 1b), we will do a rarefaction depth of **8386 reads/samples**, so that we keep all of the female-inestinal samples.
* **Any Feedback?** 

### 2. Combining Previous Aims 1 and 3
* We decided to combine these aims as we talked about in our last meeting, so our 3 new aims are:
* **Aim 1:** Investigating the Microbiome Composition in Intestinal vs Diffuse Gastric Cancer and Evaluating the Influence of Biological Sex
* **Aim 2:** Identifying the Specific Bacteria that Distinguish Intestinal vs Diffuse Cancers
* **Aim 3:** Functional Analysis and Profiling of Microbial Communities Within Subtypes

### 3. Addressing Proposal Feedback
* **Title Revision:** We have to update our current title to better reflect the functional analysis we are going to perform.
* **Introduction & Background:** We need to expand upon existing studies to ensure the research gap is explicitly clear.
* **Experimental Aims:** Our aims don't have any specific analysis, so we have to add the ones we're going to use for our research question.
* **Proposed Approach Table:** We have to go over our proposed approach table and make sure the statistical analysis/tools and R functions align logically with our aims.

### 4. Reorganized our GitHub
* We also reorganized our GitHub and created individual files for:
  * Our meeting agenda/notes (we split it by day)
  * Our lab notebook (we split it by sub-aims)
  * Our final script (including QIIME2 and R aims)
* All of the files are accessible from our README file, in (hopefully) a clean and organized manner
* **Any Feedback?**



# Meeting Minutes - March 4

- Proposal feedback
  - Main revisions: expand references in the introduction; in the aims, clearly state the analyses we will run.
  - Workshop points: unclear if points are individual vs team-based
  - Plan: do the proposal revision to recover the additional points.

## Key decisions
- Use Wilcoxon (non-parametric) instead of t-test for Shannon/alpha-diversity comparisons
  - Reason: big disagreement between t-test and Wilcoxon suggests non-normality; Shannon metrics are not safe to assume normal.
  - Example discussed: t-test p = 0.01315 vs Wilcoxon p = 0.2219 
- Rarefaction depth must be consistent across aims.
  - Decision: rarefy at 8,000 reads for all aims.
  - Reason: original depth (10,000) caused Aim 3 to collapse to ~1 female sample; lowering to 8,000 keeps the female samples.
  - Rarefaction curve appears near plateau by ~8,000 (so 8,000 is likely adequate coverage).
  - Lowest sample depth discussed was ~8,260–8,386 (confirm exact minimum from table).

## Results reviewed (current figures)

### Alpha diversity (Shannon) – subtype × gender (cross-subtype within each sex)
- Current visualization is confusing (dot indicates median/mean); change to boxplots.
  - Make two separate plots: one for males, one for females (same axes, same formatting).
- Pattern observed:
  - Males: intestinal vs diffuse shows a clear/significant difference (direction noted as intestinal higher Shannon than diffuse in discussion; confirm in plot/stats).
  - Females: trend appears opposite (significant decrease discussed), but statistical power is low (only ~3 female samples)
- Keep this figure (it is interesting and supports a sex-dependent pattern), but clearly label sample sizes (n) on the plot.

### Beta diversity – Aim 1 (diffuse vs intestinal,)
- Bray-Curtis PCoA:
  - Clear separation between diffuse and intestinal.
  - PERMANOVA p-value reported in meeting: p = 0.0018 (confirm R² from output).
- UniFrac:
  - Unweighted UniFrac: not significant (p reported ~0.38; confirm).
  - Weighted UniFrac: not significant (confirm exact p).
- Interpretation discussed:
  - Bray-Curtis detects differences driven by community composition/abundance.
  - Lack of UniFrac significance suggests limited evidence for phylogeny-based separation between subtypes (at least under current filtering/rarefaction).
