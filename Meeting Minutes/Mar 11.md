# Agenda - March 11

## 1. Aim 1 – Revisions Presentation

Alpha Diversity (Shannon)
- Boxplots, separated by sex (males vs females).
- Wilcoxon p-values and sample sizes labelled.
- Discuss sex-dependent patterns and statistical power.

Beta Diversity (Bray-Curtis & UniFrac)
- Bray-Curtis PCoA plots: subtype separation, subtype × sex patterns.
- PERMANOVA results: p-values + R².
- UniFrac results (weighted/unweighted).

Questions:
- Are figures clear and publication-ready?
- Any suggestions on labelling, axes, or plot layout?
- Cross-subtype comparisons as discussed previously? Male vs female within each subtype

## 2. Aim 2 – Key Figures Review
- Heatmap: Displaying relative abundances of key taxa.
- Volcano Plot: Highlighting significant taxa between subtypes.
- Bar Plot: Showing composition differences for top taxa.

Feedback:
- Are plots interpretable and clearly labelled?
- Should we adjust colours, grouping, or significance annotations?
- Confirm alignment with proposed analyses in the proposal.

# Meeting Minutes - March 11


## Aim 1 Feedback:

**Add Faith’s PD for another alpha diversity analysis**

Make sure that when we summarize all of our plots, we also summarize what the plots mean and represent 
- (eg: shannon diversity tells us richness and abundance)

Shannon diversity by gender within intestinal subtype shows no significant difference, we can bring up significant test limitations within our manuscript.

Weighted UniFrac PcoA by gender for both the diffuse and intestinal subtype results, we can conclude:
- Shannon diversity is significant
- Beta diversity is not significant

- Go from a male angle when describing the data “decrease in Shannon diversity and increase in beta diversity ???”


## Aim 2:
For the Top 15 Genera by Subtype plot, we should use a customized colour palette so that we can easily differentiate the genus in our plot cause right now the colours blend together because they are rainbow and transparent
- We could use colour-blind friendly colours
- Don’t worry about this taxa bar plot

We probably won't include this figure in our final manuscript anyway (no significant power)
- If we see the same trend in the other plots, the other plots are likely more valuable, but we can leave this in the supplemental figures

For the Differentially Abundant Genera Across Gastic Cancer Subtypes plot
- Shows relatively abundance, visualized by a heatmap
- “Okay, thats pretty cool”
- Cluster on the right side has a higher abundance than Heliobacteria
  - WE should look into those genera of bacteria 
- Significant value is from the indicated species analysis
  - The right side shows ancom-bc plot
- The colour density is not reflecting the log change, it might be showing the cube value
  - Since they're all significant, we should use colour to represent the log-fold change
  - This could highlight other genera of interest that are abundant
- In that case, we do see the same as the taxonomic bar plot, so we can use this figure instead of the last one.

Between the heatmap and volcano plot:
- Volcano plot is cleaner
- Separates out the genera so it's easier to interpret
- Do we want to do a cutoff for the log2 fold change
  - Bigger than 1 or smaller than -1?
    - **We should go with a filter cutoff of anything less than -1 to highlight and make the higher genera pop out more**

Do we want to look at ISA?
- Different abundance between the subtypes, we should consider other metrics that consider more than abundance (eg: subtype-specific indicator taxa)
  - We already have this
- The differential abundance results are good enough, so it won't add more to this, so just these are good

For ISA:
- We didn't look at gender specific differences, which we should do because it's a pairwise analysis, so we should look at gender and subtypes.
- Our stats are really good; the stat value cutoff is 0.7, so ours are very good, very significant “good work.”
- If we filter, the cutoff is 0.7 for the stat value and 0.05 for the q value.

Diffuse is generally associated with a more severe disease progression, since we don’t have an actual severity index value, then our previous volcano plot is good!


Saman proposed focusing on just Heliobacteria:
- Hypothesis: If we go back to aim 1 and just focus on heliobacteria, since the abundance increases for severe disease progressions, then 
- The previous group didn’t find any heliobacteria, but since the original paper did and we did, then we can look at the diffuse subtype only and look at the abundance difference between males and females. 
  - We don't have to consider the previous group since the Heliobacteria presence is annotated.
- The diffuse subtype samples should relatively be positive for heliobacteria presence
- Saman proposed to control for the heliobacteria-positive and look at differences within the diffuse subtype.
  - Our sex-based route is more novel and better
- **TLDR: look at relative abundance differences of Heliobacteria within the diffuse subtype between male and female**

## Three Main Takeaways:
- Add Faith’s PD as another alpha diversity analysis in our project.
- For the heatmap, we should do a filter and a cutoff of anything less than -1.
- We should perform another analysis where we look at relative abundance differences of Heliobacteria within the diffuse subtype, comparing the male and female groups
