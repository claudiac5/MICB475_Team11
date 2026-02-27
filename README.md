# MICB475_Team11

## Agenda and Meeting Minutes 
### January 
[Jan 28](./Meeting%20Minutes/Jan%2028.md) 

### February 
[Feb 4](./Meeting%20Minutes/Feb%204.md) | [Feb 11](./Meeting%20Minutes/Feb%2011.md) | [Feb 25](./Meeting%20Minutes/Feb%2025.md) 

### March 
[Mar 4](./Meeting%20Minutes/Mar%204.md) | [Mar 11](./Meeting%20Minutes/Feb%2011.md) | [Mar 18](./Meeting%20Minutes/Feb%2018.md) | [Mar 25](./Meeting%20Minutes/Feb%2025.md) 

## Lab Notebook 
### Data Processing 
- [P01](./Notebook/P01.Rmd) - QIIME 2 Demultiplexing
- [P02](./Notebook/P02.Rmd) - QIIME 2 Denoising and Clustering
- [P03](./Notebook/P03.Rmd) - QIIME 2 Taxonomic Analysis
- [P04](./Notebook/P04.Rmd) - QIIME 2 Alpha-rarefaction
- [P05](./Notebook/P05.Rmd) - QIIME 2 Diversity Analysis

### Data Analysis  
- [P06](./Notebook/P06.Rmd) - Aim 1: Investigating the Microbiome Composition in Intestinal vs Diffuse Gastric Cancer and Evaluating the Influence of Biological Sex
- [P07](./Notebook/P07.Rmd) - Aim 2: Identify the Specific Bacteria that Best Distinguish Intestinal vs Diffuse Cancers
- [P08](./Notebook/P08.Rmd) - Aim 3: Functional Analysis and Profiling of Microbial Communities Within Subtypes

## Final Script 
### QIIME 2
- Data Processing  

### R
- Aim 1:
- Aim 2:
- Aim 3: 



# Agenda - February 11
-  Discuss literature research on the different types/  why we think the microbiome will drive different cancer subtypes?
- Discuss our controls (male vs female, etc.)
- Discuss next steps for proposal

## Research question: How do gastric microbiome diversity and overall taxonomic composition differ between intestinal and diffuse gastric cancer subtypes?
- Healthy vs Cancer: Compare each subtype against healthy controls  to identify general tumor-associated shifts.
- Compare subtypes pairwise To test hypotheses of subtype-specific microbiomes.

## Meeting Notes — February 11

**Note:** Missed the first ~5 minutes of the meeting.

### Research Question (RQ)
How do **gastric microbiome diversity** and **overall taxonomic composition** differ between **intestinal** and **diffuse** gastric cancer subtypes?  
- Considering adding **sex** as a key factor.

### Working Hypothesis
We hypothesize that there will be **compositional and functional differences** between gastric cancer subtypes, and that these differences will be **affected by sex**.

### Proposed Analyses (Initial Plan)
- Run **PCoA** analyses **separately for males and females**.
- Focus on **beta diversity** comparisons.
- There are **3 sampling sites** available:
  - If one site appears very different, consider **merging the other two** sites (site handling to be decided after exploring clustering/differences).

### Background / Rationale Discussed
- Need to summarize **what research has been done** on microbiome differences between diffuse vs intestinal gastric cancers.
- Key points raised:
  - **Diffuse carcinomas** tend to be **less differentiated** and may have **lower microbial diversity**.
  - Another paper suggests **intestinal tumors** are associated with **different microbial groups**.
- Project novelty discussed:
  - Including **3 subtypes**: **diffuse**, **intestinal**, and **mixed**.
  - Team is leaning toward emphasizing **sex** and how it impacts the microbiome in relation to cancer.
  - Potential mechanism angle: **sex hormones (e.g., estrogen)** and links to gastric cancer.

### Proposal / Rubric Notes
- Avoid being **conclusive at the beginning** of the proposal.
- Provide clear **logic** for expected outcomes (explain *why* we think results may look a certain way).
- Explicitly describe **differences between diffuse**


# Agenda - February 25

## Restate research question:

- How do gastric microbiome diversity and taxonomic composition differ between intestinal and diffuse gastric cancer subtypes, and does biological sex modify these subtype-associated microbiome patterns?

Confirm final sample subset:

- GC only

- Intestinal vs Diffused type

- N = 45 (15 intestinal, 30 diffuse)

## QIIME2 Preprocessing

- Demultiplex summary + quality inspection

- DADA2 denoise-single (no truncation)

Generated:

- Feature table (ASVs)
- Denoising stats
- Assigned taxonomy using SILVA
- Removed mitochondria/chloroplast
-  Built a rooted phylogenetic tree

Filtered to:

- Group = Gastric cancer (GC)

- Histopathology = Intestinal type or Diffused type


## Dataset Summary
- ASVs after denoising: ~20,750
- Final subset size: 45 samples
- Rarefaction depth chosen: 10,222
- 43/45 samples retained after rarefaction


### 1. Forward Read Quality Profile

The forward-read quality plot was inspected before denoising.  
Quality scores remained high across the read length, so no truncation was applied during DADA2 denoising.
![Demultiplexing summary](images/Demux.png)

### 2. Alpha Rarefaction Curve

Alpha rarefaction curves were generated on the filtered GC intestinal vs diffuse subset.  
A sampling depth of **10,222 reads** was selected based on curve stabilization and sample retention (43/45 samples retained).

![Alpha rarefaction curve](images/Rarefaction.png)

## next steps
- Aim 1 – Community Diversity (Alpha diversity and beta diversity) 
- Aim 2 – Differential Taxa (Agglomerate to genus level) 
- Aim 3 – Sex Interaction (Alpha diversity and beta diversity: subtype * sex)
- Aim 4 – Functional Prediction (PICRUSt2 and Map to MetaCyc pathways)





