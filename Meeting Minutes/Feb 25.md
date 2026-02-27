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


# Meeting Minutes Feb 25th, 2026

### Agenda for the day:
-	Go over research question
-	Go through data processing
-	QIIME2
-	Rarefaction depth
o	No chance to go over depth with Bessie
o	Depth chosen only lost 2 samples, not maximum for features but very close so therefore chose 10,222

### By next week try to do alpha/beta diversity metrics
-	What stats are we going to perform?
This week expect feedback
-	Other TA will grade and give feedback on project proposal

No data to look at this week
-	If clear on data to analyze can move forward

Ideally one phyloseq object generated and 

Rarefaction phyloseq only needed for diversity analyses
Not needed for others so use non-rareified phyloseq for those

Aims will be split, with each person tackling one
-	Aim 4 will be split into 2

For samples we are left with 15 intestinal and 30 diffuse
-	After rarefaction we lost 2 samples and we’re left with 43
-	Both samples lost were intestinal

### Going over proposal aims:
For Aim 1:
-	Maybe justify why you use these metrics
-	Possible point of comment for TA, only address if it comes up

For team proposal after feedback we will receive TA comments
-	Address TA comments for full marks
-	Revised proposal to be resubmitted and marked by Evelyn

For Aim 2:
-	In general in our aims we don’t mention the analyses
-	Should try to relate what you’re going to do and why to what it tells you

For Aim 3:
-	For aim 3 add in diversity and composition

For Aim 4:
-	This is fine, no real feedback

### In proposed approach
For Aim 1:
-	In stats tests no linear model is possible, not a continuous variable
For Aim 2:
-	For ISA beyond P value may need some adjusted stats value
-	What does your bar plot represent? Is it for ISA?
-	Normally ISA just outputs stats and table and adjusted P-value
-	Unless you want to go back and calculate relative abundance of ASVs
-	Again for heatmap you need to use a measurement e.g. abundance
For Aim 3:
-	Maybe merge with aim 1? Like Aim 1A/AB?
-	Since no large difference
-	Will basically be very similar to Aim 2, just within each subtype male and female
For Aim 4:
-	Looks good

### Future meetings will have more data so we can discuss troubleshooting
