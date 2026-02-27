# Meeting Minutes - January 28
### Introduction to Weekly Team Meetings 
- We will lead discussion using your agenda
- Take notes, discussion points, action items etc to put on github 
- PART OF COURSE EVAL

### Before every meeting 
- Write down anything from meeting
- Questions
- Concerns abt script
- Can ask about analyzing results 
Complete tasks by night before meeting!!
- Upload finished analysis 
- Can always email for help outside meetings 

### Github 
- consult Module 18 
- One public github repo
- MICB475_Team11
- Organize by type of analysis 
Each folder - have subdirectories 
Figures, codes, etc

### README File  
- Hve bullet point discussed during meeting
- And what want to talk abt after meeting 
- PART OF COURSE EVAL 
- Organized and help writing 
- HAVE THIS SET UP BEFORE NEXT MEETING (Agenda and what to do)
- @ bessie and evelyn 

### Course timeline
- Pick one or two datasets next week- At end of next week - commit to 1 or 2 dataset 
- And preliminary questions want to ask 
Week 4 - brainstorm potential topics 
Week 6 - proposal due

## Pitches 
### Gastric cancer: Lauren histological subtypes and the microbiome

Gap / novelty
- Wang et al. grouped all gastric cancers together and did not stratify by Lauren classification
- Metadata includes intestinal, diffuse, and mixed subtypes
- Very few studies compare gastric microbiomes across these subtypes
- Wang et al. observed lower H. pylori in intestinal-type tumors vs diffuse-type but did not analyze overall community differences
  
Key background (Lauren classification)
- Intestinal type: gland-forming, inflammation-driven, strongly linked to H. pylori, more localized
- Diffuse type: poorly organized cells, infiltrative growth, weaker H. pylori link, more genetically driven (e.g. CDH1)

Core distinction:
– Intestinal = structured, chronic inflammation pathway
– Diffuse = disorganized, infiltrative pathway

Research question
- Do gastric microbiome diversity and taxonomic composition differ between intestinal and diffuse gastric cancer subtypes?

** Focus on Lauren's classification 


### Gastric cancer: tumor location (cardia vs non-cardia) and the microbiome

Gap / novelty
- Wang et al. did not analyze microbiome differences by tumor location
-Metadata labels tumors as cardia vs non-cardia
- Cardia and non-cardia cancers are epidemiologically distinct but rarely compared microbiome-wise
- Limited prior work suggests higher Helicobacter abundance in cardia cancers, but evidence is sparse
  
Key background (tumor location)
- Cardia: upper stomach near esophagus, linked to reflux and obesity, weaker H. pylori association, greater oral/esophageal influence
- Non-cardia: body/antrum, strong H. pylori and chronic gastritis association, classic gastric cancer pathway
  
Research question
- Does tumor location (cardia vs non-cardia) associate with distinct gastric mucosal microbiome community structure and taxa?

** Focus on tumor location

### Multiple sclerosis: NSAID use as a microbiome confounder

Gap / novelty
- Many MS patients use NSAIDs for pain
- NSAIDs are known to alter gut microbiota
- The original study did not control for common non-MS medications
  
Rationale
- NSAIDs have been linked to shifts in Acidaminococcaceae, Bacteroides, and Enterobacteriaceae
- Could explain part of the MS-associated dysbiosis
  
Research question
- Is NSAID (Non-steroidal anti-inflammatory drug) use associated with distinct gut microbiome composition in MS patients, and does it contribute to observed differences vs healthy controls?


### Alcohol Consumption

Background 
- Drinking alcohol is associated with variation in human oral microbiome (Fan et al., 2018)
- Smoking alters oral microbiome diversity and composition
    - lower alpha diversity (Yu et al., 2017)

 Gap 
 - Thomas et al. (2014) investigated how alcohol and tobacco consumption affects bacterial richness in oral cavity mucosa biofilms
 - small sample size of 22 subjects
 - investigated 7 chronic and heavy users of both substances and 6 active smokers (rest were controls)
 - found a signifcant decrease in species richness in only smokers and smokers and drinkers compared to control
 - could not find a paper investigating drinkers and how smoking could have an affect
 - few studies investigated whether smoking modifies or changes the strength of associations between alcohol intake and oral microbiome 

Research question
- Does smoking status modify the association between alcohol consumption and oral microbiome diversity and composition? (ie is the effect of alcohol different in smokers vs non-smokers)

Research question
- Does alcohol consumption affect the oral microbiome differently in men versus women?
- 

### ART and microbiome

Possible Project Ideas (Andy)

10.1186/s40168-023-01718-4
HIV dataset

Does being on certain anti-retrovirals correlate with microbiome differences after controlling for BMI and antibiotics?
-	BMI
-	Antibiotics
-	Cohort_short
-	Current_art
Possible confounders: Arm_short, Age, Gender

### Gastric cancer and metabolic markers

https://doi.org/10.3389/fmicb.2020.00997
Gastric Cancer dataset

Do metabolic markers (like glucose, cholesterol triglycerides) explain microbiome variations independent of cancer stage?
-	Group
-	BMI
-	Blood Glucose
-	Total Cholesterol
-	Triglyceride
Possible confounders: Age, Gender, H. pylori test, Site


### Possible Neuro-Imbalances Within Different Disorders: 

Background / Novelty
- Most of the previous UJEMI papers focused on a single disease, but for this project, we will look at two different neurological conditions/disorders and compare the metabolic signatures found in them. 
- Instead of focusing on what makes up and categorizes each disorder, we will instead look at the shared 'neuro-inflammatory' metabolic signatures.
- We will use the "Depression" dataset as a baseline for mental health and see if those same markers appear in the Parkinsons dataset
- We will focus on certain "pro-inflammatory" metabolic signals that are common in both (eg: copper and potassium), which we hope will suggest that the microbiome doesn't just reflect a specific disease, but a general state of neurological stress.

Research Question
- Across patients diagnosed with either Parkinson’s Disease and Depression, are there conserved "neuro-inflammatory" metabolic signals or distinct differences across different neurological conditions?

### Infant Microbiome

This project uses an infant health and exposure metadata dataset to ask microbiology-focused questions about how early-life factors shape the gut microbiome. Although the dataset is primarily metadata (e.g., delivery mode, feeding type, antibiotic exposure, symptoms, and growth/appetite measures), each sample can be linked to microbiome sequencing outputs (e.g., 16S rRNA or shotgun profiles) using `#SampleID`. By combining metadata with microbial community data, we can test how specific exposures are associated with changes in microbial diversity, community composition, and the relative abundance of key taxa.

Key research directions include:
- **Early-life exposures -> microbiome:** Do delivery mode (C-section vs vaginal) and antibiotic exposure predict reduced microbial diversity or shifts in taxa such as *Bifidobacterium*?
- **Diet/feeding -> microbiome:** Does feeding type relate to the dominance of milk-adapted microbes and overall community structure?
- **Symptoms/illness -> microbiome signatures:** Are gastrointestinal symptoms (colic, diarrhea) associated with distinct microbiome profiles or dysbiosis-like patterns?
- **Allergy/immune outcomes -> microbiome:** Do infant or maternal allergy indicators correspond to differences in microbiome diversity or immunomodulatory taxa?
- **Growth/body composition -> microbes:** Is infant growth (e.g., weight-for-length z-scores) associated with microbiome “maturity,” diversity, or specific microbial patterns?
- **Family/shared environment:** Do related individuals (same household) have more similar microbiomes than unrelated individuals, suggesting microbial sharing?

Overall, the goal is to use real-world infant exposure and health variables to generate testable hypotheses about microbial ecology in early life and how microbiome development may relate to health outcomes.

### Infant Microbiome
Does early antibiotic exposure influence infant gut microbiome maturation and appetite-related behaviours, and do these effects differ between male and female infants?

