# Agenda - March 18
### Aim 3 Figures Review:
- Heatmap: Review the functional signatures, and talk about the diversity of the Diffused subtype, the heterogeneity found in the Diffused subtype (there's an active cluster, an inactive cluster, and a mixed cluster)
- Volcano Plot Validation: Talk about why the Diffused type appears significantly more metabolically active than the Intestinal type, and if q < 0.05 is good.
- Talk about which pathways are upregulated, and if there are any recognizable patterns.
- Do we want to include sex? (I remembered this last minute, and we didn't talk about it in our proposal, but it could be a good idea to note the different sexes on the heatmap and volcano plot)
  - We could use different symbols to denote the different sexes in the volcano plot, but how can we differentiate in the heatmap?

### Reviewing Any Updated Figures From Last Week
- Volcano plot: show stronger difference by adding a +/- 1 log2fc
- Box plot + Volcano plot: helicobacter relative abundance in diffuse gastric cancer
- Heatmap: change colors to represent log2fc?

### Next Steps?
- Prepare to write our manuscript draft
- Prepare for presentation (?)
- Anything else?

# Meeting Minutes - March 18

## Aim 3 
- Functional analysis btw diffuse and intestinal subtypes 
- Which pathways are upregulated 

### Heatmap 
- Q: Did we use code provided from the module? No 
- So in the code we would need to annotate all the terms but its fine if we can extract names later
Diffuse group 
Observations: 3 different subgroups
- One subgroup to left very active
- Most pathways upregulated
- Near center - not very upregulated 

### In the proposal, we didn’t differentiate btw male and female but should we look at sex as a factor?
Males might have more upreg??
- On the very left side of heatmap there are 2 pathways enriched
- Could sex play a role in this? 
To incorporate sex: can add different packages 
- Don’t actually need to change anything from heat map / create 2 different heatmaps (one for male and one for female)
- Just send code to bessie and she can add legend for us 
### HOWEVER, from this heatmap, there are no clear clustering / separation
- She doubts that the current clustering is due to sex 
- Diffuse and intestinal seem kinda mixed together 

### Volcano plot 
What data was log transformed? The q values?
- Didn’t do log2fold 
- Apparently already have logfold from the pie crust script 
- One of the data tables had log2fold = from ancomb? 
### DOUBLE CHECK THE DIRECTION!
- Based on heatmap, we think its intestinal (right) vs diffuse (left)
- Bc we have more pathways upregulated in diffuse 
- Currently colored by q values, cant tell what is up/ down regulated 
### The directionally matters 
- If flip the group, it will change the values 
- If do differential abundance of pathways: Diffuse vs intestinal OR Intestinal vs diffuse

#### Results of Aim 3 show
- Diffuse subtypes - more upregulated and more metabolically active
- Makes sense because diffuse have stronger disease progression so being more upregulated makes sense 

### Next Step 
Look into pathway, does it agree with our hypothesis?
- What are the pathways?
- Current results there are no descriptions, need to search them up?
- Google these database, map these pathway names
- Search “How map names to the descriptions”

Shows the relative abundance 
- Back in Aim 1, when we compare shannon diversity within diffuse
- Male and female different significantly
We found high abundance in males = richness decrease

Therefore, if it makes more sense/ will be more interesting to look at the functional composition in diffuse between male and female

## Aim 2 
- Focus on general diffuse and intestinal

### Box plot
- Found enrichment in heliobacter and it differs by sex
-   Can do a functional analysis for this comparative group in Aim 3 

### Volcano plot 
Lactobacillae - produces acid, lowers immune system 
- Not found in females
- Supports previous finidngs 

## Next Week
What information do you need to show what we’ve seen so far
- organize figures by making slides
- show what we found:

### Start off in Aim 1 Diversity Analysis 
We compare the 2 subtypes of gastric cancer and looked at the overall diveristy differences btw sex
Diffuse condition: shannon is not significant
- don't need to do faith PD analysis btw 

### Move into Aim 2 to look at the compositional difference
- Diffuse condition is more enriched with helioabacter = makes sense that overall abundance will decrease
- Our boxplot shows that there is a sex difference where heliobacter is more enriched in males (this corresponds to the depleted alpha diversity in males observed in Aim 1) 

### Therefore, it would makes sense to find differences in functional  
- Then go overall functional difference
- Then focus on the functional composition specifically in DIFFUSE between MALE and FEMALES 
