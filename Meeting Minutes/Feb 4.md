# Agenda - February 4
- Present our chosen dataset 
- Propose potential research questions and areas of exploration

# Meeting Minutes - February 4

### Proposed Datasets and Feedback (Bolded)
<ins>Alcohol dataset:</ins>
- How smoking status affects alcohol consumption and the microbiome
**ALREADY DONE LAST TERM**

<ins>Alcohol affecting microbiome in males vs females</ins>
- 80 men vs 200 females
**Sample sizes are okay!**
- Alcohol consumption is annotated as a continuous variable (grams)
  - We would do different levels / amount of alcohol
**We could do male vs female for smoking and non-smoking**
  - More of a yes or no

**WE COULD COMBINE WITH THE VAPING DATASET!**

- **Small sample size = anything less than 3 (anything above this allows us to do statistical significance)**


<ins>Cardia vs Non-cardia (gastric cancer)</ins>
- Different parts of the stomach
**ALREADY DONE LAST TERM**

<ins>Infant dataset</ins>
- Different modes of delivery
**ALREADY BEEN DONE**
- Use of antibiotics
**SOMEONE TRIED, VERY FEW SAMPLES**
  - Interested in seeing if babies develop antibiotic-resistant genes
- Diet
**ALREADY BEEN DONE**
  - Weight change over the first week and microbiome
**COULD BE INTERESTING**
**CONTROL FOR BREAST-FED VS FORMULA**

<ins>HIV (therapy impact on microbiome)</ins>
- Do specific drugs have an impact on the microbiome
**ALREADY BEEN DONE**


<ins>OUR THREE SELECTED ONES:</ins>
- Alcohol difference in microbiome (male vs female in smoking vs non-smoking)
- Use of antibiotics in babies
- How differences in weight impact the microbiome

### General Feedback
- FUNCTIONAL ANALYSIS (for antibiotics and weight)
  - We have to learn ourselves
   - Everyone is doing it, so not that bad
  - Worried if we get enough statistical data analysis

- Binning might be an issue for the weight dataset
- They are Not sure of how binning would work
  - We could bin them by breast-fed vs formula
    - We DONT want to focus on feeding habits (compounding variable)
- We could graph time vs weight
  - Two lines (breast-fed vs formula)
- For alpha diversity (shannon)
  - We expect different profiles
- We don’t want to focus primarily on breast-fed vs formula (ALREADY BEEN DONE)
  - It’s tricky to bin, we could calculate how much each baby grew, and bin them by % change in weight.
- We can focus on the pre-existing bin of male vs female
- OR bin them on underweight, normal weight, overweight

-A lot of things to focus on (male vs female, breast-fed vs formula)
  - We can go forward with both in parallel 

- SUGGESTION: Look at first few weeks, then look at % change in weight over 12 months (lots of math)

### OUR FINAL SELECTED DATASET/RESEARCH QUESTION
<ins>EXPERIMENTAL AIMS/DESIGN</ins>
- Saman brought up the gastric vs intestinal cancer:
- They collected the lesions from the body and the antrum of the stomach
- We have four comparison groups (that are already defined)
  - Mixed
  - Diffused
  - Intestinal
  - Healthy

Diversity metrics (how many shared microbes vs distinct)
- Measure within each site
  - Then type
    - Look for effective sites (?)
- Core microbiome
- Indicator taxa (which species are highly associated with each group)
- DESeq (between shared, which are rich, depleted between certain conditions)
- Functional analysis  (which things are up vs downregulated)
  - Picrust2

**H. pylori defines where the cancer lesion is**

Keep in mind:
- Dr. Evelyn will send us the server credentials that we will use to run our work on
  - Its faster than working in our own individual containers
- We can start working on data processing 
  - We need this done before our proposal
  - If we hit a point about filtering/trimming
    - Email Bessie with the graph to make sure it's good
  - Make a habit of making a screen

<ins>FOR NEXT MEETING:</ins>
- Read over the proposal before the next meeting
- Find literature research on the different types
  - Literature on why we think the microbiome will drive different cancer subtypes?
  - Better understanding of tissue sites and whatnot for our project
- Take a good look at what we are controlling for (male vs female, etc.)

