**Significance:** Gastric cancer is one of the most common cancers and the 5th in terms of mortality
(https://gco.iarc.who.int/media/globocan/factsheets/cancers/7-stomach-fact-sheet.pdf)

**Lauren**
Most commonly, they are gastric adenocarcinomas, which are graded based on the Lauren Histological Scores.

Lauren Histology essentially divides them into 2 major subtypes, intestinal and diffuse-type.

These have huge implications on actual clinical outcomes, with diffuse being associated with peritoneal spread and poorer outcomes. 
https://pmc.ncbi.nlm.nih.gov/articles/PMC4840723/#b1-ol-0-0-4337

Lauren histology score also has implicatons for how it is treated, and possible targetted chemotherapies.

"Diffuse adenocarcinomas are defined by their growth pattern as tumours infiltrating the stroma as discohesive tumour cells arranged singly and in small clusters. The intestinal type is defined by its cytoarchitecture and characterized by cohesive cells which form gland-like structures. Mixed tumours have both an intestinal and diffuse component, while indeterminate types include most of the undifferentiated tumors"
https://link.springer.com/article/10.1186/s12957-017-1187-3

**Feasability**

Wang et al. (paper we're basing dataset on) profiled gastric mucosal microbiota across histologic stages and reported diversity loss across progression, and categorized Laurén subtype differences

Sarhadi et al. already showed diffuse adenocarcinoma had the lowest alpha diversity, with Enterobacteriaceae increased in both intestinal and diffuse subtypes, and Bifidobacteriaceae decreased specifically in diffuse disease.
https://pmc.ncbi.nlm.nih.gov/articles/PMC7888145/

It's also very possible that microenvironmental differences between intestinal and diffuse tumours would create different ecological niches for microbial colonization.

Account for confounders like H. pylori and age, and there should be a thread to pull on.

(Chronic H. Pylori is a huge risk factor and is carcinogenic to humans)
https://www.cancer.gov/about-cancer/causes-prevention/risk/infectious-agents/h-pylori-fact-sheet

**Why would microbiome and cancer subtypes be associated?**

Diffuse adenocarcinomas undifferntiated cell types, and aggressiveness may be linked to lower microbiota diversity. It also has much poorer immune cell invasion.
(https://pmc.ncbi.nlm.nih.gov/articles/PMC7888145/)

Intestinal tumors are arranged in tubular/glandular structures, and are often associated with intestinal metaplasia (esophagus/stomach tissue remodeled to look like intestine, often a precancerous condition).
Associated with lymph/vascular invasion and scattered lesions. Also enriched in elderly men, and affects the gastric antrum. 
https://pmc.ncbi.nlm.nih.gov/articles/PMC4840723/#b1-ol-0-0-4337

Because of these patterns, differences in how they arise, and the local tumour microenvironment's impact on the microbial niches, it would be expected that some association could be found. 

**Further Analysis**
Another study by Yang Et al. confirms the use of the Lauren classification to describe the differences in the microbial composition of gastric cancer. However, they additionally used a classification known as the ZJU index to describe metabolic differences between different gastric cancer subtypes. The ZJU index is a NAFLD (Nonalcoholic fatty liver disease) biomarker that uses fasting plasma glucose (FPG), TG, ALT/AST ratios and BMI to determine metabolic identity (this was developed in China). Our dataset contains all of these values except for the ALT/AST ratios, and we don't have fasting plasma glucose, only blood glucose levels, so we won't be able to apply this classification, but its good to know that it exists.
(https://pmc.ncbi.nlm.nih.gov/articles/PMC12598303/)
