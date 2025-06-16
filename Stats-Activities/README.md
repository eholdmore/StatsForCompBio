# Activities

This folder contains the materials for the five activities included in the workshop.

### Installation Options

**Recommended**: Use `renv` for exact package versions:
```r
install.packages("renv")
renv::restore()

Alternative: If you have issues with `renv`, use the setup script:
```r
source("setup.R")

## Activity 1: Sample Size & Power Analysis

Power analysis calculations in R: How many samples do you need to detect a certain effect? What power do you have to detect signal with given a certain sample size?

## Activity 2: Probability Distributions

Interactive simulation from [Seeing Theory](https://seeing-theory.brown.edu/probability-distributions/index.html) demonstrating probability concepts.

1. After following the link above, click on the "Discrete and Continuous" box on the right.
2. Use the buttons and dropdown menu to explore the forms of different probability distributions. What types of data might be suitable for each?

ProbStats.org also offers an [interactive visualzation of probability distributions](https://probstats.org/) that you may find useful. 

## Activity 3: GLMs in Action

Simulate gene expression count data and fit both Poisson and negative binomial generalized linear models (GLMs) to it. You'll compare model fits, assess overdispersion, and interpret coefficients, with a focus on log fold change. This activity introduces concepts foundational to RNAseq differential expression tools like DESeq2.

## Activity 4: Correcting for Multiple Comparisons

Hands-on multiple comparisons analysis in R.
Adapted from [Adapted from Quantitative Understanding in Biology: 1.6 Multiple Hypothesis Testing and Non-parametric Tests](https://physiology.med.cornell.edu/people/banfelder/qbio/resources_2021/2021_1.6_multiple_hypotheses_and_non_parametric_tests.pdf) by Jason Banfelder

## Activity 5: Principles of Data Visualization

Design and creation of data visualizations using real genomic datasets in R, followed by group discussion on best practices and interpretation.

## Activity 6: Ethical Dilemmas 

Case studies and group discussions exploring ethical dilemmas and practical challenges faced in computational biology research.

The case study below is a summary based on Case Study 1 from the [Princeton Dialogues on AI and Ethics Case Studies](https://aiethics.princeton.edu/case-studies/).

### Summary

Type 2 diabetes is a chronic condition causing elevated blood glucose levels, leading to serious medical issues. The disease affects 30.3 million people in the US, with low-income individuals, Native Americans, Black people, and Hispanic people being most affected. While type 2 diabetes is manageable with proper care, treatment can be expensive and onerous. Individuals with diabetes must test their blood several times per day and self-administer insulin injections based on complex calculations. Many diabetics eventually become delinquent in their testing and insulin shots. Because the calculations are complicated and imprecise, even compliant patients often miscalculate their optimal insulin dosages. To address this, researchers from St. Marcarius-Alexandria University Hospital developed a multiplatform application called Charlie, which uses artificial intelligence technologies to make diabetic care easier, more holistic, and more accessible. Charlie uses smartwatches' biosensors to test blood glucose through the skin, calculating the optimal level and type of insulin for each user. It also includes a data collection platform, providing accurate, individualized insulin dosage recommendations and personalized reminders. Charlie also contains a forum for information sharing and social networking, aiming to counteract misinformation about diabetes and healthcare. The app uses natural language processing techniques to analyze emerging discourse and improve customized treatments.

Charlie passed the university IRB process with minimal risk to individual subjects. The technology uses existing frameworks to provide more efficient, accurate, and personalized care. After IRB approval, Charlie was rolled out in a clinical trial at the University Hospital, showing positive results for users. However, racial minorities did not experience the same positive results as white users. The social networking forum was also mixed, with conflicting reports and hostile arguments. The research team formed a "Benevolence Team" to address issues of inequality and speech on the forum and explore other ways to improve Charlie. After polling users and hosting open discussion threads, the researchers decided on several initiatives.

First, they wanted to know why some users were more compliant with the app than others, and what could be done about it. Using its natural language processing and machine learning capabilities, they discovered a correlation between those who read articles about disagreements in the scientific community, the theory being that this made them less trusting of medical advice in general. Charlie’s research team concluded that they could improve compliance by inundating high risk users with scientific information that tied type 2 diabetes to lifestyle choice and minimizing their exposure to disagreements within the scientific community.

Second, to combat hostile discussion and echo chambers, the team will revise the app's content-based recommendations (automated content moderation algorithms or ACMAs) on its social media forum to prioritize information that is most generally accepted or approved by users. They will also introduce individualized content filtering based on the app’s extensive data collection and analysis to create a more pleasurable, personally relevant experience. 

Third, to determine the most effective approaches to type 2 diabetes management, the team will perform controlled mini-experiments on a subset of users using a "multi-armed bandits" (MAB) (aka “explore vs exploit”) model. This approach allows researchers to test different treatment options and collect data about which interventions work best under different circumstances. Rather than provide each user with the treatment protocols that appear most appropriate for them given the limited data available at the outset, researchers purposely provide some users with sub-optimal solutions to “explore” outcomes and gather additional data to enrich Charlie’s algorithms. While acknowledging that this approach may not have been ideal for individuals in each instance, they justified the experiments by explaining that they would increase understanding about healthcare protocols and produces better outcomes in the long-run.

A new version of Charlie, released to previous research subjects, showed improvements in healthcare metrics, particularly among minority and at-risk users. The revisions to Charlie's ACMAs, which favored general advice and provided individualized content, alleviated the "cacophony of voices" on the social networking forum. However, when the research team published their results, concerns about censorship, special treatment for at-risk users, and the MAB approach arose, leading to debates.


American Statistical Association's [Ethical Guideline for Statistical Practice](https://www.amstat.org/docs/default-source/amstat-documents/ethicalguidelines.pdf?Status=Master&sfvrsn=bdeeafdd_6/) 
