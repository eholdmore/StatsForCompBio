### Activities for Statistics for Computational Biology Projects Workshop ###
### Erica M. Holdmore ###
### Last updated: June 16, 2025 ###

#### Setup Environment ####
install.packages("renv")
renv::restore()

#### Activity 1: Power Analysis ####

# load package
#install.packages("pwr")
library(pwr)
# calculate effect size 
pwr.anova.test(k=2,            # compare 2 groups (treatment and control)
               n=20,          # sample size is 20 (per group)
               sig.level=.05,  # significant level = .05
               power=.8)       # power = .8

# Exercise: determine how many participants you would need in each group (sample size) 
# to have a power of 80% and a moderate effect size of 25% for each of the following tests.

# one-way ANOVA
# See example above.
?pwr.anova.test()
#pwr.anova.test(k=,f=,sig.level=,power=)

# GLM
# Note: instead of groups (k), this command requires degrees of freedom for the
# numerator (in our case u = 2-1 = 1); v will be the min. sample size.
?pwr.f2.test()
#pwr.f2.test(u=, f2=, sig.level=, power=)

# paired t-test (two tailed)
# Note: for this test, you will need to enter some additional info about
# the type of test (type="paired") and alt. hypothesis (alternative="two.sided")
?pwr.t.test()
#pwr.t.test(d=, sig.level=, power=, type=, alternative=)

# independent t-test (one tailed - "greater")
# Note: this one is more tricky - since it is an unbalanced test, it requires
# sample size for each group (n1 and n2). 
# Assuming you have 400 participants total, play around with different n1 and n2 
# values to see how you can achieve a power of 80%.
# What happens when you decrease the total sample size (n1+n2)?
?pwr.t2n.test()
#pwr.t2n.test(d=, n1=, n2=, sig.level=, alternative=)

# chi-squared test
# Note: for this test, will need to include the numerator degrees of freedom
# as you did with the GLM power analysis (see above)
?pwr.chisq.test()
#pwr.chisq.test(w=, p=, df=1, sig.level=)

#### Activity 3: GLMs in Action ####
# load packages
#install.packages(MASS)
library(MASS)

# Simulate simple biological count data (e.g., gene expression counts)
set.seed(123)
n <- 30  # samples per group
group <- factor(rep(c("control", "treatment"), each = n))
mu <- ifelse(group == "control", 10, 20)  # mean expression level
counts <- rpois(2 * n, lambda = mu)

# Quick look at the data
table(group)
head(counts)

# Combine into a data frame
dat <- data.frame(group = group, counts = counts)

# Plot raw counts
boxplot(counts ~ group, data = dat, main = "Simulated Count Data", ylab = "Gene Expression (counts)")

# Fit a Poisson GLM
pois_mod <- glm(counts ~ group, family = "poisson", data = dat)
summary(pois_mod)

# Check for overdispersion (a key concept in DESeq2)
dispersion <- sum(residuals(pois_mod, type = "pearson")^2) / df.residual(pois_mod)
dispersion  # >1 suggests overdispersion

# Fit a Negative Binomial GLM
nb_mod <- glm.nb(counts ~ group, data = dat)
summary(nb_mod)

# Compare models
AIC(pois_mod, nb_mod)

# Extract estimated effect size (log fold change)
coef(nb_mod)

# Interpret:
# What does the coefficient for 'grouptreatment' mean?
# Tip: Use exp(coef(nb_mod)["grouptreatment"]) to get the fold-change from the log scale.
# "The treatment group has ____ times the expected expression compared to the control group."

# Visualize model fits
library(ggplot2)
dat$fit_pois <- fitted(pois_mod)
dat$fit_nb <- fitted(nb_mod)

ggplot(dat, aes(x = group, y = counts)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", color = "blue", size = 3) +
  geom_point(aes(y = fit_pois), color = "orange", shape = 1) +
  geom_point(aes(y = fit_nb), color = "green", shape = 2) +
  labs(title = "Observed vs Fitted Counts", y = "Counts") +
  theme_minimal()

# Discussion prompts:
# 1. Why does Poisson underestimate variability?
# 2. What role does the link function play in GLMs?
# 3. How does this relate to modeling gene expression in DESeq2?

# Extra challenge:
# Modify the simulation to add a batch effect or covariate.
# dat$batch <- rep(c("A", "B"), times = n)
# Try re-fitting: glm.nb(counts ~ group + batch, data = dat)

#### Activity 4: Multiple Comparisons ####

# Adapted from Quantitative Understanding in Biology: 
# 1.6 Multiple Hypothesis Testing and Non-parametric Tests, by Jason Banfelder
# We'll simulate some data where we have 3000 genes
# and assume we have 6 BRCA1 and 6 BRCA2 replicates.

sample.size <- 6
effect.size <- 3.0

# A quick recap of power analysis:
# Let's also assume that the effect size is quite large (d = 3)
pwr.t.test(n=sample.size, d=effect.size, sig.level=0.05, type="paired", alternative="two.sided")
# A power of 0.999 means that we have a >99% chance of detecting a true positive.

# To evaluate how well each method works, we need to make some assumptions about what is true.
# Assume that 200 of the 3000 genes have significantly different expression.
# The code below sets up the data.
true.count <- 200
gene.count <- 3000
x <- list()
y <- list()
for (idx in 1:gene.count) {
  x[[idx]] <- rnorm(sample.size, mean = 0)
  y[[idx]] <- rnorm(sample.size, mean = ifelse(idx <= true.count, effect.size, 0.0))
}

# Then, we can calculate the raw (unadjusted) p-values...
p.raw <- numeric(length = gene.count)
for (idx in 1:gene.count) {
  p.raw[idx] <- t.test(x[[idx]], y[[idx]])$p.value
}
# and plot them using a histogram, where the left bar is all the significant values.
hist(p.raw,
     breaks = seq(from = 0, to = 1, by = 0.05),
     ylim = c(0, 400))
# Notice there are over 300 significant outcomes, which we expect.
# 0.05 * (3000-200) = 140 false positives in addition to 200 true positives

# Now we can try different correction methods and see how many/which outcomes are significant.
# Method 1: Bonferroni - controls FWER; very stringent
sum(p.adjust(p.raw, method = 'bonferroni') < 0.05)
which(p.adjust(p.raw, method = 'bonferroni') < 0.05)
# How many significant outcomes are observed using a Bonferroni correction?
# Are all the significant outcomes those that we know are true positives (first 200)?
# Method 2: Benjamini-Hochberg - controls FDR; less conservative
sum(p.adjust(p.raw, method = 'BH') < 0.05)
which(p.adjust(p.raw, method = 'BH') < 0.05)
# How many significant outcomes are observed using a Benjamini-Hochberg correction?
# Are all the significant outcomes those that we know are true positives (first 200)?

#### Activity 5: Data Visualization ####

# Adapted from Genomic Data Visualization and Interpretation. https://genviz.org/
#install.packages("ggplot2")
library(ggplot2)

## Choose your own adventure! Pick ONE of the example plots below.

#### Gene Coverage Plots ####
# https://genviz.org/module-03-genvisr/0003/04/01/gencov_GenVisR/
# load the dataset
covData <- read.delim("http://genomedata.org/gen-viz-workshop/GenVisR/STAT1_mm9_coverage.tsv")

# We need to reformat the data a bit.
# rename the columns
colnames(covData) <- c("chromosome", "start", "end", "TAC245", "TAC265")

# create a function to split the dataframe into lists of dataframes 
# one for each sample ("TAC245" and "TAC265")
samples <- c("TAC245", "TAC265")
a <- function(x, y){
  col_names <- c("chromosome", "end", x)
  y <- y[,col_names]
  colnames(y) <- c("chromosome", "end", "cov")
  return(y)
}
covData <- lapply(samples, a, covData)

#changes the names to match the samples
names(covData) <- samples

# install and load the BSgenome package and list available genomes
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BSgenome")
library("BSgenome")
available.genomes()

# install and load the mm9 BSgenome from UCSC
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm9")
library("BSgenome.Mmusculus.UCSC.mm9")
genomeObject <- BSgenome.Mmusculus.UCSC.mm9

# Install and load a TxDb object for the mm9 genome
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
BiocManager::install("TxDb.Mmusculus.UCSC.mm9.knownGene")
library("TxDb.Mmusculus.UCSC.mm9.knownGene")
TxDbObject <- TxDb.Mmusculus.UCSC.mm9.knownGene

# get the chromosome, and the start and end of our coverage data
chromosome <- as.character(unique(covData[[1]]$chromosome))
start <- as.numeric(min(covData[[1]]$end))
end <- as.numeric(max(covData[[1]]$end))

# define the genomic range
grObject <- GRanges(seqnames=c("chr1"), ranges=IRanges(start=start, end=end))

# Here's a place to start plotting:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
BiocManager::install("GenVisR")
library("GenVisR")
genCov(x=covData, txdb=TxDbObject, gr=grObject, genome=genomeObject, cov_plotType="line")

# Check out what else genCov() can do
?genCov()

# Try playing around with options to make changes to the plot above (color, labels, etc.)

#### Variant Allele Frequency Distributions ####
# https://genviz.org/module-02-r/0002/03/03/ggplot2_exercises/
#install.packages("ggplot2")
library(ggplot2)
# load the dataset
variantData <- read.delim("http://genomedata.org/gen-viz-workshop/intro_to_ggplot2/ggplot2ExampleData.tsv")
variantData <- variantData[variantData$dataset == "discovery",][1:300,]

# Here's a place to start:
ggplot() + geom_violin(data=variantData, aes(x=Simple_name, y=tumor_VAF)) + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  xlab("Sample") + ylab("Variant Allele Fraction")

# Check out what else genCov() can do
?geom_violin()

# Try playing around with options to make changes to the plot above (color, labels, etc.)

#### Copy Number Frequency Plots ####
# https://genviz.org/module-03-genvisr/0003/05/01/cnFreq_GenVisR/
# install.packages("stringr")
library("stringr")

# get locations of all cn files
files <- Sys.glob("cc2_files/*cc2.tsv")

# create function to read in and format data
a <- function(x){
  # read data and set column names
  data <- read.delim(x, header=FALSE)
  colnames(data) <- c("chromosome", "start", "end", "probes", "segmean")
  
  # get the sample name from the file path
  sampleName <- str_extract(x, "H_OM.+cc2")
  sampleName <- gsub(".cc2", "", sampleName)
  data$sample <- sampleName
  
  # return the data
  return(data)
}

# run the anonymous function defined above
cnData <- lapply(files, a)

# turn the list of data frames into a single data frame
cnData <- do.call("rbind", cnData)

# Here's a place to start plotting:
cnFreq(cnData, genome="hg19")
# Alternatively, focus on one gene
layer1 <- geom_vline(xintercept=c(39709170))
cnFreq(cnData, genome="hg19", plotChr="chr17", plotLayer=layer1)

# Check out what else genCov() can do
?geom_violin()

# Try playing around with options to make changes to the plot above (color, labels, etc.)

