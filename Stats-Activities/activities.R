### Activities for Statistics for Computational Biology Projects Workshop ###
### Erica M. Holdmore ###
### Last updated: March 14, 2024 ###

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

#### Activity 3: Multiple Comparisons ####

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
sum(p.adjust(p.raw, method = 'bonferroni') < 0.05)
which(p.adjust(p.raw, method = 'bonferroni') < 0.05)
# How many significant outcomes are observed using a Bonferroni correction?
# Are all the significant outcomes those that we know are true positives (first 200)?
sum(p.adjust(p.raw, method = 'BH') < 0.05)
which(p.adjust(p.raw, method = 'BH') < 0.05)
# How many significant outcomes are observed using a Benjamini-Hochberg correction?
# Are all the significant outcomes those that we know are true positives (first 200)?

#### Activity 4: Data Visualization ####

# Adapted from Genomic Data Visualization and Interpretation. https://genviz.org/
#install.packages("ggplot2")
library(ggplot2)

## Choose your own adventure! Pick ONE of the example plots below.

#### Gene Coverage Plots ####
# load the dataset
covData <- read.delim("http://genomedata.org/gen-viz-workshop/GenVisR/STAT1_mm9_coverage.tsv")

# We need to reformat the data a bit.
# rename the columns
colnames(covData) <- c("chromosome", "start", "end", "TAC245", "TAC265")

# create a function to split the dataframe into lists of dataframes for each sample
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
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("BSgenome", version = "3.8")
library("BSgenome")
available.genomes()

# install and load the mm9 BSgenome from UCSC
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm9", version = "3.8")
library("BSgenome.Mmusculus.UCSC.mm9")
genomeObject <- BSgenome.Mmusculus.UCSC.mm9

# Install and load a TxDb object for the mm9 genome
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("TxDb.Mmusculus.UCSC.mm9.knownGene", version = "3.8")
library("TxDb.Mmusculus.UCSC.mm9.knownGene")
TxDbObject <- TxDb.Mmusculus.UCSC.mm9.knownGene

# get the chromosome, and the start and end of our coverage data
chromosome <- as.character(unique(covData[[1]]$chromosome))
start <- as.numeric(min(covData[[1]]$end))
end <- as.numeric(max(covData[[1]]$end))

# define the genomic range
grObject <- GRanges(seqnames=c("chr1"), ranges=IRanges(start=start, end=end))

# Here's a place to start plotting:
genCov(x=covData, txdb=TxDbObject, gr=grObject, genome=genomeObject, cov_plotType="line")

#### Variant Allele Frequency Distributions ####
# load the dataset
variantData <- read.delim("http://genomedata.org/gen-viz-workshop/intro_to_ggplot2/ggplot2ExampleData.tsv")
variantData <- variantData[variantData$dataset == "discovery",][1:300,]

# Here's a place to start:
ggplot() + geom_violin(data=variantData, aes(x=Simple_name, y=tumor_VAF)) + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  xlab("Sample") + ylab("Variant Allele Fraction")

ggplot(data=variantData, aes(x=Simple_name, y=tumor_VAF)) + 
  geom_violin(aes(fill=Simple_name)) + geom_jitter(width=.1, alpha=.5) + 
  theme(axis.text.x=element_text(angle=90, hjust=1), legend.position="none") + 
  xlab("Sample") + ylab("Variant Allele Fraction")

#### Copy Number Frequency Plots ####

# install.packages("stringr")
library("stringr")

# get locations of all cn files
files <- Sys.glob("~/Desktop/*cc2.tsv")

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

# focus on one gene using ggplot2
layer1 <- geom_vline(xintercept=c(39709170))
cnFreq(cnData, genome="hg19", plotChr="chr17", plotLayer=layer1)
