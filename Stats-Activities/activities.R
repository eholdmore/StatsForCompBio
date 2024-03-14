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
