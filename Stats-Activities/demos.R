### Demos for Statistics for Computational Biology Projects Workshop ###
### Erica M. Holdmore ###
### Last updated: June 16, 2025 ###

#### Setup Environment ####
install.packages("renv")
renv::restore()

#### Outlier Analysis ####
# install.packages("ggplot2")
library(ggplot2)

# Generate normal data
set.seed(123) # Setting seed for reproducibility
normal_data <- rnorm(100, mean = 0, sd = 1)
# Attach an outlier on the right
right_outlier <- c(normal_data, 8, 9)
# Plot simple chart
plot(right_outlier, ylab = "Log2 CN Ratio")
boxplot(right_outlier, horizontal = TRUE, xlab = "Log2 CN Ratio")

## IQR and Tukey method
# computing IQR and outliers
Q1 <- quantile(right_outlier, 0.25)
Q3 <- quantile(right_outlier, 0.75)
IQR <- IQR(right_outlier)
# define "acceptable" range for values
multiplier <- 1.5
lower_bound <- Q1 - (IQR * multiplier)
upper_bound <- Q3 + (IQR * multiplier)

# Create a dataframe for ggplot
data_frame <- data.frame(value = right_outlier)
# Create a scatter plot to show outliers
ggplot(data_frame, aes(x = seq_along(value), y = value)) +
  geom_point(aes(color = (value < lower_bound | value > upper_bound)), size = 3) +
  scale_color_manual(values = c("black", "red"), name = "Outlier") +
  geom_hline(yintercept = lower_bound, linetype = "dashed", color = "red") +
  geom_hline(yintercept = upper_bound, linetype = "dashed", color = "red") +
  labs(x = "Index", y = "Log2 CN Ratio") +
  theme_minimal()

## Significance Tests (Grubbs' and Rosner)
# perform significance tests
library(outliers)
library(EnvStats)
grubbs.test(right_outlier)
# k= number of suspected outliers
rosnerTest(right_outlier,k = 5)

## Unsupervised Clustering (DBSCAN)
# Adapted from D. Diachkov https://blog.devgenius.io/outlier-detection-in-r-with-dbscan-dd54ebc5404d
#install.packages("dbscan")
library(dbscan)

# create data
set.seed(1) # Ensures reproducibility
time <- seq(1, 100, by = 1)
temperature <- sin(time / 8) * 10 + rnorm(100, mean = 20, sd = 2) # Simulated temperature readings
fake_temperature <- c(temperature[1:49],0,temperature[51:99],50)
data <- data.frame(time, fake_temperature)

# plot data
plot(data$time, data$fake_temperature, 
     main = "Temperature Readings Over Time", 
     xlab = "Time (hours)", 
     ylab = "Temperature (Celsius)", 
     pch = 20)

# run DBSCAN
# eps: How densely clustered do points have to be
#      to be considered "neighbors"?
# minPts: How many "neighbors" does a point need to
#         have in order to be a "core point"?
dbscan_result <- dbscan(data, eps = 6, minPts = 4)
dbscan_result

# plot clusters
par(mar=c(5.0, 4.0, 4.0, 8.0), xpd=TRUE)
plot(data$time, data$fake_temperature, col = c("red","black")[dbscan_result$cluster +1], 
     main = "DBSCAN Clustering Results", 
     xlab = "Time (hours)", 
     ylab = "Temperature (Celsius)", 
     pch = 20)

legend("topleft", 
       legend = c("Outlier", "Normal"),
       fill = c("red","black"))




#### Zero-Inflated Data ####
# Adapted from Fukami and Henderson https://fukamilab.github.io/BIO202/04-C-zero-data.html
#install.packages("lattice")
#install.packages("MASS")
#install.packages("pscl")
#install.packages("lmtest")
library(lattice)
library(MASS)
require(pscl) # can also use package ZIM for zero-inflated models
library(lmtest)

# simulate poisson data
set.seed(123) # Setting seed for reproducibility
poisson_data <- rpois(600, lambda=10)
# add some zeros
zero_inflated <- c(poisson_data, rep(0,200))
hist(zero_inflated, breaks=20)
# check what percent of your data are zeros
100*sum(zero_inflated == 0)/length(zero_inflated)

# simulate a predictor
predictor <- zero_inflated + rnorm(length(zero_inflated),mean=3,sd=1)
data1 <- as.data.frame(cbind(zero_inflated,predictor))

## Poisson GLM
M1 <- glm(zero_inflated ~ predictor,
          family = 'poisson',
          data = data1)
summary(M1)

## Zero-inflated Poisson GLM
M2 <- zeroinfl(zero_inflated ~ predictor | ## Predictor for the Poisson process
                 predictor, ## Predictor for the Bernoulli process;
               dist = 'poisson',
               data = data1)
summary(M2)

# Use likelihood ratio test to see which model fits data better
lrtest(M1, M2)
#### Regression (Linear & Logistic) ####
#install.packages("ggplot2")
library(ggplot2)

## Linear Regression
# simulate normal data
set.seed(123) # Setting seed for reproducibility
normal_data <- rnorm(100, mean = 100, sd = 10)
# simulate a predictor
predictor <- normal_data + rnorm(length(normal_data),mean=10,sd=5)
data1 <- as.data.frame(cbind(normal_data,predictor))

# fit linear model
fit1 <- lm(normal_data ~ predictor, data = data1)
summary(fit1)
fit1$coefficients

# plot data and model
# Note that ggplot2 has a built in function for plotting lm
ggplot(data1, aes(x = predictor, y = normal_data)) + 
  geom_point(alpha = 0.75) + 
  geom_smooth(method="lm", color = 'red') +
  xlab("Histone Modification Score") + ylab("Gene Expression")

## Logistic Regression
# simulate binary data
set.seed(123) # Setting seed for reproducibility
binom_data <- rbinom(100, 1, prob = 0.4)
# simulate a predictor
predictor <- binom_data + rnorm(length(binom_data),mean=2,sd=1)
data1 <- as.data.frame(cbind(binom_data,predictor))

# fit logistic model
fit2 <- glm(binom_data ~ predictor, data=data1)
summary(fit2)
fit2$coefficients

# plot data and model
# Note that ggplot2 has a built in function for plotting glms of many families
ggplot(data1, aes(x = predictor, y = binom_data)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), color = "red") +
  xlab("HRD Score") + ylab("BRCA Status")

