# Statistics notes in R
# Personal learning notes from HarvardX PH525.1x (Statistics and R).
#
# Topics: sampling distributions, the CLT, t-tests, Monte Carlo simulation,
# power, permutation tests, robust statistics, rank tests and 2x2 tables.
#
# Data: ChickWeight is built into R. The mouse and birth-weight data are
# downloaded from the public genomicsclass/dagdata repository on GitHub.
# Only base R is needed. Run the script from top to bottom.


# ---- 0. Data ----------------------------------------------------------------

base_url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/"
mice   <- read.csv(paste0(base_url, "femaleMiceWeights.csv"))
babies <- read.table(paste0(base_url, "babies.txt"), header = TRUE)

chow <- mice$Bodyweight[mice$Diet == "chow"]   # control diet
hf   <- mice$Bodyweight[mice$Diet == "hf"]     # high-fat diet

bwt_ns <- babies$bwt[babies$smoke == 0]        # birth weight (oz), non-smoking mothers
bwt_s  <- babies$bwt[babies$smoke == 1]        # birth weight (oz), smoking mothers

# Population standard deviation (divides by n, unlike sd(), which divides by n - 1)
popsd <- function(x) sqrt(mean((x - mean(x))^2))


# ---- 1. Sampling distribution of the mean ---------------------------------

# A sample mean changes from sample to sample. Its spread (the standard error)
# should shrink like sigma / sqrt(n).
set.seed(101)
sample_means <- function(x, n, B = 5000) replicate(B, mean(sample(x, n)))

m10 <- sample_means(bwt_ns, 10)
m40 <- sample_means(bwt_ns, 40)

par(mfrow = c(1, 2))
rng <- range(c(m10, m40))
hist(m10, xlim = rng, main = "Means of samples of 10", xlab = "Sample mean")
hist(m40, xlim = rng, main = "Means of samples of 40", xlab = "Sample mean")
par(mfrow = c(1, 1))

# Simulated spread vs. the formula sigma / sqrt(n)
c(simulated_10 = sd(m10), formula_10 = popsd(bwt_ns) / sqrt(10),
  simulated_40 = sd(m40), formula_40 = popsd(bwt_ns) / sqrt(40))

# Is the distribution of the means roughly normal? Points near the line say yes.
qqnorm(m10); qqline(m10)


# ---- 2. Two-sample t-test by hand -------------------------------------------

# Difference in mean weight: high-fat minus chow
diff_obs <- mean(hf) - mean(chow)
se       <- sqrt(var(hf) / length(hf) + var(chow) / length(chow))
t_obs    <- diff_obs / se

# Two-sided p-value from the normal approximation (CLT)
p_normal <- 2 * (1 - pnorm(abs(t_obs)))

# Two-sided p-value from the t-distribution (better for small samples)
tt       <- t.test(hf, chow)
p_t      <- tt$p.value

c(t = t_obs, p_normal = p_normal, p_t = p_t)

# 95% confidence interval using the t-distribution
diff_obs + c(-1, 1) * qt(0.975, df = tt$parameter) * se


# ---- 3. Monte Carlo check of the t-distribution -----------------------------

# If the data are normal, t = sqrt(n) * mean(x) / sd(x) follows a t-distribution
# with n - 1 degrees of freedom. Check this by simulation.
set.seed(202)
n  <- 8
ts <- replicate(10000, {
  x <- rnorm(n)
  sqrt(n) * mean(x) / sd(x)
})

# Proportion of simulated t-values above 2.5 vs. t and normal theory
c(simulated = mean(ts > 2.5),
  t_theory  = 1 - pt(2.5, df = n - 1),
  normal    = 1 - pnorm(2.5))

hist(ts, breaks = 60, freq = FALSE, xlim = c(-6, 6),
     main = "Simulated t-statistics (n = 8)", xlab = "t")
curve(dt(x, df = n - 1), add = TRUE, lwd = 2)
curve(dnorm(x), add = TRUE, lty = 2)


# ---- 4. Power by simulation -------------------------------------------------

# Power = how often we reject the null when a real difference exists.
# Here the "real difference" is smoking vs. non-smoking birth weight.
power_sim <- function(N, alpha = 0.05, B = 2000) {
  rejected <- replicate(B, {
    a <- sample(bwt_ns, N)
    b <- sample(bwt_s,  N)
    t.test(b, a)$p.value < alpha
  })
  mean(rejected)
}

set.seed(303)
Ns <- c(10, 25, 50, 100)
power_table <- data.frame(
  N          = Ns,
  alpha_0.05 = sapply(Ns, power_sim, alpha = 0.05),
  alpha_0.01 = sapply(Ns, power_sim, alpha = 0.01)
)
power_table
# Larger N raises power. A stricter alpha lowers it.


# ---- 5. Permutation test ----------------------------------------------------

# No distribution assumptions: shuffle the group labels many times to build the
# null distribution of the difference in means.
set.seed(404)
N <- 12
a <- sample(bwt_s,  N)
b <- sample(bwt_ns, N)
obs <- mean(a) - mean(b)

pooled <- c(a, b)
B <- 2000
null <- replicate(B, {
  shuffled <- sample(pooled)
  mean(shuffled[1:N]) - mean(shuffled[(N + 1):(2 * N)])
})

# Two-sided p-value, counting the observed result itself
(sum(abs(null) >= abs(obs)) + 1) / (B + 1)

hist(null, breaks = 40, main = "Permutation null distribution",
     xlab = "Difference in means")
abline(v = c(-abs(obs), abs(obs)), lty = 2, lwd = 2)


# ---- 6. Robust statistics ---------------------------------------------------

# One extreme value can wreck the mean and the standard deviation,
# but the median and the MAD barely move.
w10 <- ChickWeight$weight[ChickWeight$Time == 10]

summaries <- function(x) c(mean = mean(x), sd = sd(x), median = median(x), mad = mad(x))
rbind(clean        = summaries(w10),
      with_outlier = summaries(c(w10, 1000)))

# Correlation: Pearson is sensitive to outliers, Spearman uses ranks
wide <- reshape(as.data.frame(ChickWeight), idvar = c("Chick", "Diet"),
                timevar = "Time", direction = "wide")
wide <- na.omit(wide)
early <- wide$weight.10
late  <- wide$weight.21

c(pearson_clean     = cor(early, late),
  pearson_outlier   = cor(c(early, 2000), c(late, 2000)),
  spearman_clean    = cor(early, late, method = "spearman"),
  spearman_outlier  = cor(c(early, 2000), c(late, 2000), method = "spearman"))


# ---- 7. t-test vs. rank test (Wilcoxon) with an outlier ---------------------

d1 <- ChickWeight$weight[ChickWeight$Time == 10 & ChickWeight$Diet == "1"]
d4 <- ChickWeight$weight[ChickWeight$Time == 10 & ChickWeight$Diet == "4"]

# wilcox.test() warns about ties and uses a normal approximation. That is fine here.
p_values <- function(x, y) {
  c(t_test    = t.test(x, y)$p.value,
    wilcoxon  = suppressWarnings(wilcox.test(x, y)$p.value))
}

rbind(clean        = p_values(d1, d4),
      with_outlier = p_values(c(d1, 500), d4))
# The t-test is thrown off by the outlier. The rank test is much more stable.


# ---- 8. Association in a 2x2 table ------------------------------------------

# Made-up example counts: genotype (rows) by case/control status (columns)
tab <- matrix(c(30, 20, 18, 32), nrow = 2,
              dimnames = list(genotype = c("AA", "Aa"),
                              status   = c("case", "control")))
tab
chisq.test(tab)    # large-sample test (uses a continuity correction in the 2x2 case)
fisher.test(tab)   # exact test, safer with small counts
