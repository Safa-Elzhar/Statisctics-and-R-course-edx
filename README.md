# Statistics Notes in R

My own learning notes and code while working through HarvardX **PH525.1x: Statistics and R**
(Data Analysis for the Life Sciences), taught by Rafael Irizarry.

This repository is a set of worked examples on public data, organized by topic. It is not a set of
course answers. The examples use their own seeds, sample sizes and parameters.

## What is covered

| Topic | Idea |
|---|---|
| Sampling distributions | How a sample mean varies, and how its spread shrinks with sample size |
| t-test by hand | Standard error, t-statistic, p-values (normal vs. t), confidence interval |
| Monte Carlo simulation | Checking the t-distribution by repeating a random experiment |
| Power | How sample size and the significance level affect the chance of detecting an effect |
| Permutation test | A p-value built by shuffling group labels, with no distribution assumptions |
| Robust statistics | Mean and SD vs. median and MAD, Pearson vs. Spearman correlation |
| Rank tests | t-test vs. Wilcoxon test when an outlier is present |
| 2x2 tables | Chi-square and Fisher's exact test |

## How to run

1. Install [R](https://cran.r-project.org/) and, optionally, [RStudio](https://posit.co/download/rstudio-desktop/).
2. Open `statistics_notes.R` and run it from top to bottom (it needs an internet connection to download the data).

Only base R is used, so no extra packages are required.

## Data

- `ChickWeight`: built into R.
- Mouse weights (`femaleMiceWeights.csv`) and birth weights (`babies.txt`): from the public
  [genomicsclass/dagdata](https://github.com/genomicsclass/dagdata) repository.

## Notes

- The script sets a seed before each simulation so results are reproducible.
- Simulation results will change slightly if you change a seed or the number of repeats.

## Acknowledgements

Course: HarvardX PH525.1x, Rafael Irizarry, Harvard T.H. Chan School of Public Health.
Textbook: [Data Analysis for the Life Sciences](https://leanpub.com/dataanalysisforthelifesciences).
