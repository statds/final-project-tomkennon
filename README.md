# Finding an Ultimate Limit for an NBA Player's Shooting Percentage
This is a repository for my Undergraduate Honors Thesis at UConn.  This analysis was also used for a graduate course I took at UConn called STAT 6494: Data Science in Action.  The repository includes a project proposal, interim project presentation, final project report and corresponding final project report presentation.


## Brief Description of Project

Extreme value theory is used to estimate the ultimate upper or lower limit for an NBA season's league leading player's shooting percentage (free throw, 2 point field goal, and 3 point field goal).  The limits are found using the generalized extreme value distribution with parameters optimized using the Nelder-Mead method maximizing the loglikelihood.  Two different techniques are applied in this project to finding an optimal generalized extreme value distribution location parameter $\mu$ including a constant and a Gompertz curve.  A 95% bootstrap confidence interval is calculated for these limits.  Kolmogorov-Smirnov tests and a Score test are run as well through the bootstrapped datasets to evaluate goodness of fits.


## Files
The main files in this repository are as follows:
- Final Report Source: [RMD](https://github.com/statds/final-project-tomkennon/blob/master/Final%20Project%20Report/Final%20Paper/template.Rmd)
- Final Report Paper: [PDF](https://github.com/statds/final-project-tomkennon/blob/master/Final%20Project%20Report/Final%20Paper/template.pdf)

## File Types
This {.PDF} paper is generated from an R markdown {.Rmd} file.  The datasets {.csv} and image files {.png} used are included in the repository to allow for reproducible code.
