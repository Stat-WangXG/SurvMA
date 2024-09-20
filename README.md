# *SurvMA* R Package

## Model Averaging Prediction of Personalized Survival Probabilities

This is an R package for providing predictions of personalized **Surv**ival probabilities based on **M**odel **A**veraging approaches.
- The underlying key statistical techniques we used are "B-spline" and "Brier Score".
- Two scenarios of candidate submodels are allows:
-- *Scenario 1*: **partial linear** Cox model;
-- *Scenario 2*: **time-varying coefficient** model.

## Package description and included main functions

Installation of this package can be done locally after downloading the package manually from this github website. We will also upload this package to the Comprehensive R Archive Network (CRAN) so that it can be downloaded as a standard R package. Currently, it can be loaded using R command
```R
devtools::install_github("Stat-WangXG/SurvMA")
library(SurvMA)
```

The main functions included in our R package is *SurvMA.Fit()* and *SurvMA.Predict()*. 
They can be called via the following R commands:
- **SurvMA.Fit()**: obtaining parameters for candidate submodels and optimal model averaging weights
```R
SurvMA.Fit(formula, sdata, submodel = c("PL","TVC"), continuous = NULL, control)
```
- **SurvMA.Predict()**: predicting personalized survival probabilities
```R
SurvMA.Predict(object, covariates, times)
```
We refer to its help page for more detailed explanations of the corresponding arguments (typing *?SurvMA.Fit()*). 

