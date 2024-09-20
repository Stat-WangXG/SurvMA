# *SurvMA* R Package

## Model Averaging Prediction of Personalized Survival Probabilities

This is an R package for providing predictions of personalized **Surv**ival probabilities based on **M**odel **A**veraging.
- **KEY IDEA**: approximate the conditional survival function using *a weighted average of multiple candidate submodels*. 
- Two scenarios of **Candidate Submodels** are allowed:
  - *Scenario 1*: **Partial Linear** Cox Model;
  - *Scenario 2*: **Time-Varying Coefficient** Cox Model.
- The underlying key statistical techniques we used are "B-spline" and "Brier Score".
- *This R package was contributed by **Mengyu Li**, **Jie Ding** and **Xiaoguang Wang** (Dalian University of Technology).*

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
Visually, the usages of these key functions and part of the their key arguments are summarized via the following figure:
![image](https://github.com/user-attachments/assets/81116d19-85d3-4992-a6a6-57f234279cab)

## Numerical illustrations

### Illustration via a simulated dataset: from partial linear additive Cox model

Load the dataset
```R
data(SimData.APL)
head(SimData.APL,2)
```

Split the data into training and test datasets
```R
set.seed(1)
train.index <- sort(sample(1:200,0.75*200))
sdata.train <- SimData.APL[train.index,]
sdata.test  <- SimData.APL[-train.index,]
```

Fit the data using provided R function *SurvMA.Fit()*
```R
set.seed(1)
sol.SurvMA.PL <- SurvMA.Fit(
  formula = Surv(time,delta) ~ X + U1 + U2 + U3 + U4 + U5 + U6,
  sdata = SimData.APL, submodel = "PL", continuous = 2:7
)
print(sol.SurvMA.PL$weights)
```

Do prediction using provided R function *SurvMA.Predict()*
```R
predict.SurvMA.PL <- SurvMA.Predict(
  object = sol.SurvMA.PL, 
  covariates = sdata.test[,-c(1,2)],
  times = round(quantile(sdata.test$time,c(0.25,0.50,0.75)),2)
)
head(predict.SurvMA.PL$sprobs,2)
```

### Illustration via real dataset: a breast cancer dataset

Load the dataset
```R
data(RealData.ROT)
summary(RealData.ROT$time)
table(RealData.ROT$delta)
```

Plot the Kaplan-Meier curve 
```R
plot(
  survfit(Surv(time,delta) ~ 1, data = RealData.ROT),
  mark.time = TRUE, conf.int = TRUE, lwd=2,
  xlim = c(0,3200), ylim=c(0.4,1),
  xlab="Time (in Days)", ylab="Estimated Survival Probability"
)
```

Test time-varying effects
```R
TVC.Test <- cox.zph(coxph(Surv(time, delta)~., data = RealData.ROT))
print(TVC.Test)
plot(
  TVC.Test, resid = FALSE, lwd = 2,
  xlab = "Time (in Days)",
  ylab = paste("Coefficient for",colnames(RealData.ROT)[1:6])
)
```

Split the data into training and test datasets
```R
set.seed(1)
n <- nrow(RealData.ROT)
train.index <- sort(sample(1:n,0.75*n))
sdata.train <- RealData.ROT[train.index,]
sdata.test  <- RealData.ROT[-train.index,]
```

Fit the data using provided R function *SurvMA.Fit()*
```R
set.seed(1)
sol.SurvMA.ROT <- SurvMA.Fit(
  formula = Surv(time, delta) ~ age + meno + pgr + er + hormon + chemo,
  sdata = sdata.train, submodel = "TVC", continuous = NULL
)
print(sol.SurvMA.ROT$weights)
```

Do prediction using provided R function *SurvMA.Predict()*
```R
predict.SurvMA.ROT <- SurvMA.Predict(
  object = sol.SurvMA.ROT, covariates = 
    sdata.test[,!(colnames(sdata.test) %in% c("time","delta"))],
  times = round(quantile(sdata.test$time,c(0.25,0.50,0.75)))
)
head(predict.SurvMA.ROT$sprobs,2)
```
