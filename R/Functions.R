

###############################
# General Package Description #
###############################

#' SurvMA: Model Averaging Prediction of Personalized Survival Probabilities using R Package SurvMA.
#'
#' This package provides model averaging-based approaches that can be used to predict personalized survival probabilities.
#'
#' @importFrom methods is
#' @importFrom stats as.formula extractAIC na.omit quantile runif
#'
#' @docType package
#' @name SurvMA
NULL

###################################
# Statements of Included Datasets #
###################################

#' @title SimData.APL: A simulated dataset based on a pre-specified partly linear additive Cox model.
#'
#' @aliases SimData.APL
#'
#' @keywords datasets
#'
#' @examples
#'
#' \donttest{
#' # An example of illustrating this dataset can be found in the help page of our
#' #    function \code{SurvMA.Fit()} by typing \code{?SurvMA.Fit()}.
#' }
#'
"SimData.APL"


#' @title SimData.TVC: A simulated dataset based on a pre-specified time-varying coefficients Cox model.
#'
#' @aliases SimData.TVC
#'
#' @keywords datasets
#'
#' @examples
#'
#' \donttest{
#' # An example of illustrating this dataset can be found in the help page of our
#' #    function \code{SurvMA.Fit()} by typing \code{?SurvMA.Fit()}.
#' }
#'
"SimData.TVC"


#' RealData.ROT
#'
#' @title RealData.ROT: A simulated dataset based on a pre-specified time-varying coefficients Cox model.
#'
#' @aliases RealData.ROT
#'
#' @keywords datasets
#'
#' @examples
#'
#' \donttest{
#' # An example of illustrating this dataset can be found in the help page of our
#' #    function \code{SurvMA.Fit()} by typing \code{?SurvMA.Fit()}.
#'
#' # It was originally extracted from the dataset rotteram in R package survival
#' # The specific extractions can be done using the following R commands
#' library(survival)
#' RealData.ROT <- na.omit(rotterdam[
#'   rotterdam$year %in% c(1992,1993),-c(1,2,5,6,7,12,13)
#' ])
#' rownames(RealData.ROT) <- NULL
#' colnames(RealData.ROT)[c(7,8)] <- c("time","delta")
#' }
#'
"RealData.ROT"



##################################################
# Key Functions: SurvMA.Fit() + SurvMA.Predict() #
##################################################


# = Model Fitting
#' @title Model averaging prediction of personalized survival probabilities (model fitting)
#'
#' @param formula a formula expression, of the form \code{response ~ predictors}.
#'   The \code{response} is a \code{Surv} object (from R package "survival") with right censoring.
#'   It is used to specify the included covariates (risk factors).
#'   See the documentation for \code{survreg} and \code{Surv} in R package \code{survival} for details.
#'   The expression to the right of the "~" specifies the covariates.
#' @param sdata a survival dataset (dataframe) in which to interpret the variables named in the \code{formula} and the \code{cureform}.
#' @param submodel a character string defining the groups of candidate submodels, as introduced.
#'   It can be \code{"PL"} for partial linear Cox submodels or \code{"TVC"} for time varying coefficient Cox submodels.
#' @param continuous a vector of integers representing the positions of continuous covariates within \code{predictors} specified in \code{formula}.
#'   If \code{submodel="TVC"} is set, this argument is redundant and the default value \code{NULL} is sufficient.
#' @param control indicates more detailed control of the underlying model averaging fitting procedures.
#'   It is a list of the following three arguments:
#'   \code{K.set} specifies the range of the number of spline basis functions, with the default being \code{K.set=c(5:10)};
#'   \code{criterion} is a character string that specifies the information criterion for choosing the optimal number of B-spline basis functions and
#'      it can be either the default Akaike Information Criterion (\code{criterion="AIC"}) or the Bayesian Information Criterion (\code{criterion = "BIC"});
#'   \code{method} determines the approach to estimate the survival function of censoring time, which can be \code{method="KM"} to
#'      estimate it via the Kaplan-Meier estimator or \code{method = "Cox"} to estimate it via the Cox proportional hazards model.
#'
#' @details
#'   This is a function used to conduct model averaging prediction (model fitting) of personalized survival probabilities.
#'   For obtaining specific predictions of personalized survival probabilities, see another function \code{SurvMA.Predict()}.
#'   The underlying methods are based on the paper titled "Semiparametric model averaging method for survival probability predictions of patients", which has been published in Mengyu Li and Xiaoguang Wang (2023) <doi:10.1016/j.csda.2023.107759>.
#'
#' @return
#'   A list of fitted results that contain not only parameter estimates for all candidate submodels, but also optimal averaging weights (\code{weights}).
#'
#' @examples
#'
#' \donttest{
#' #----------------------------------------------------------#
#' # Basic preparations before running subsequent examples ####
#' #----------------------------------------------------------#
#'
#' rm(list=ls(all=TRUE))
#'
#' ## library necessary packages
#' library(SurvMA)
#' library(survival)
#' #'
#' #--------------------------------------------------------------#
#' # Simulated dataset: from partial linear additive Cox model ####
#' #--------------------------------------------------------------#
#'
#' ## Pre-process the dataset
#'
#' # - load the dataset
#' data(SimData.APL)
#' head(SimData.APL,2)
#'
#' # - split the data into training and test datasets
#' set.seed(1)
#' train.index <- sort(sample(1:200,0.75*200))
#' sdata.train <- SimData.APL[train.index,]
#' sdata.test  <- SimData.APL[-train.index,]
#'
#' ## Fit the dataset via our model averaging method
#'
#' # - fit the data using provided R function SurvMA.Fit
#' set.seed(1)
#' sol.SurvMA.PL <- SurvMA.Fit(
#'   formula = Surv(time,delta) ~ X + U1 + U2 + U3 + U4 + U5 + U6,
#'   sdata = SimData.APL, submodel = "PL", continuous = 2:7
#' )
#' print(sol.SurvMA.PL$weights)
#'
#' # - do prediction using provided R function SurvMA.Predict
#' predict.SurvMA.PL <- SurvMA.Predict(
#'   object = sol.SurvMA.PL,
#'   covariates = sdata.test[,-c(1,2)],
#'   times = round(quantile(sdata.test$time,c(0.25,0.50,0.75)),2)
#' )
#' head(predict.SurvMA.PL$sprobs,2)
#'
#'
#' #-----------------------------------------------------------#
#' # Real dataset: using time-varying coefficient Cox model ####
#' #   - the breast cancer data originally from survival package
#' #-----------------------------------------------------------#
#'
#' ## Pre-process the dataset
#'
#' # - load the dataset
#' data(RealData.ROT)
#' summary(RealData.ROT$time)
#' table(RealData.ROT$delta)
#'
#' # - plot the Kaplan-Meier curve
#' plot(
#'   survfit(Surv(time,delta) ~ 1, data = RealData.ROT),
#'   mark.time = TRUE, conf.int = TRUE, lwd=2,
#'   xlim = c(0,3200), ylim=c(0.4,1),
#'   xlab="Time (in Days)", ylab="Estimated Survival Probability"
#' )
#'
#' # - test time-varying effects
#' TVC.Test <- cox.zph(coxph(Surv(time, delta)~., data = RealData.ROT))
#' print(TVC.Test)
#' par(mfrow=c(2,3))
#' plot(
#'   TVC.Test, resid = FALSE, lwd = 2,
#'   xlab = "Time (in Days)",
#'   ylab = paste("Coefficient for",colnames(RealData.ROT)[1:6])
#' )
#'
#' # - split the data into training and test datasets
#' set.seed(1)
#' n <- nrow(RealData.ROT)
#' train.index <- sort(sample(1:n,0.75*n))
#' sdata.train <- RealData.ROT[train.index,]
#' sdata.test  <- RealData.ROT[-train.index,]
#'
#' ## Fit the dataset via our model averaging method
#'
#' # - fit the data using provided R function SurvMA.Fit
#' set.seed(1)
#' sol.SurvMA.ROT <- SurvMA.Fit(
#'   formula = Surv(time, delta) ~ age + meno + pgr + er + hormon + chemo,
#'   sdata = sdata.train, submodel = "TVC", continuous = NULL
#' )
#' print(sol.SurvMA.ROT$weights)
#'
#' # - do prediction using provided R function SurvMA.Predict
#' predict.SurvMA.ROT <- SurvMA.Predict(
#'   object = sol.SurvMA.ROT, covariates =
#'     sdata.test[,!(colnames(sdata.test) %in% c("time","delta"))],
#'   times = round(quantile(sdata.test$time,c(0.25,0.50,0.75)))
#' )
#' head(predict.SurvMA.ROT$sprobs,2)
#'
#'
#' #----------------------------------------------------------------#
#' # Simulated dataset: from time-varying coefficients Cox model ####
#' #----------------------------------------------------------------#
#'
#' ## Pre-process the dataset
#'
#' # - load the dataset
#' data(SimData.TVC)
#' head(SimData.TVC,2)
#'
#' # - split the data into training and test datasets
#' set.seed(1)
#' train.index <- sort(sample(1:150,0.75*150))
#' sdata.train <- SimData.TVC[train.index,]
#' sdata.test  <- SimData.TVC[-train.index,]
#'
#' ## Fit the dataset via our model averaging method
#'
#' # - fit the data using provided R function SurvMA.Fit
#' set.seed(1)
#' sol.SurvMA.TVC <- SurvMA.Fit(
#'   formula = Surv(time,delta) ~ Z1 + Z2 + Z3 + Z4 + Z5 + Z6,
#'   sdata = sdata.train, submodel = "TVC", continuous = NULL
#' )
#' print(sol.SurvMA.TVC$weights)
#'
#' # - do prediction using provided R function SurvMA.Predict
#' predict.SurvMA.TVC <- SurvMA.Predict(
#'   object = sol.SurvMA.TVC,
#'   covariates = sdata.test[,-c(1,2)],
#'   times = round(quantile(sdata.test$time,c(0.25,0.50,0.75)),2)
#' )
#' head(predict.SurvMA.TVC$sprobs,2)
#' }
#’
#' @export SurvMA.Fit
SurvMA.Fit <- function(
    formula, sdata, submodel = c("PL","TVC"), continuous = NULL,
    control = list(
      K.set = c(5:10),
      criterion = "AIC",
      method = "KM"
    )
){

  ## fit the model with different setups of submodels
  if(submodel=="PL"){

    # re-shape the inputted dataset
    sdata.reshape <- data.frame(
      sdata[,all.vars(formula)[-c(1,2)][continuous],drop=FALSE],
      sdata[,all.vars(formula)[-c(1,2)][-continuous],drop=FALSE],
      time = sdata[,all.vars(formula)[1]],
      delta = sdata[,all.vars(formula)[2]]
    )

    # fit
    sol.ma <- MACOX(
      data = sdata.reshape,
      p = length(continuous),
      K_set = control$K.set,
      compute.S = FALSE,
      method = control$method,
      criterion = control$criterion
    )

    # rename and reformulate some of the elenemts
    sol.ma$weights <- round(as.vector(sol.ma$Weight_IBS),4)
    names(sol.ma$weights) <- paste("submodel ",1:length(continuous),sep="")


  }else if(submodel=="TVC"){

    # re-shape the inputted dataset
    sdata.reshape <- data.frame(
      sdata[,all.vars(formula)[-c(1,2)],drop=FALSE],
      time = sdata[,all.vars(formula)[1]],
      delta = sdata[,all.vars(formula)[2]]
    )

    # fit
    sol.ma <- MAtvcCOX(
      data = sdata.reshape,
      K_set = control$K.set,
      test.plot = FALSE,
      compute.S = FALSE,
      method = control$method,
      criterion = control$criterion
    )

    # rename and reformulate some of the elenemts
    sol.ma$weights <- round(as.vector(sol.ma$MA_weights),4)
    names(sol.ma$weights) <- paste("submodel ",1:(ncol(sdata.reshape)-2),sep="")


  }


  ## include inputted information
  sol.ma$inputs <- list(
    type = "Fit",
    formula = formula,
    submodel = submodel,
    continuous = continuous,
    control = control
  )

  ## output the results
  return(sol.ma)

}


# = Prediction
#' @title Model averaging prediction of personalized survival probabilities (prediction)
#'
#' @param object a list of all outputted results from another main function named \code{SurvMA.Fit()}.
#' @param covariates is a data.frame with rows representing individuals and columns containing the necessary covariates used in the \code{formula} argument of \code{SurvMA.Fit()}.
#' @param times specifies the time points at which survival probabilities will be calculated.
#'
#' @details
#'   This is a function used to conduct model averaging prediction (prediction) of personalized survival probabilities.
#'   For preliminary model fitting process, see another function \code{SurvMA.Fit()}.
#’
#' @return
#'   A list of fitted results that contain, for example, predicted values of personalized survival probabilities.
#'
#' @examples
#'
#' \donttest{
#'
#' # Examples of illustrating the usages of this function can be found in the help
#' #    page of our another function \code{SurvMA.Fit()} by typing \code{?SurvMA.Fit()}.
#'
#' }
#'
#' @export SurvMA.Predict
SurvMA.Predict <- function(object,covariates,times){

  ## fit the model with different setups of submodels
  submodel <- object$inputs$submodel
  if(submodel=="PL"){

    # re-shape the inputted dataset
    covariates.reshape <- data.frame(
      covariates[,all.vars(object$inputs$formula)[-c(1,2)][object$inputs$continuous],drop=FALSE],
      covariates[,all.vars(object$inputs$formula)[-c(1,2)][-object$inputs$continuous],drop=FALSE]
    )

    # predict
    predict.ma <- MApredict(MACOX.object=object,newdata=covariates,t_star=times)

    # rename and reformulate some of the elenemts
    predict.ma$sprobs <- predict.ma$pMA_t_star
    colnames(predict.ma$sprobs) <- paste("time=",times,sep="")
    rownames(predict.ma$sprobs) <- NULL

  }else if(submodel=="TVC"){

    # re-shape the inputted dataset
    covariates.reshape <- data.frame(
      covariates[,all.vars(object$inputs$formula)[-c(1,2)],drop=FALSE]
    )

    # predict
    predict.ma <- MAtvcpredict(MAtdcCOX.object=object,newdata=covariates,t_star=times)

    # rename and reformulate some of the elenemts
    predict.ma$sprobs <- predict.ma$S_MA_t_star
    colnames(predict.ma$sprobs) <- paste("time=",times,sep="")
    rownames(predict.ma$sprobs) <- NULL

  }


  ## include inputted information
  predict.ma$inputs <- list(
    type = "Predict",
    object = object,
    covariates = covariates,
    times=times
  )

  ## output the results
  return(predict.ma)

}


#################################################################
# Function: MACOX() - Partly Linear Additive PH Candidate Model #
#################################################################

MACOX <- function(data, p, K_set=(5:10), compute.S=TRUE, method = "KM", criterion = "AIC") {

  data=na.omit(data)
  data=data[data[,ncol(data)-1]>0,]

  t_ties=as.numeric(names(table(data[,ncol(data)-1]))[table(data[,ncol(data)-1])>1])
  if (length(t_ties)>0)
  {
    row_ties=which(data[,ncol(data)-1] %in% t_ties)
    time_after=data[,ncol(data)-1]
    while (1)
    {
      time_after[row_ties]=data[row_ties,ncol(data)-1]+data[row_ties,ncol(data)-1]*runif(length(row_ties),0.0001,0.001)
      model_c_ties=survival::survfit(survival::Surv(time_after,1-status)~1,data = data.frame(cbind(time_after,status=data[,ncol(data)] )))
      if ( length(model_c_ties$surv) == length(unique(time_after)) )
        break
    }
  }
  data[,ncol(data)-1]=time_after

  n=nrow(data) # sample size
  e=ncol(data)-2 # number of covariates
  z=as.matrix(data[,1:e])
  time=data[,ncol(data)-1] # observed time
  status=data[,ncol(data)] # censoring indicator
  censoring_rate=1-mean(status)
  status=status[order(time)]
  z=z[order(time),]
  time=time[order(time)]
  train.data=data.frame(cbind(z,time,status))

  submodel_all=as.list(rep(0,p)) #candidate models
  names(submodel_all)=paste("submodel",1:p,sep="")
  Coef_m<-as.list(rep(0,p))
  names(Coef_m)=paste("submodel",1:p,sep="")
  Est.eff <- matrix(0,nrow =n,ncol = p)
  K_value<-matrix(0,p,1)
  for (m in 1:p)
  {
    Cre_mavalue<-matrix(0,length(K_set),1)
    Cre_ma=1e+100
    K_ma=K_set[1]
    for ( i_ma in 1:length(K_set) )
    {
      Bspline_ma=splines::bs(train.data[,m],df=K_set[i_ma],degree = 3)
      data_bs<-data.frame(train.data[,1:e],p_ma=Bspline_ma,time=train.data$time,status=train.data$status)
      formula_ma=paste("Surv(time,status)~",paste(colnames(data_bs)[1:e][-m],sep="",collapse = "+"),"+",paste(colnames(data_bs)[(e+1):(e+K_set[i_ma])],sep="",collapse = "+"),sep = "")
      formula_ma=as.formula(formula_ma)
      fit_ma<-survival::coxph(formula_ma,data=data_bs)

      if (criterion == "AIC") {
        Cre_tempma <- extractAIC(fit_ma)[2]
      } else if (criterion == "BIC") {
        Cre_tempma <- extractAIC(fit_ma)[2] + extractAIC(fit_ma)[1] * (log(n) - 2)
      } else {
        stop("criterion must be 'AIC' or 'BIC'")
      }
      Cre_mavalue[i_ma,]<-Cre_tempma
      if ( Cre_tempma < Cre_ma )
      {
        K_ma=K_set[i_ma]
        Cre_ma=Cre_tempma
      }
    }

    pma=splines::bs(train.data[,m],df=K_ma,degree = 3)    ### change the type of splines
    for (i in 1:K_ma) {
      colnames(pma)[i] <- paste0("bs.U.", m, "_", i)
    }

    data_ma<-data.frame(train.data[,1:e],pma,time=train.data$time,status=train.data$status)

    formula_m=paste("Surv(time,status)~",paste(colnames(data_ma)[1:e][-m],sep="",collapse = "+"),"+",paste(colnames(data_ma)[(e+1):(e+K_ma)],sep="",collapse = "+"),sep = "")
    formula_m=as.formula(formula_m)
    fit_m<-survival::coxph(formula_m,data=data_ma)

    coef_m<-as.matrix(fit_m$coefficients)
    colnames(coef_m) <- "Est"

    rc_m=as.matrix(cbind(data_ma[,1:e][,-m],data_ma[,(e+1):(e+K_ma)]))
    Est.eff[,m]<-rc_m%*%coef_m
    K_value[m,]<-K_ma
    Coef_m[[m]]<-coef_m
    submodel_all[[m]]=fit_m
  }


  St <- function(tmm)
  {
    return(SF(tmm,Est.eff[,i],time,status))
  }


  ESSP <- list()
  for(i in 1:p)
  {
    essp <- matrix(unlist(lapply(time,St)), nrow=n,ncol=n)
    ESSP[[i]] <- essp
  }

  if (method == "Cox") {
    fit_cens <- survival::coxph(Surv(time, 1-status) ~., data = train.data,x=TRUE,y=TRUE)
  } else if (method == "KM") {
    fit_cens <- survival::survfit(Surv(time, 1-status) ~ 1, data = train.data)
  } else {
    stop("Unknown method")
  }

  SP_cens <- pec::predictSurvProb(fit_cens,newdata=train.data,times=time)###predict the probablity of censoring time in the train set

  ####model averaging  weights
  GYZ <- diag(SP_cens)

  D_ma <- list()
  d_ma <- list()

  for( im in 1:n){

    It <- 1*(train.data$time>time[im])
    WGt <- ((1-It)*train.data$status)/GYZ+It/SP_cens[,im]
    WGt[is.nan(WGt)] <- 0

    esSP <- matrix(0,n,p)
    for (im1 in 1:p) {
      esSP1 <- ESSP[[im1]][,im]
      esSP[,im1] <- esSP1
    }

    Dmat1_ma <- t(esSP)%*%diag(WGt)%*%esSP
    dvec1_ma <- as.vector(t(esSP)%*%diag(WGt) %*% It )
    D_ma[[im]]<-Dmat1_ma
    d_ma[[im]]<-dvec1_ma

  }

  S1_ma <- 0
  H_ma <- 0

  for(i2 in 1:n){

    S1_ma <- S1_ma+D_ma[[i2]]
    H_ma <- H_ma+d_ma[[i2]]

  }

  Dmat_ma <- S1_ma
  dvec_ma <- H_ma
  Amat_ma <- cbind(rep(1,p), diag(p))
  bvec_ma <- c(1, rep(0,p))

  sol_ma=suppressWarnings(try(quadprog::solve.QP(Dmat_ma,dvec_ma,Amat_ma,bvec_ma,meq=1),silent = T))
  if (is(sol_ma, "try-error") == TRUE)
  {
    cat("There is no solution to the quadratic programming problem!","\n")
    return(NULL)
  }

  w_ma <- as.matrix(sol_ma$solution)  ###IBS weight


  ### maxlik weight
  A <- diag(p)
  B <- rep(0,p)

  res_weight <- maxLik::maxNM(full_cox_loglik, start=c(rep(1/p,p)), constraints=list(ineqA=A, ineqB=B),
                              control=list(printLevel=10^(-8)),
                              new=Est.eff,time=time,status=status)

  w_maxlik<-as.matrix(res_weight$estimate)


  MASim.SP=NULL # each column denotes a MA estimator of the conditional survival function
  MASim.SPMAX=NULL # each column denotes a MA estimator of the conditional survival function
  if (compute.S==T)
  {
    MASim.SP<-matrix(0,nrow = n,ncol = length(time))
    MASim.SPMAX<-matrix(0,nrow = n,ncol = length(time))

    for( k in 1:length(time)){

      t=time[k]

      Est.sp <- matrix(0,nrow =n,ncol = p)       ###the estimated survival probability in the train set given the time point
      for (m2 in 1:p) {
        est.sp <- SF(t,Est.eff[,m2],time,status)
        Est.sp[,m2] <- est.sp
      }

      ESP_ma <- Est.sp%*%w_ma
      ESP_maxlik<-SFM(t,Est.eff,w_maxlik,time,status)

      MASim.SP[,k]<-ESP_ma
      MASim.SPMAX[,k]<-ESP_maxlik

    }
  }
  result <- list(p=p, K_value=K_value, e=e, time=time, status=status, Est.eff=Est.eff, Coef=Coef_m, Weight_IBS=w_ma, Weight_maxlik=w_maxlik, ESP_IBS=MASim.SP, ESP_maxlik=MASim.SPMAX )

  return(result)
}



######################################################################
# Function: MAtvcCOX() - Time-Varying Coefficient PH Candidate Model #
######################################################################

MAtvcCOX<-function(data,K_set=c(5:10),test.plot=F,compute.S=T, method = "KM", criterion = "AIC")
{

  data=na.omit(data)
  data=data[data[,ncol(data)-1]>0,]

  #deal with tied events
  t_ties=as.numeric(names(table(data[,ncol(data)-1]))[table(data[,ncol(data)-1])>1])
  if (length(t_ties)>0)
  {
    row_ties=which(data[,ncol(data)-1] %in% t_ties)
    time_after=data[,ncol(data)-1]
    while (1)
    {
      time_after[row_ties]=data[row_ties,ncol(data)-1]+data[row_ties,ncol(data)-1]*runif(length(row_ties),0.0001,0.001)
      model_c_ties=survival::survfit(Surv(time_after,1-delta)~1,data = data.frame(cbind(time_after,delta=data[,ncol(data)] )))
      if ( length(model_c_ties$surv) == length(unique(time_after)) )
        break
    }
  }
  data[,ncol(data)-1]=time_after

  n=nrow(data) # sample size
  p=ncol(data)-2 # number of covariates
  submodel_all=as.list(rep(0,p)) #candidate models
  names(submodel_all)=paste("submodel",1:p,sep="")
  z=as.matrix(data[,1:p])
  time=data[,ncol(data)-1] # observed time
  delta=data[,ncol(data)] # censoring indicator
  censoring_rate=1-mean(delta)
  delta=delta[order(time)]
  z=z[order(time),]
  time=time[order(time)]
  data_original=data.frame(cbind(z,time,delta))


  ###############################
  vfit = survival::coxph(Surv(time,delta) ~.,data=data_original,iter.max = 40)   ### test the model with time-varying and time-invarying
  if ( test.plot==T )
  {
    plot(survival::cox.zph(vfit),resid=F,col = "red")
  }

  if (method == "Cox") {
    model_c <- survival::coxph(Surv(time,1-delta) ~., data = data_original,x=TRUE,y=TRUE)
    G <- pec::predictSurvProb(model_c,newdata=data_original,times=time)###predict the probablity of censor data in the train set
  } else if (method == "KM") {
    model_c <- survival::survfit(Surv(time,1-delta) ~ 1, data = data_original)
    G=matrix(rep(model_c$surv,n),byrow = F,ncol = n)
  } else {
    stop("Unknown method")
  }

  if (n-nrow(G)!=0)
  {
    G=rbind(G,matrix(rep(G[nrow(G),],n-nrow(G)),nrow = n-nrow(G),byrow = T))
  }
  #############################


  K_MA=rep(0,p)
  Coef=as.list(rep(0,p)) #candidate models
  names(Coef)=paste("submodel",1:p,sep="")
  Bspline_all=lapply(1:p,function(x) x)
  data_new_all=lapply(1:p,function(x) x)
  for (m in 1:p)
  {
    K_set=unique(sort(K_set))
    K=K_set[1]
    AIC_MA=1e+100
    Candidate_error=FALSE
    for (i in 1:length(K_set))
    {
      Bspline=splines::bs(time,df=K_set[i],knots =quantile(time[delta==1],probs = seq(1:(K_set[i]-4))/(K_set[i]-3)),intercept=T,degree = 3)
      data_new<-data.expand(delta, time, z, Bspline,ncol(Bspline))
      data_new<-data.frame(data_new)
      index=which(data_new$time_start==data_new$time_stop)
      if (length(index)!=0)
        data_new<-data_new[-index,]

      formula_m=paste("Surv(time_start,time_stop,delta_new)~",paste(colnames(data_new)[1:p][-m],sep="",collapse = "+"),"+",paste(rep(colnames(data_new)[m],ncol(Bspline)),":",colnames(data_new)[(p+1):(p+ncol(Bspline))],sep="",collapse = "+"),sep = "")
      formula_m=as.formula(formula_m)
      r_temp=suppressMessages(try(survival::coxph(formula_m,data = data_new,control = survival::coxph.control(timefix = F),singular.ok=F,iter.max = 40),silent = T))
      if (is(r_temp, "try-error") == TRUE)
      {break} else
      {
        submodel_m=survival::coxph(formula_m,data = data_new,control = survival::coxph.control(timefix = F),singular.ok=F,iter.max = 40)
      }

      if (criterion == "AIC") {
        AIC_temp=extractAIC(submodel_m)[2]
      } else if (criterion == "BIC") {
        AIC_temp <- extractAIC(submodel_m)[2] + extractAIC(submodel_m)[1] * (log(n) - 2)
      } else {
        stop("criterion must be 'AIC' or 'BIC'")
      }

      if ( AIC_temp < AIC_MA )
      {
        K_MA[m]=K_set[i]
        AIC_MA=AIC_temp
        submodel_all[[m]]=submodel_m
        Bspline_all[[m]]=Bspline
        data_new_all[[m]]=data_new
        Coef[[m]]<-as.matrix(submodel_m$coefficients)
      }
    }

    if ( is(r_temp, "try-error") == TRUE & i==1 )
    {
      Candidate_error=TRUE
      break
    }

  }

  if (Candidate_error==TRUE)
  {
    cat("Candidate model error!","\n")
    return(NULL)
  }

  Dmat=0
  dvec=0
  term_C=0
  for (i in 1:n)
  {
    S_i=vapply(1:p
               ,FUN = function(x)   Breslow.S.m(n,delta,z,Bspline_all[[x]],x,submodel_all[[x]][["coefficients"]][1:(p-1)],submodel_all[[x]][["coefficients"]][p:(p+K_MA[x]-1)],z[i,])
               ,FUN.VALUE = as.double(1:n))
    W_i=diag(vapply(1:n,FUN =function(x) w.t.i.G(time[x],i,time,delta,G),FUN.VALUE = as.numeric(1)))
    I_i=time<time[i]
    Dmat=Dmat+2*t(S_i)%*%W_i%*%S_i
    dvec=dvec+2*t(I_i)%*%W_i%*%S_i
    term_C=term_C+t(I_i)%*%W_i%*%I_i
  }
  Amat=as.matrix(rbind(rep(1,p),diag(p)))
  bvec=c(1,rep(0,p))
  r_temp=suppressWarnings(try(quadprog::solve.QP(Dmat,dvec,t(Amat),bvec=bvec,meq=1),silent = T))
  if (is(r_temp, "try-error") == TRUE)
  {
    cat("There is no solution to the quadratic programming problem!","\n")
    return(NULL)
  }
  answer=quadprog::solve.QP(Dmat,dvec,t(Amat),bvec=bvec,meq=1)
  w_best=answer$solution # the weight vector of candidate models
  S_MA_train=NULL # each column denotes a MA estimator of the conditional survival function
  if (compute.S==T){
    S_MA_train=vapply(1:n,function(i) vapply(1:p,FUN = function(x)
      Breslow.S.m(n,delta,z,Bspline_all[[x]],x,submodel_all[[x]][["coefficients"]][1:(p-1)],submodel_all[[x]][["coefficients"]][p:(p+K_MA[x]-1)],z[i,]),FUN.VALUE = as.double(1:n))%*%w_best,as.double(1:n))
  }
  return( list(data_train=data,z.order=z,Bspline=Bspline_all,time.order=time,delta.order=delta,n=n,p=p,K_value=K_MA,test=survival::cox.zph(vfit),candidate_models=submodel_all,Coef=Coef,MA_weights=w_best) )
}


#####################################################################
# Function: MApredict() - Partly Linear Additive PH Candidate Model #
#####################################################################

MApredict<-function(MACOX.object,newdata,t_star) #newdata can be a p-dimensional vector, n*p dimensional matrix or data frame
{
  if (sum(is.na(newdata))>0)
  {
    cat("Warning:There is NA in newdata!","\n")
    return(NULL)
  }

  p<-MACOX.object$p
  K_value<-MACOX.object$K_value
  e<-MACOX.object$e
  Est.eff<-MACOX.object$Est.eff
  coef_m<-MACOX.object$Coef
  time<-MACOX.object$time
  status<-MACOX.object$status
  Weight_IBS<-MACOX.object$Weight_IBS
  Weight_maxlik<-MACOX.object$Weight_maxlik
  test.data<-newdata

  pEst.eff <- matrix(0,nrow =nrow(test.data),ncol = p)

  for (m in 1:p) {

    pma=splines::bs(test.data[,m],df=K_value[m],degree = 3)

    test.data_ma<-data.frame(test.data[,1:e],pma)

    prc_m=as.matrix(cbind(test.data_ma[,1:e][,-m],test.data_ma[,(e+1):(e+K_value[m])]))
    pEst.eff[,m]<-prc_m%*%coef_m[[m]]

  }

  pMASim.SP<-matrix(0,nrow = nrow(test.data),ncol = length(t_star))
  pMASim.SPMAX<-matrix(0,nrow = nrow(test.data),ncol = length(t_star))

  for( k in 1:length(t_star)){

    t=t_star[k]

    pEst.sp <- matrix(0,nrow =nrow(test.data),ncol = p)      ###the estimated survival probability effects in the test set given the time point
    for (pm2 in 1:p) {

      pest.sp <- preSF(t,Est.eff[,pm2],pEst.eff[,pm2],time,status)
      pEst.sp[,pm2] <- pest.sp
    }

    pESP_ma <- pEst.sp%*%Weight_IBS
    pESP_maxlik<-preSFM(t,Est.eff,pEst.eff,Weight_maxlik,time,status)

    pMASim.SP[,k]<-pESP_ma
    pMASim.SPMAX[,k]<-pESP_maxlik

  }

  ### t=time
  pMASim.SP_test<-matrix(0,nrow = nrow(test.data),ncol = length(time))
  pMASim.SPMAX_test<-matrix(0,nrow = nrow(test.data),ncol = length(time))

  for( k in 1:length(time)){

    t=time[k]

    pEst.sp_test <- matrix(0,nrow =nrow(test.data),ncol = p)      ###the estimated survival probability effects in the test set given the time point

    for (pm2 in 1:p) {

      pest.sp_test <- preSF(t,Est.eff[,pm2],pEst.eff[,pm2],time,status)
      pEst.sp_test[,pm2] <- pest.sp_test
    }

    pESP_ma_test <- pEst.sp_test%*%Weight_IBS
    pESP_maxlik_test<-preSFM(t,Est.eff,pEst.eff,Weight_maxlik,time,status)

    pMASim.SP_test[,k]<-pESP_ma_test
    pMASim.SPMAX_test[,k]<-pESP_maxlik_test

  }

  return( list(n_test=nrow(test.data),t=t_star,obstime=time,pMA_t_star=pMASim.SP,pMA.maxlik_t_star=pMASim.SPMAX,MAXpMA_test=pMASim.SP_test,pMA.maxlik_test=pMASim.SPMAX_test) )

}


##########################################################################
# Function: MAtvcpredict() - Time-Varying Coefficient PH Candidate Model #
##########################################################################

MAtvcpredict <- function(MAtdcCOX.object,newdata,t_star)
{
  if (sum(is.na(newdata))>0)
  {
    cat("Warning:There is NA in newdata!","\n")
    return(NULL)
  }

  t_star=t_star[order(t_star)]
  time=MAtdcCOX.object$time.order
  location_t_star=vapply(1:length(t_star),function(x) sum(t_star[x]>=time),FUN.VALUE = as.numeric(1) )
  n=MAtdcCOX.object$n
  p=MAtdcCOX.object$p
  delta=MAtdcCOX.object$delta.order
  z=MAtdcCOX.object$z.order
  Bspline=MAtdcCOX.object$Bspline
  K=MAtdcCOX.object$K_value
  submodel_all=MAtdcCOX.object$candidate_models
  w_best=MAtdcCOX.object$MA_weights

  if ( is.vector(newdata) & !is.list(newdata) )
  {
    S_MA_test=vapply(1:p,FUN = function(x)
      Breslow.S.m(n,delta,z,Bspline[[x]],x,submodel_all[[x]][["coefficients"]][1:(p-1)],submodel_all[[x]][["coefficients"]][p:(p+K[x]-1)],newdata),FUN.VALUE = as.double(1:n))%*%w_best
    if (max(location_t_star)==0)
    {
      S_MA_t_star=rep(1,length(t_star))
    }
    else
    {
      S_MA_t_star=c(rep(1,sum(location_t_star==0)),S_MA_test[location_t_star[ (sum(location_t_star==0)+1):length(t_star) ]])
    }
    S_MA_t_star=as.matrix(S_MA_t_star)
    rownames(S_MA_t_star)=t_star

    return( list(n_test=1,t_star.order=t_star,S_MA_t_star=S_MA_t_star,time.order=time,delta.order=delta,S_MA_test=as.vector(S_MA_test)) )
  }

  else if ( is.matrix(newdata) | is.data.frame(newdata)  )
  {
    newdata=as.matrix(newdata)
    n_test=nrow(newdata)
    S_MA_test=vapply(1:n_test,function(i) vapply(1:p,FUN = function(x)
      Breslow.S.m(n,delta,z,Bspline[[x]],x,submodel_all[[x]][["coefficients"]][1:(p-1)],submodel_all[[x]][["coefficients"]][p:(p+K[x]-1)],newdata[i,]),FUN.VALUE = as.double(1:n))%*%w_best,as.double(1:n))
    if(max(location_t_star)==0)
    {
      S_MA_t_star=matrix(1,nrow = length(t_star),ncol=n_test)
    }
    else
    {
      S_MA_t_star=rbind( matrix(1,nrow=sum(location_t_star==0),ncol=n_test) , S_MA_test[location_t_star[ (sum(location_t_star==0)+1):length(t_star) ],] )
    }
    rownames(S_MA_t_star)=t_star

    return( list(n_test=n_test,t_star.order=t_star,S_MA_t_star=t(S_MA_t_star),time.order=time,delta.order=delta,S_MA_test=t(S_MA_test)) )
  }

  else
  {
    cat("Warning:The newdata is not a vector or matrix or data frame!","\n")
    return(NULL)
  }
}

#############################
# Other Auxiliary Functions #
#############################

SF <- function(tm,Phi,time,status){

  Phiw <- as.matrix(Phi)
  H0t <- sum(status*sapply(time,function(timei){sum((time>=timei)*exp(Phiw))})^(-1)*(time<=tm))
  Stx <- exp( -H0t*exp(Phiw) )
  return(Stx)

}

# St <- function(tmm)
# {
#   return(SF(tmm,Est.eff[,i],time,status))
# }


full_cox_loglik <- function(w_maxlik,new,time,status){

  Phiw <- as.matrix(new)%*%w_maxlik
  Lw <- sum( status*( Phiw - sapply(time,function(timei){log(sum((time>=timei)*exp(Phiw)))}) ) )
  return(Lw)

}


SFM <- function(tm,hratio,w,time,status){
  Phiw <- as.matrix(hratio)%*%w
  H0t <- sum(status*sapply(time,function(timei){sum((time>=timei)*exp(Phiw))})^(-1)*(time<=tm))
  Stx <- exp( -H0t*exp(Phiw) )
  return(Stx)
}

preSF <- function(tm,Phi,prePhi,time,status){

  Phiw <- as.matrix(Phi)
  H0t <- sum(status*sapply(time,function(timei){sum((time>=timei)*exp(Phiw))})^(-1)*(time<=tm))
  prePhiw <- as.matrix(prePhi)
  Stx <- exp( -H0t*exp(prePhiw) )
  return(Stx)

}

preSFM <- function(tm,hrcov,prehrcov,w,time,status){
  Phiw <- as.matrix(hrcov)%*%w
  H0t <- sum(status*sapply(time,function(timei){sum((time>=timei)*exp(Phiw))})^(-1)*(time<=tm))
  prePhiw<-as.matrix(prehrcov)%*%w
  Stx <- exp( -H0t*exp(prePhiw) )
  return(Stx)
}


Breslow.S.m<-function(n,delta,z, B_spline,m,beta_m,theta_m,covariate) # the Breslow estimator of S_m(t|covariate)
{
  Surv.m.covariate=rep(0,n)
  if (ncol(z)==2)
  {
    hazard0<-vapply(c(1:n),FUN =function(x) sum( exp( z[x:n,-m]*beta_m + z[x:n,m]*as.numeric(theta_m%*%B_spline[x,]) )),FUN.VALUE =numeric(1) )
  }
  else
  {
    hazard0<-vapply(c(1:n),FUN =function(x) sum( exp( z[x:n,-m]%*%beta_m + z[x:n,m]*as.numeric(theta_m%*%B_spline[x,]) )),FUN.VALUE =numeric(1) )
  }
  hazard0<-delta/hazard0
  hazard0[delta==0]=0
  hazard0[hazard0==Inf]=1e+300
  Surv.m.covariate<-exp(-cumsum( as.vector(hazard0) *  as.vector( exp( as.numeric(covariate[-m]%*%beta_m) + covariate[m]*B_spline[,]%*%theta_m ) ) ))
  Surv.m.covariate
}


w.t.i.G<-function(t,i,time,delta,G)
{
  if ( time[i]<=min(time) & t<min(time) )
    weight.t.i=((time[i]<=t)*delta[i])/(1)+(time[i]>t)/(1)
  if ( time[i]<=min(time) & !t<min(time) )
    weight.t.i=((time[i]<=t)*delta[i])/(1)+ifelse(time[i]<=t,0,(time[i]>t)/(G[rev(which(t>=time))[1],i]))
  if ( !time[i]<=min(time) & t<min(time) )
    weight.t.i=((time[i]<=t)*delta[i])/(G[rev(which(time[i]>time))[1],i])+(time[i]>t)/(1)
  if ( !time[i]<=min(time) & !t<min(time) )
    weight.t.i=((time[i]<=t)*delta[i])/(G[rev(which(time[i]>time))[1],i])+ifelse(time[i]<=t,0,(time[i]>t)/(G[rev(which(t>=time))[1],i]))
  return(weight.t.i)
}


data.expand<- function(delta2, time2, z2, bs7_2,K) {

  n=length(time2)
  delta_new = numeric(n*(n+1)/2)
  time_start = numeric(n*(n+1)/2)
  time_stop = numeric(n*(n+1)/2)
  z_new = matrix(0, n*(n+1)/2,ncol(z2))
  bs_new = matrix(0, n*(n+1)/2,K)
  index = numeric(n*(n+1)/2)
  key=0
  i=0
  repeat {
    i=i+1

    j=0
    repeat{
      j=j+1

      delta_new[key+j] = delta2[i]*(i==j)
      if (j==1) {
        time_start[key+j]=0
      }

      if (j>1) {
        time_start[key+j]=time2[j-1]
      }

      time_stop[key+j]=time2[j]
      z_new[key+j,]=z2[i,]
      bs_new[key+j,]=bs7_2[j,]
      index[key+j]=j
      if(j==i) break
    }

    key=key+j
    if(i==n) break
  }
  list(z_new=z_new, bs_new=bs_new, delta_new=delta_new, time_start=time_start, time_stop=time_stop, index=index)
}








