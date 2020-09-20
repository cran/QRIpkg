#' @title Quantile Regression Index Score
#' @description
#' The QRI_func() performs quantile regression analysis using age and sex as predictors to calculate
#' the Quantile Regression Index (QRI) score for each individual’s regional brain imaging metrics
#' and then averages across the regional scores to generate an average tissue specific score for each subject.
#' The QRI indicates individual deviations from the expected aging trajectory.
#' Positive QRI indicates accelerated vs. expected aging trajectory while negative QRI indicates delayed aging.
#' The expected aging trajectory is modeled based on sample of controls.
#' @param ID a column name of subject ID in data
#' @param DXcontrol The expected aging trajectory should only be calculated from the controls(i.e. DXcontrol='control==0').
#' If DXcontrol=NULL, the expected aging trajectory will be calculated from the full data.
#' @param predictors column names of predictors (i.e. age, sex)
#' @param resp.range a column range of responses
#' @param data a data frame contains predictors and responses in the quantile regression model.
#' @details
#' The QRI score can be used as an alternative to BrainAge to assess accelerated brain aging by determining an individuals' placement
#' on the expected aging trajectory.A study by Ryan et al (2020) demonstrated that QRI and BrainAge share up to 80% of the variance in
#' both patients and controls. The typical function usage involves calling the QRI function with the following parameters (age, sex) on
#' a list of tissue-specific neuroimaging traits such as regional white matter fractional anisotropy, regional gray matter cortical
#' thickness, or gray matter subcortical volumes. Quantile regression is performed using the controls (DXcontrol='control==0') to generate
#' the normative curves for the 5th, 50th, and 95th percentiles. Then each patient (DXcontrol='control==1') and control’s individual
#' (DXcontrol='control==0') data is compared to the expected aging trajectory. Each regional measure is assigned a score based upon its
#' location: values > 95% of the expected age data are assigned a value of “-1”; values < 5% receive a value of “1”; all others are
#' assigned “0”. The function then averages across the regional data to generate a tissue-specific QRI score (i.e. white matter QRI).
#' @return This function returns the average tissue-specific QRI scores for all subjects.
#' @note
#' The QRI_func() function is developed at the Maryland Psychiatric Research Center, Department of Psychiatry,
#' University of Maryland School of Medicine. This project is supported by NIH R01 EB015611 grant. Please cite our funding if
#' you use this software.
#'
#' Meghann C. Ryan, L. Elliot Hong, Kathryn S. Hatch, Shuo Chen, Krystl Haerian, Jingtao Wang, Eric L. Goldwaser, Xiaoming Du,
#' Bhim M. Adhikari, Heather Bruce, Stephanie Hare, Mark D. Kvarta, Neda Jahanshad, Thomas E. Nichols, Paul M. Thompson,
#' Peter Kochunov. The Additive Impact of Metabolic Disorders and Psychiatric Illnesses on Accelerated Brain Aging. In Review
#' @references
#' Roger Koenker (2020). quantreg: Quantile
#' Regression. R package version 5.61.
#' https://CRAN.R-project.org/package=quantreg
#'
#' R Core Team (2020). R: A language and environment for statistical computing. R
#' Foundation for Statistical Computing, Vienna, Austria. URL
#' https://www.R-project.org/.
#' @examples
#' QRI <- QRI_func(ID='ID', DXcontrol='Control==0', predictors=c('Age','Sex'), resp.range=c(5:6),
#' data=QRIpkg::subcortical)
#' @export

QRI_func <- function(ID, DXcontrol, predictors, resp.range, data) {
  resp.data <- data.frame(data[,resp.range])
  resp.names <- names(data)[resp.range]
  colnames(resp.data) <- resp.names

  if(is.null(DXcontrol)==T){
    control.data <- data}
  else{
    control.data <- subset(data,eval(parse(text=DXcontrol)))}

  QRI.score <- data.frame(ID=data[,which(names(data)==ID)])

  for (j in 1:ncol(resp.data)) {
    rq.formula <- stats::as.formula(paste(resp.names[j], ' ~ ', paste('poly(', paste(predictors, collapse = "+"), ',2,raw=TRUE',')' )))
    model1 <- quantreg::rq(rq.formula,data=control.data, tau=0.05)
    model2 <- quantreg::rq(rq.formula,data=control.data, tau=0.95)

    newdat <- data.frame(data[,match(predictors,colnames(data))])
    colnames(newdat) <- predictors
    CI.pred5  <- stats::predict(model1, newdata = newdat, interval = 'confidence')
    CI.pred95 <- stats::predict(model2, newdata = newdat, interval = 'confidence')

    tmp <- data.frame(tmp.var=resp.data[,j], CI.pred5.lower = CI.pred5[,2], CI.pred95.higher =  CI.pred95[,3])
    for (i in 1:nrow(tmp)) {
      if(is.na(tmp$tmp.var[i]) == TRUE ){tmp$score[i] <- NA}
      else if(tmp$tmp.var[i] < tmp$CI.pred5.lower[i]){tmp$score[i] <- 1}
      else if(tmp$tmp.var[i] > tmp$CI.pred95.higher[i]){tmp$score[i] <- -1}
      else{tmp$score[i] <- 0}
    }

    QRI.score <- cbind(QRI.score,tmp$score)
    }

  colnames(QRI.score) <- c(ID,paste0(resp.names,'.','QRI'))
  average.QRI <- rowMeans(data.frame(QRI.score[,-1]), na.rm = T)
  QRI.score$average.QRI.score <- average.QRI

  return(QRI.score)
}
