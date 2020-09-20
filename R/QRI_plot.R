#' @title Quantile Regression Index plot
#' @description
#' The QRI_plot() is used to plot Quantile Regression Index (QRI) and generate the normative curves for individual's regional brain imaging metrics.
#' @param x the x coordinate for the QRI plot
#' @param y the y coordinate for the QRI plot
#' @param xlab the label for the x coordinate
#' @param ylab the label for the y coordinate
#' @param DXcontrol the expected aging trajectory. It should only be calculated from the controls(i.e. DXcontrol='control==0').
#' If DXcontrol=NULL, the expected aging trajectory will be calculated from the full data.
#' @param data a data frame contains the predictor(x coordinate), response(y coordinate) and control(DXcontrol) in the quantile regression model.
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
#' @return This function returns a plot for individual measurements.
#' @note
#' The QRI_plot() function is developed at the Maryland Psychiatric Research Center, Department of Psychiatry,
#' University of Maryland School of Medicine. This project is supported by NIH R01 EB015611 grant. Please cite our funding if
#' you use this software.
#'
#' Meghann C. Ryan, L. Elliot Hong, Kathryn S. Hatch, Shuo Chen, Krystl Haerian, Jingtao Wang, Eric L. Goldwaser, Xiaoming Du,
#' Bhim M. Adhikari, Heather Bruce, Stephanie Hare, Mark D. Kvarta, Neda Jahanshad, Thomas E. Nichols, Paul M. Thompson,
#' Peter Kochunov. The Additive Impact of Metabolic Disorders and Psychiatric Illnesses on Accelerated Brain Aging. In Review
#' @references
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
#'
#' Roger Koenker (2020). quantreg: Quantile
#' Regression. R package version 5.61.
#' https://CRAN.R-project.org/package=quantreg
#'
#' R Core Team (2020). R: A language and environment for statistical computing. R
#' Foundation for Statistical Computing, Vienna, Austria. URL
#' https://www.R-project.org/.
#' @examples
#' QRIplot <- QRI_plot(x='Age',y='Ventricle', xlab='Age', ylab='Ventricle', DXcontrol='Control==0',
#' data=QRIpkg::subcortical)
#' @export

QRI_plot <- function(x, y, xlab, ylab, DXcontrol, data){
  column_num <- match(c(x,y),colnames(data))
  mydata <- data.frame(data[,column_num])
  colnames(mydata) <- c('x','y')

  if(is.null(DXcontrol)==T){control_data <- data[,column_num]}
  else{control_data <- subset(data,eval(parse(text=DXcontrol)))[,column_num]}
  colnames(control_data) <- c('x','y')

  rq_formula <- stats::as.formula(paste('y', ' ~ ', paste('poly(','x', ',2,raw=TRUE',')' )))
  model1 <- quantreg::rq(rq_formula,data=control_data, tau=0.05)
  model2 <- quantreg::rq(rq_formula,data=control_data, tau=0.5)
  model3 <- quantreg::rq(rq_formula,data=control_data, tau=0.95)

  newdat=data.frame(x=mydata[,1])
  CI_5  <- data.frame(stats::predict(model1, newdata = newdat, interval = 'confidence'),x=newdat)
  CI_50 <- data.frame(stats::predict(model2, newdata = newdat, interval = 'confidence'),x=newdat)
  CI_95 <- data.frame(stats::predict(model3, newdata = newdat, interval = 'confidence'),x=newdat)

  mydata <- data.frame(mydata,
                       CI_5_lower=CI_5[,2],CI_5_higher=CI_5[,3],
                       CI_50_lower=CI_50[,2], CI_50_higher=CI_50[,3],
                       CI_95_lower=CI_95[,2],CI_95_higher=CI_95[,3])

  for (i in 1:nrow(mydata)) {
    if(is.na(mydata[i,2]) == T  ){mydata$color[i] <- NA}
    else if(mydata[i,2] < mydata[i,3] | mydata[i,2] > mydata[i,8]){mydata$color[i] <- '1'}
    else{mydata$color[i] <- '2'}}

  myplot <- ggplot2::ggplot(data=mydata) +
    ggplot2::geom_point(ggplot2::aes(x=x,y=y,colour = color),size=2) +
    ggplot2::geom_line(ggplot2::aes(x=x,y=fit),data=CI_5,colour = c('red'), size=1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = CI_5_lower, ymax = CI_5_higher, x = x),
                         alpha = 0.2,linetype=2, colour = "red", fill = "red") +
    ggplot2::geom_line(ggplot2::aes(x=x,y=fit),data=CI_50,colour = c('blue'), size=1)  +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = CI_50_lower, ymax = CI_50_higher, x = x),
                         alpha = 0.2,linetype=2, colour = "blue", fill = "blue") +
    ggplot2::geom_line(ggplot2::aes(x=x,y=fit),data=CI_95,colour = c('green'), size=1)  +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = CI_95_lower, ymax = CI_95_higher, x = x),
                         alpha = 0.2,linetype=2, colour = "green", fill = "green") +
    ggplot2::scale_colour_manual(values = c('black','grey')) +
    ggplot2::labs(x = xlab,y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.text.x = ggplot2::element_text(hjust = 0.4,size = 20),
                   axis.text.y = ggplot2::element_text(hjust = 0.4,angle = 90,size = 20),
                   axis.title.x= ggplot2::element_text(vjust = -2,size=22),
                   axis.title.y= ggplot2::element_text(vjust = 4,size=22),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.margin = ggplot2::unit(c(2,1,1,1),"cm"),
                   legend.position = "none")
  return(myplot)
}

utils::globalVariables(c("fit", "color", "CI_5_lower", "CI_5_higher", "CI_50_lower",
                         "CI_50_higher", "CI_95_lower", "CI_95_higher"), package="QRIpkg", add=F)
