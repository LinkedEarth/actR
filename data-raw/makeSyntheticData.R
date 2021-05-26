#library(arima)
#library(ggpubr)
library(ggthemes)
library(geoChronR)
library(ggplot2)

#' Create a linear ramp around the midpoint of a time series
#'
#' @param lngth length of the time axis
#' @param width the width of the ramp. width  = 0 results in a step jump at the midpoint. width = lngth/2 results in a smooth linear ramp
#' @importFrom fBasics Ramp
#' @return values of the function (between 0 and 1)
#' @export

linear.ramp <-function(lngth,width){
  if (width>lngth/2) {stop("Ramp width cannot exceed half of the length")}
  xs = seq(lngth)/lngth-1/2
  ramp = 4*fBasics::Ramp(xs,a = -width/lngth)
  ramp[xs>width/lngth] = 1.0
  return(ramp)
}



ar1.noise <- function(lngth,g,sigma){
  ar1=sigma/sqrt(1-g^2)*arima.sim(model=list(g,0,0),n=lngth)
  return(ar1)
}


# global variables
theme_set(ggthemes::theme_hc(style = "darkunica"))
lngth <- 400
time <- 1.0*seq(lngth)

# case 1: abrupt shift
signal1 <- linear.ramp(lngth,0)
noise_lo <- ar1.noise(lngth,0.8,0.2)
noise_hi <- ar1.noise(lngth,0.8,0.8)

df1 <- data.frame(signal=signal1,noise_hi,noise_lo)
c1a <- ggplot(data=df1) +
  geom_line(aes(x=time,y=signal1+noise_lo),color='orange') +
  geom_line(aes(x=time,y=signal1),color='white') +
  ylab("y(t)") + xlab(NULL) + ylim(-3,4) +
  ggtitle(expression('(1a) Abrupt Shift + AR(1), '*sigma*' = 0.2)'))
  # TODO: put y in LiPD-TIBBLE


c1b <- ggplot(data=df1) +
  geom_line(aes(x=time,y=signal1+noise_hi),color='orange') +
  geom_line(aes(x=time,y=signal1),color='white') +
  ggtitle(expression('(1b) Abrupt Shift + AR(1), '*sigma*' = 0.8)')) +
  ylab('y(t)') + xlab(NULL) + ylim(-3,4)
# TODO: put y in LiPD-TIBBLE


# case 2: noisy jump with trend
signal2 <- linear.ramp(lngth,50)
df2 <- data.frame(signal=signal2,noise_hi,noise_lo)
c2a <- ggplot(data=df2) +
  geom_line(aes(x=time,y=signal2+noise_lo),color='orange') +
  geom_line(aes(x=time,y=signal2),color='white') +
  ylab("y(t)") + xlab(NULL) + ylim(-3,4) +
  ggtitle(expression("(2a) Gradual shift + AR(1), "*sigma*" = 0.2"))
# TODO: put y in LiPD-TIBBLE

c2b <- ggplot(data=df2) +
  geom_line(aes(x=time,y=signal2+noise_hi),color='orange') +
  geom_line(aes(x=time,y=signal2),color='white') +
  ylab("y(t)") + xlab(NULL) + ylim(-3,4) +
  ggtitle(expression("(2b) Gradual shift + AR(1), "*sigma*" = 0.8"))

# TODO: put y in LiPD-TIBBLE

#plot22 <- ggarrange(c1a,c1b,c2a,c2b)
#show(plot22)

# case 3: noisy jump with age errors
bam.model <- list(ns = 1000, name = "bernoulli", param = 0.05)
timeMat <- geoChronR::simulateBam(X = matrix(1,nrow = length(time)),
                                  t = as.matrix(time),
                                  model = bam.model,
                                  ageEnsOut = TRUE)$ageEns

c3a <- geoChronR::plotTimeseriesEnsRibbons(X = timeMat, Y = signal1+noise_lo, color.line = 'orange') +
  geom_line(aes(x=time,y=signal1),color='black') + geoChronRPlotTheme(style="darkunica") +
  ggtitle('3a) 1a + age uncertainties') +
  ylab("y(t)") + xlab("t") + ylim(-3,4)

# TODO: put Y in LiPD-TIBBLE

c3b <- geoChronR::plotTimeseriesEnsRibbons(X = timeMat, Y = signal1+noise_hi, color.line = 'orange') +
  geom_line(aes(x=time,y=signal1),color='black') + geoChronRPlotTheme(style="darkunica") +
  ggtitle('3b) 1b + age uncertainties') +
  ylab("y(t)") + xlab("t") + ylim(-3,4)

# TODO: put Y in LiPD-TIBBLE

#show(c3a)
grobs = rbind(c1a,c1b,c2a,c2b,c3a,c3b)
plots <- gridExtra::grid.arrange(c1a,c1b,c2a,c2b,c3a,c3b, nrow = 3)




#Create lipd-ts-tibble
#this will be a 4 row tibble, for the four examples, that includes age uncertainty in all.
row1 <- tibble::tibble(time = list(time),
                       timeUnits = "CE",
                       timeEnsemble = list(timeMat),
                       timeVariableName = "year",
                       archiveType = "pseudoproxy",
                       dataSetName = "pseudoShiftData_400yr",
                       geo_latitude = NA,
                       geo_longitude = NA,
                       geo_elevation = NA,
                       paleoData_variableName = "Abrupt shift + AR(1) (sigma = 0.2)",
                       paleoData_units = "unitless",
                       paleoData_values = list(as.matrix(df1$signal + df1$noise_lo)),
                       paleoData_TSid = "pseudoAS02")

row2 <- row1
row2$paleoData_values <- list(as.matrix(df1$signal + df1$noise_hi))
row2$paleoData_variableName = "Abrupt shift + AR(1) (sigma = 0.8)"
row2$paleoData_TSid = "pseudoAS08"

row3 <- row1
row3$paleoData_values <- list(as.matrix(df2$signal + df2$noise_lo))
row3$paleoData_variableName = "Gradual shift + AR(1) (sigma = 0.2)"
row3$paleoData_TSid = "pseudoGS02"


row4 <- row1
row4$paleoData_values <- list(as.matrix(df2$signal + df2$noise_hi))
row4$paleoData_variableName = "Gradual shift + AR(1) (sigma = 0.8)"
row4$paleoData_TSid = "pseudoGS08"



pseudoShifts <- dplyr::bind_rows(row1,row2,row3,row4)

pseudoShifts <- structure(pseudoShifts,class = c("lipd-ts-tibble",class(tibble::tibble())))


usethis::use_data(pseudoShifts, overwrite = TRUE)

