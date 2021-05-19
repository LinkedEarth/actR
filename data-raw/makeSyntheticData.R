library(arima)
library(ggpubr)
library(ggthemes)
library(latex2exp)


midpoint.jump <-function(lngth,hght){
  "Generate a Heaviside jump of height `hght` at the midpoint of a series of length `lngth`  "
  jump = rep(1.,lngth)
  jump[1:round(lngth/2,0)] = 0.
  return(jump)
}

linear.trend <- function(lngth,hght){
  return(seq(lngth)/lngth*hght)
}

ar1.noise <- function(lngth,g,sigma){
  ar1=sigma/sqrt(1-g^2)*arima.sim(model=list(g,0,0),n=lngth)
  return(ar1)
}


# global variables
theme_set(theme_hc(style = "darkunica"))
lngth <- 400
time <- 1.0*seq(lngth)

# case 1: noisy jump
signal1 <- midpoint.jump(lngth,1)
noise_lo <- ar1.noise(lngth,0.8,0.2)
noise_hi <- ar1.noise(lngth,0.8,0.8)

df1 <- data.frame(signal=signal1,noise_hi,noise_lo)
c1a <- ggplot(data=df1) + geom_line(aes(x=time,y=signal),color='white') +
  geom_line(aes(x=time,y=signal+noise_lo),color='orange') +
  ylab("y(t)") + xlab("t") + ylim(-1,2) +
  ggtitle(TeX(r'(Jump + AR(1), $\sigma=0.2$)'))
show(c1a)

c1b <- ggplot(data=df1) + geom_line(aes(x=time,y=signal),color='white') +
  geom_line(aes(x=time,y=signal+noise_hi),color='orange') +
  ggtitle(TeX(r'(Jump + AR(1), $\sigma=0.8$)')) +
  ylab('y(t)') + xlab("t") + ylim(-1,2)


# case 2: noisy jump with trend
signal2 <- signal1 + linear.trend(lngth,1)
df2 <- data.frame(signal=signal2,noise_hi,noise_lo)
c2a <- ggplot(data=df2) + geom_line(aes(x=time,y=signal),color='white') +
  geom_line(aes(x=time,y=signal+noise_lo),color='orange') +
  ylab("y(t)") + xlab("t") + ylim(-1,2) +
  ggtitle(TeX(r'(Jump + trend + AR(1), $\sigma=0.2$)'))
c2b <- ggplot(data=df2) + geom_line(aes(x=time,y=signal),color='white') +
  geom_line(aes(x=time,y=signal+noise_hi),color='orange') +
  ylab("y(t)") + xlab("t") + ylim(-1,2) +
  ggtitle(TeX(r'(Jump + trend + AR(1), $\sigma=0.8$)'))

plot22 <- ggarrange(c1a,c1b,c2a,c2b)
show(plot22)

# case 3: noisy jump with age errors

