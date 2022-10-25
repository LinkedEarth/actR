# this explains what this function is doing


# define signal characteristics
N <- 100
start <- 50
amp <- 1

# import relevant functions
source("helper_functions.R")
source("bayes_chgpt.R")

library(pracma)
library(dplyr)

# generate the signal 
signal <- make_transition(length=N, amp=amp, start=start) # later we will change this to allow for different types of jumps

# ROC parameters
M <- 200 # number of noise realizations 
tau <- 5  # tolerance for change detection
nprob = 10  # number of probability thresholds
prob_thresh <- linspace(0.1, 0.9, nprob)  # vector of probability thresholds 

sigma <- c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.5) # noise levels
#sigma <- c(1) 
nsigma <- length(sigma)
gamma <- 0.9 # noise autocorrelation 

# for debugging only
#M = 10

# define output arrays
true_pos <- array(0, dim=c(nprob,nsigma,M))
true_neg <- array(0, dim=c(nprob,nsigma,M))
false_pos <- array(0, dim=c(nprob,nsigma,M))
false_neg <- array(0, dim=c(nprob,nsigma,M))

# define cgangepoint detection parameters
X <- matrix(c(rep(1,N), seq(1,N)), ncol=2)  
m <- dim(X)[2]   
t = seq(1,N)
#d_min <- 5 
#k_0 <- c(0.001, 0.1)    
beta0 <- rep(0,m)
#v_0 <- 1 
k_max <- 2     
#num_samp <- 1000 

# define window for change detection
upper = start + tau
lower = start - tau

for (i in nsigma){ # can include a completion bar
  noise <- ar1.noise(n=N, ncol = M, gamma=gamma, sigma=sigma[i])  # define noise
  
  for (j in 1:M){
    noisy_signal <- signal + noise[,j]
    Y = noisy_signal$values
    year = noisy_signal$time
    #sig_0=var(Y)  # match prior variance to actual variance
    # apply changepoint detection
    result = Bayes_Chgpt(Y, X, t, k_max)
    chgpt_loc = result$chgpt_loc
    
    for (k in 1:nprob) {
      chgpts = chgpt_find(chgpt_loc = chgpt_loc)
      for (v in length(chgpts)){
        if (between(chgpts[v],lower, upper) == TRUE & chgpt_loc[[chgpts[v]]] >= prob_thresh){
          true_pos[i,j,k] = true_pos[i,j,k] + 1
        }
        else if (between(chgpts[v],lower, upper) == TRUE & chgpt_loc[[chgpts[v]]] <= prob_thresh) {
          false_neg[i,j,k] = false_neg[i,j,k] + 1
        }
        else if (between(chgpts[v],lower, upper) == FALSE & chgpt_loc[[chgpts[v]]] >= prob_thresh ) {
          false_pos[i,j,k] = false_pos[i,j,k] + 1
        }
        else if (between(chgpts[v],lower, upper) == FALSE & chgpt_loc[[chgpts[v]]] <= prob_thresh){
          true_neg[i,j,k] = true_neg[i,j,k] + 1
        }
      # see how many of these chgpts are between upper and lower, using the function "between": https://www.geeksforgeeks.org/check-if-a-numeric-value-falls-between-a-range-in-r-programming-between-function/
      
      # update true_pos[i,j,k] and false_pos[i,j,k]  (values can be more than 1; if 3 changepoints are found, you put 3 down)
    }
    
    
    # count True/False Positive/Negatives
    
    
    # if (result$chgpts <= upper & result$chgpts >= lower) {
    #   if (result$chgpt_loc[[result$chgpts]] >= prob_thresh) {  # do we need this?
    #     true_pos[i,j] = 1
    #     }
    #   else { 
    #     false_neg[i,j] = 1
    #     }
    # } else {
    #   if (result$chgpt_loc[[result$chgpts]] >= prob_thresh) {
    #     false_pos[i,j] = 1
    #     }
    #   else {
    #     true_neg[i,j] = 1
    #     }
    # } 
  }
}

# next steps:
# - generalize your code so that it can handle more than one changepoint 
# - debug  and loop over 8 values of prob_thresh 
# - obtain the true positive and false positive rates for each value of prob_thresh
# - chart ROC curve

