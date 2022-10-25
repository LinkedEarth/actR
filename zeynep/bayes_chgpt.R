# Step 1: Calculate the probability of the data for every possible substring
probability <- function(i,j, Y, X, d_min, k_0, sig_0, beta0)
{
  Py = -Inf
  if (j-i>=d_min)
  {
    # Calculate beta_hat
    J = t(X[i:j,])%*%X[i:j,] + diag(k_0, nrow = length(k_0))
    # A diagonal matrix with the vector k_0 as the entries
    beta_hat = solve(J)%*%(k_0*beta0 + t(X[i:j,])%*%Y[i:j])
    # Should be matrix multiplication, but R does vector multiplication by components
    
    # Calculate v_n and s_n
    a = j-i+1    # "a" is the number of data points, not the distance between them, which becomes
    # relevant if the data points are not equally spaced
    v_n = v_0 + a  # a is the number of data points in the segment
    s_n = v_0*sig_0 + sum(k_0*(beta0-beta_hat)^2) + sum((Y[i:j]-X[i:j,]%*%beta_hat)^2)
    #Again, you get away with this because of how R does vector multiplication
    
    # Calculate the probability density of the data, stored in log form
    Py = v_0*log(sig_0*v_0/2)/2 -v_n*log(s_n/2)/2 +log(prod(k_0))/2 +
      lgamma(v_n/2) - lgamma(v_0/2) -a*log(2*pi)/2 - log(det(J))/2
    
  }
  return(Py)
}

# Step 2: Forward Recursion
partition_fn <- function(Py, k_max, N)
{
  # The Forward Recursion step of the Bayesian Change Point Algorithm.
  # P(k,j) is the probability that the first j data points contain k change points [Equation 4]
  
  P = matrix(-Inf, nrow=k_max, ncol=N)        # -Inf b/c starts in log form
  
  #k is the number of change points
  k=1            # First row is different from the rest, as you add together two homogeneous segments
  
  for (j in (k+1):N)     # Note: Several of these terms will be -INF, due to d_min parameter
  {
    temp = unlist(lapply(1:(j-1), function(x) Py[1,x]+Py[x+1,j]))
    
    #NOTE: TO AVOID UNDERFLOW, USE:
    M_temp = max(temp)
    if (M_temp>-Inf)
    {
      temp = temp - M_temp
      P[k,j]=log(sum(exp(temp))) +M_temp
      # Equation (4) - Marginalize over all possible placements of the change point.
    }
    else
    {
      P[k,j] = -Inf
    }
    
  }
  
  for (k in 2:k_max)
  {
    for (j in (k+1):N)  # Note: Several of these terms will be -INF as well
    {
      
      temp = unlist(lapply(1:(j-1), function(x) P[k-1,x]+Py[x+1,j]))
      
      # NOTE: TO AVOID UNDERFLOW, USE:
      M_temp = max(temp)
      if (M_temp>-Inf)
      {
        temp = temp - M_temp
        P[k,j]=log(sum(exp(temp))) +M_temp
        # Equation (4) - Marginalize over all possible placements of the change point.
      }
      else
      {
        P[k,j] = -Inf;
      }
      
    }
  }
  return(P)
}         # of function partition_fn

#' @export
#' @name bayes_chgpt_compute
#'
#' Detects changepoints in series Y based on design matrix X
#' @param Y a vector of length N
#' @param X design matrix,  a matrix of dimensions N x m, where m is the number of predictors
#' @param t time axis common to Y and X (is it needed at all?)
#' @param k_max maximum number of change points allowed 
#' @param d_min minimum distance between adjacent change points
#' @param k_0 hyper parameter for the prior on the regression coefficients
#' @param v_0 prior location parameter of the Scaled − Inverse $\chi^2$ distribution for \sigma
#' @param sig_0 prior scale parameter of the Scaled − Inverse $\chi^2$ distribution for \sigma
#' @param num_samp number of sampled solutions
#' @return a list containing k, \beta, chgpt_loc, samp_holder, model, size, year, and Y

# The main function
bayes_chgpt_compute <- function(Y, X, t=seq(1,length(Y)), k_max=10, d_min=5, k_0 = c(0.001, 0.1), v_0=1, sig_0=std(Y), num_samp = 1000)
{
  mean.Y <- mean(Y)     # For numerical stability
  Y <- Y - mean.Y
  N <- length(Y)    # Total number of observations (a global variable)
  m <- dim(X)[2]    # Number of predictor variables (a global variable)
  
  I <-  diag(m)
  beta0 <- rep(0,m) # Mean of multivariate normal prior on regression coefficients
  
  
  #***** Step 1: Calculate the Probability of the Data for Every Possible Substring
  
  Py = matrix(-Inf, nrow=N, ncol=N)
  
  for (i in 1:N)
  {
    for (j in i:N)
    {
      Py[i,j] = probability(i,j, Y, X, d_min, k_0, sig_0, beta0)
    }
  }
  
  #*********** Step 2: Forward Recursion ************************************
  
  P = partition_fn(Py, k_max, N)
  
  
  #*************** Step 3: Stochastic Backtrace *****************************
  k = P[,N]
  
  for (i in 1:k_max)
  {
    if (N - (i+1)*d_min+i >= i)   # i.e. if a plausible solution exists
    {
      k[i] = k[i] + log(0.5) - log(k_max) - log(choose(N-(i+1)*d_min + i, i))
      # Above item is N_k
    }
  }
  k = c(Py[1,N]+log(0.5), k)      # Adding in the possibility of zero change points
  
  # NOTE: TO AVOID UNDERFLOW, USE:
  temp = max(k)
  k = k - temp
  #
  
  total_k = log(sum(exp(k)))      # Adds logariths and stores answer in log form
  k = exp(k-total_k)              # Normalize the vector - Used in Equations 4 and 7
  
  
  # Variables used in the sampling procedure:
  samp_holder = matrix(0, nrow = num_samp, ncol = k_max)
  chgpt_loc = rep(0,N)            # The number of times a change point is identified at each data point
  BETA = matrix(0, nrow = m, ncol = N)
  
  for (i in 1:num_samp){
    #  (3.1) Sample a Number of Change Points
    num_chgpts = sample(0:k_max, size=1, prob=k)
    
    back = N
    if (num_chgpts>1){
      # (3.2) Sample the Locations of the Change Points
      for(kk in num_chgpts:2){  # Start at the end of the time series and work backwards
        temp = unlist(lapply(1:(back-1), function(x) P[kk-1,x]+Py[x+1,back]))
        M_temp = max(temp); temp = temp - M_temp   # Used to avoid underflow
        
        total = log(sum(exp(temp)))
        # Equation (4) - Marginalize over all possible placements of the change point.
        temp = exp(temp-total)        # Normalize the vector
        changepoint = sample(1:(back-1), size = 1, prob=temp)  # Sample the location of the change point
        chgpt_loc[changepoint] = chgpt_loc[changepoint]+1   # Keep track of change point locations
        samp_holder[i,kk] = changepoint
        
        # (3.3) Sample the Regression Parameters
        # Regression Coefficients (Beta)
        # Calculate beta_hat
        J = solve(t(X[(changepoint+1):back,])%*%X[(changepoint+1):back,] + k_0*I)
        beta_hat = J%*%(k_0*beta0 + t(X[(changepoint+1):back,])%*%Y[(changepoint+1):back])
        
        v_n = v_0 + back-changepoint  # the number of data points in the segment
        s_n = v_0*sig_0 + k_0*sum((beta0-beta_hat)^2) + sum((Y[(changepoint+1):back]-X[(changepoint+1):back,]%*%beta_hat)^2)
        samp_sigma = 1/(rchisq(1, v_n))*s_n
        Z = rnorm(m)
        betas = Z%*%chol(samp_sigma*J)+t(beta_hat)
        
        BETA[,(changepoint+1):back] = BETA[,(changepoint+1):back] + rep(betas, (back-changepoint))
        
        back = changepoint
      }
    }
    if(num_chgpts >0){
      # The Final Changepoint
      # (3.2) Sample the Locations of the Change Points
      kk=1
      temp = unlist(lapply(1:(back-1), function(x) Py[1,x]+Py[x+1,back]))
      M_temp = max(temp); temp = temp - M_temp   # Used to avoid underflow
      
      total = log(sum(exp(temp)))
      # Equation (4) - Marginalize over all possible placements of the change point.
      temp = exp(temp-total)        # Normalize the vector
      changepoint = sample(1:(back-1), size=1, prob=temp)  # Sample the location of the change point
      chgpt_loc[changepoint] = chgpt_loc[changepoint]+1   # Keep track of change point locations
      samp_holder[i,kk] = changepoint
      
      # (3.3) Sample the Regression Parameters
      # Regression Coefficients (Beta)
      # Calculate beta_hat
      J = solve(t(X[(changepoint+1):back,])%*%X[(changepoint+1):back,] + k_0*I)
      beta_hat = J%*%(k_0*beta0 + t(X[(changepoint+1):back,])%*%Y[(changepoint+1):back])
      
      v_n = v_0 + back-changepoint  # the number of data points in the segment
      s_n = v_0*sig_0 + k_0*sum((beta0-beta_hat)^2) + sum((Y[(changepoint+1):back]-X[(changepoint+1):back,]%*%beta_hat)^2)
      samp_sigma = 1/(rchisq(1, v_n))*s_n
      
      Z = rnorm(m)
      betas = Z%*%chol(samp_sigma*J)+t(beta_hat)
      
      BETA[,(changepoint+1):back] = BETA[,(changepoint+1):back] + rep(betas, (back-changepoint))
      
      # The final sub-interval
      J = solve(t(X[1:changepoint,])%*%X[1:changepoint,] + k_0*I)
      beta_hat = J%*%(k_0*beta0 + t(X[1:changepoint,])%*%Y[1:changepoint])
      
      v_n = v_0 + changepoint  # the number of data points in the segment
      s_n = v_0*sig_0 + k_0*sum((beta0-beta_hat)^2) + sum((Y[1:changepoint]-X[1:changepoint,]%*%beta_hat)^2)
      samp_sigma = 1/(rchisq(1, v_n))*s_n
      
      Z = rnorm(m)
      betas = Z%*%chol(samp_sigma*J)+t(beta_hat)
      
      BETA[,1:changepoint] = BETA[,1:changepoint] + rep(betas, changepoint)
    } else { # zero change points, so a single homogeneous segment
      # (3.3) Sample the Regression Parameters
      # Regression Coefficients (Beta)
      J = solve(t(X)%*%X + k_0*I)
      beta_hat = J%*%(k_0*beta0 + t(X)%*%Y)
      v_n = v_0 + N  # the number of data points in the segment
      s_n = v_0*sig_0 + k_0*sum((beta0-beta_hat)^2) + sum((Y[1:N]-X[1:N,]%*%beta_hat)^2)
      samp_sigma = 1/(rchisq(1, v_n))*s_n
      
      Z = rnorm(m)
      betas = Z%*%chol(samp_sigma*J)+t(beta_hat)
      BETA = BETA + rep(beta_hat, N)
    }
  }
  
  BETA = BETA/num_samp      # Average regression coefficient at each data point
  chgpt_prob=chgpt_loc/num_samp # Posterior probability of a change point
  
  ### (4) Export output
  
  model = unlist(lapply(1:N, function(x) X[x,]%*%BETA[,x]))
  size = abs(max(Y) - min(Y))
  # restore mean
  Y = Y+mean.Y
  model = model+mean.Y
  
  return(list(k=k, BETA=BETA, chgpt_prob=chgpt_prob, samp_holder=samp_holder, model=model, size=size, year=t, Y=Y))
}

#' bayes_chgpt_plot
#'
#' Plots the results of the bayesian change point algorithim
#' @param year y-axis 
#' @param Y x-axis values of data
#' @param model model fit for noisy data
#' @param chgpt_prob posterior probability of a changepoint
#' @param size scaling parameter for y-axis  
#' @param p1 ggplot object that plots signal + noise, signal, and model
#' @param p2 ggplot object that plots the changepoint probability over time
#' @return a list containing p1 and p2 

bayes_chgpt_plot <- function(year, Y, model, chgpt_prob, size){
  dataframe = data.frame(year,Y, model, chgpt_prob)
  p1 <- ggplot(dataframe) + geom_line(mapping=aes(x=year,y=Y,colour='Signal + Noise')) + ylim(c(min(Y)-0.33*size, max(Y)))  + geom_line(mapping=aes(x=year,y=model,colour='Model')) + scale_color_manual(name = " ", breaks = c('Signal + Noise',  'Model', 'Signal'), values=c('Signal + Noise'='blue','Model'='darkgoldenrod2', 'Signal'='black')) + xlab("") + theme(legend.position="top",axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.margin=unit(c(1,1,0,1), "cm"))
  p2 <- ggplot(dataframe) + geom_line(mapping=aes(x=year,y=chgpt_prob,colour='red')) + ylab("P(Change Point)") + xlab("Year") + theme(legend.position="none", plot.margin=unit(c(0,1,1,1), "cm"))
  return((list(p1=p1,p2=p2)))
}

#' chgpt_find
#' 
#' Finds the index of specific changepoints                        
#' @param chgpt_prob posterior probability of a changepoint
#' @param chgpts index of where changepoint occurs
#' @return returns chgpts 

chgpt_find <- function (chgpt_prob){
  temp = which(diff(sign(diff(chgpt_prob)))==-2)+1
  chgpts = NULL
  for (i in temp){
    if(chgpt_prob[i] > 0.1)
    {
      chgpts = c(chgpts,i)
    }
  }
  return(chgpts)
}

#' bayes_chgpt_roc
#'
#' Computes the Receiver Operator Characteristic for the function bayes_chgpt_compute()
#' @param signal a vector of length N
#' @param noise a matrix of dimensions N x M, where M is the number of realizations
#' @param prob_thresh a vector of threshold probabilities, length P
#' @param loc location of the changepoint to identify
#' @param tol tolerance in identifying the changepoint
#' @param k_max maximum number of change points allowed 
#' @param true_pos matrix containing instances of true positives detected
#' @param false_pos matrix containing instances of false positives detected
#' @param true_neg matrix containing instances of true negatives 
#' @param false_neg matrix containing instances of false negatives
#' @return a list containing True_pos, false_pos, true_neg, and false_neg

bayes_chgpt_roc <- function (signal, noise, prob_thresh, loc = start, tol = 3, k_max = 2){
  
  # define array sizes
  N <- dim(noise)[1]
  M <- dim(noise)[2]
  P <- length(prob_thresh)
  
  # define output arrays
  true_pos <- matrix(0, nrow=P, ncol=M)
  false_pos <- matrix(0, nrow=P, ncol=M)
  true_neg <- matrix(0, nrow=P, ncol=M)
  false_neg <- matrix(0, nrow=P, ncol=M)
  
  # define changepoint detection parameters
  t <- seq(1, N)
  X <- matrix(c(rep(1, N), t), ncol = 2)
  m <- dim(X)[2]
  beta0 <- rep(0, m)
  
  d_min <- 5      # minimum distance between change points
  k_0 <- c(0.001, 0.1)     # Hyper parameter for the prior on the regression coefficients
  v_0 <- 1; sig_0=var(Y)   # You can also use: sig_0 = var(Y)  
  # Hyperparameters for the scaled-inverse chi-square distribution
  k_max <- 10        # Maximum number of change points
  num_samp <- 1000   # number of sampled solutions
  
  # define window for change detection
  upper = loc + tol
  lower = loc - tol
  
  for (j in 1:M) { #loop over noise realizations
    Y <- signal$values + noise[, j] # add noise
    sig_0 = std(Y)  # match prior standard deviation to estimated standard deviation
    # apply changepoint detection
    result = bayes_chgpt_compute(Y, X, t, k_max)
    chgpt_prob = result$chgpt_prob
    chgpts = chgpt_find(chgpt_prob = chgpt_prob)
    
    if (!is.null(chgpts)) {
      for (k in 1:P) {
        #loop over probability thresholds
        for (v in length(chgpts)) {
          # loop over changepoints
          if (between(chgpts[v], lower, upper) &
              chgpt_prob[[chgpts[v]]] >= prob_thresh[k]) {
            true_pos[k, j] = true_pos[k, j] + 1
          } else if (between(chgpts[v], lower, upper) &
                       chgpt_prob[[chgpts[v]]] < prob_thresh[k]){
            false_neg[k,j] = false_neg[k, j] + 1
          } else if (chgpt_prob[[chgpts[v]]] >= prob_thresh[k]){
            false_pos[k, j] = false_pos[k, j] + 1
          } else if (chgpt_prob[[chgpts[v]]] < prob_thresh[k]) {
            true_neg[k,j] = true_neg[k,j] + 1
          } 
        } # end v loop
      } # end k loop
    } # end non-empty check
  } # end j loop
  return(list(true_pos=true_pos,false_neg=false_neg,true_neg=true_neg,false_pos=false_pos))
}
  
