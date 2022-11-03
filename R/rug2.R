#rug2
# Binary Segmentation
BS = function(Y.fun, beginning=1, end=length(Y.fun), linear = FALSE)
{
  chgpts = c()
  found = 0

  N = length(Y.fun)
  x = (1:N)/N
  criteria = 3.37*sqrt(x*(1-x))    # detection criteria
  # Use 3.83 for alpha = 0.01, 3.37 for alpha = 0.05, and 3.13 for alpha = 0.10
  # Alternate criteria:
  #criteria = rep(1.3582, N)        # criteria = 1.629 for alpha = 0.01; criteria = 1.3582 for alpha = 0.05
  # If this criteria is used, you will need to adjust the which.max command below by eliminating the x*(1-x) part, which is currently commented out

  X = matrix(c(rep(1,N), 1:N), ncol = 2, byrow = FALSE)

  c=1

  if(linear)
  { c = 2  }

  model = lm(Y.fun[1:N]~X[,c])

  hat.beta = model$coefficients    # not actually needed
  resid = model$residuals
  sig.hat = summary(model)$sigma

  W.n = rep(0,N)
  for (n in (c+1):(N-1))
  {
    W.n[n] = sum(resid[1:n])/sig.hat/sqrt(N)
    if (abs(W.n[n])> criteria[n])
    {
      found = 1    # If W.n exceeds the threshold criteria, then a change point exists
    }
  }

  if(found ==1)
  {
    #new.chgpt = which.max(abs(W.n)/sqrt(x*(1-x)))
    # The above line of code works for a constant model, but not necessarily for a linear model
    # So we'll brute force the location of the change point instead (which is slower, but accurate)
    new.chgpt = OLS(Y.fun, c(1, N), linear)
    chgpts = c(chgpts, beginning +new.chgpt-1 )  # adjusts for our continually cutting the data set
    # Now split the data and look for more change points!
    chgpts = c(chgpts, BS(Y.fun[1:new.chgpt], beginning, beginning +new.chgpt-1, linear))
    chgpts = c(chgpts, BS(Y.fun[(new.chgpt+1):N], beginning + new.chgpt, end, linear))
  }
  return(chgpts)
}

# A helper function... also used in our implementation of the WBS algortihm below
OLS = function(Y, FM, linear)
{
  # Finds the actual position of the change point... which may be different than the max of W.n

  start = FM[1]; finish = FM[2]

  N = length(Y)
  X = matrix(c(rep(1,N), 1:N), ncol = 2, byrow = FALSE)

  c=1

  if(linear)
  { c = 2  }

  temp = rep(Inf, N)
  for (i in (start+c):(finish-c))
  {
    M1 = lm(Y[start:i]~X[start:i,c])
    M2 = lm(Y[(i+1):finish] ~ X[(i+1):finish, c])
    temp[i] = sum(M1$residuals^2) + sum(M2$residuals^2)

  }
  C = which.min(temp)  # i.e. which position will minimized the squared error if selected as the change point

  return(C)
}


# Wild Binary Segmentation
# To use the Wild Binary Segmentation algorithm, install the “wbs” package. The input to this function is a vector, Y, which contains the data. The function returns a number of different variables, among which is a list of the change points that were detected.

library(wbs)


# Modified Version of Wild Binary Segmentation
# The Wild Binary Segmentation algorithm included in the “wbs” package can only be used to detect changes in the mean. Since we are also interestd in looking for change points in linear functions, a modified version of the original “wbs” function was created (by us) that can also be applied to linear functions. Again, the input to wbs.linear() function is a vector, Y, which contains the data. By default, the parameter linear is FALSE, so the model will try to fit a constant function, making it equivalent to the original “wbs” function (in theory). Changing linear to TRUE will cause the wbs.linear() function to fit linear functions instead. The output to the function is a list of the change points that were detected, which is different than the wbs() function.

# Wild Binary Segmentation - Modified to handle linear functions
# The initial function creates the random segments and finds max(W.n) for each.
# The functions that follow find the change points and calculate the value of W.n, respectively.
wbs.linear = function(Y, linear = FALSE)
{

  N = length(Y)

  X = matrix(c(rep(1,N), 1:N), ncol = 2, byrow = FALSE)

  c=1

  if(linear)
  { c = 2  }

  num.samp = 5*N       # changing to zero makes this BS
  index = sample(1:N, num.samp*2, replace = TRUE)
  index = c(index, 1, N)     # We'll also include the full segment as one of our intervals
  FM = matrix(index, ncol = 2, byrow = TRUE)    # This creates the localized intervals

  test = rep(0, num.samp+1)


  for (i in 1:num.samp)
  {
    if (FM[i,1]>FM[i,2])    # we need to make sure the starting point is before the ending point
    {                       # otherwise swap the order of the FM's
      temp = FM[i,1]
      FM[i,1]=FM[i,2]; FM[i,2] = temp
    }
    if (FM[i,1]-FM[i,2]>c)   # i.e. enough data points to fit a model
    {
      W.n = find_CUSUM(Y, FM[i,], linear)    #find CUSUM statistic (i.e. W.n) for the interval
      test[i] = max(abs(W.n))        # Save the max value of this function to test for significance
    }
  }

  # Let's also consider the full segment
  W.n = find_CUSUM(Y, c(1,N), linear)
  test[num.samp+1] = max(abs(W.n))


  #Now search for change points
  chgpts = c()
  C = find_chgpts(Y, FM, test, 1, N, chgpts, linear)
  return(C)
}

# Helper function which detects change points based on the value of W.n for each random interval
find_chgpts = function(Y, FM, TEST, start, finish, chgpts, linear)
{
  #criteria = 1.629    #alpha = 0.01
  #criteria = 1.481    #alpha = 0.025
  criteria = 1.3582    #alpha = 0.05
  test = TEST            # so we don't overwrite the master copy
  num.samp = length(TEST) -1

  r = dim(FM)[1]

  for (i in 1:r)
  {
    if (FM[i,1] < start | FM[i,2] > finish)
    {  test[i] = 0  }
    # This ensures we are only checking intervals within our current bounds
  }

  M = max(test); pos = which.max(test)

  if (M>criteria)    # i.e. if a change point exists
  {
    C = OLS(Y, FM[pos,], linear)   #Find the actual location of the change point

    chgpts = c(chgpts, C)

    # This is the binary segmentation part
    W_n = find_CUSUM(Y, c(start, C), linear)   # We have a new "complete" segment, the left half of the data
    TEST[num.samp+1] = max(abs(W_n))
    FM[(num.samp+1),] = c(start, C)
    chgpts = find_chgpts(Y, FM, TEST, start, C, chgpts, linear)

    W_n = find_CUSUM(Y, c(C+1, finish), linear) # A new "complete" segment, the right half of the data
    TEST[num.samp+1] = max(abs(W_n))
    FM[(num.samp+1),] = c(C+1, finish)
    chgpts = find_chgpts(Y, FM, TEST, C+1, finish, chgpts, linear)
  }

  return(chgpts)
}

# A helper fucntion which calculates the value of W.n for an arbitrary sub-interval of the data set.
find_CUSUM = function(Y, FM, linear)
{
  N = FM[2]-FM[1] + 1
  W_n = rep(0,N)
  X = matrix(c(rep(1,N), 1:N), ncol = 2, byrow = FALSE)

  c=1

  if(linear)
  { c = 2  }

  if (N > c)
  {
    model = lm(Y[FM[1]:FM[2]]~X[,c])

    hat.beta = model$coefficients    # not actually needed
    resid = model$residuals
    sig.hat = summary(model)$sigma

    W.n = rep(0,N)
    for (n in (c+1):(N-c))
    {
      W.n[n] = sum(resid[1:n])/sig.hat/sqrt(N)
    }

  }
  return(W.n)
}


# The Bayesian Change Point Algorithm
# Finally, we give the functions used in the Bayesian Change Point algorithm. The first two are helper functions, and the last is the main function.

# Step 1: Calculate the probability of the data for every possible substring
probability2 <- function(i,j, Y, X, d_min, k_0, sig_0, beta0, v_0)
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
partition_fn2 <- function(Py, k_max, N)
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


# The main function
Bayes_Chgpt <- function(Y, X, t=seq(1,length(Y)), k_max=10, d_min=10, k_0 = c(0.1, 0.1), v_0=1, sig_0=var(Y), num_samp = 500)
{
  mean.Y = mean(Y)     # For numerical stability
  Y = Y - mean.Y
  N <<- length(Y)    # Total number of observations (a global variable)
  m <<- dim(X)[2]    # Number of predictor variables (a global variable)

  I = diag(m)
  beta0 <<- rep(0,m) # Mean of multivariate normal prior on regression coefficients


  #***** Step 1: Calculate the Probability of the Data for Every Possible Substring

  Py = matrix(-Inf, nrow=N, ncol=N)

  for (i in 1:N)
  {
    for (j in i:N)
    {
      Py[i,j] = probability2(i,j, Y, X, d_min, k_0, sig_0, beta0, v_0)
    }
  }

  #*********** Step 2: Forward Recursion ************************************

  P = partition_fn2(Py, k_max, N)


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

  for (i in 1:num_samp)
  {
    #***** (3.1) Sample a Number of Change Points ***************
    num_chgpts = sample(0:k_max, size=1, prob=k)

    back = N
    if (num_chgpts>1)
    {
      #*****(3.2) Sample the Locations of the Change Points *****
      for(kk in num_chgpts:2)  # Start at the end of the time series and work backwards
      {
        temp = unlist(lapply(1:(back-1), function(x) P[kk-1,x]+Py[x+1,back]))
        M_temp = max(temp); temp = temp - M_temp   # Used to avoid underflow

        total = log(sum(exp(temp)))
        # Equation (4) - Marginalize over all possible placements of the change point.
        temp = exp(temp-total)        # Normalize the vector
        changepoint = sample(1:(back-1), size = 1, prob=temp)  # Sample the location of the change point
        chgpt_loc[changepoint] = chgpt_loc[changepoint]+1   # Keep track of change point locations
        samp_holder[i,kk] = changepoint

        #********* (3.3) Sample the Regression Parameters ********
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
    if(num_chgpts >0)
    {
      # The Final Changepoint
      #*****(3.2) Sample the Locations of the Change Points *****
      kk=1

      temp = unlist(lapply(1:(back-1), function(x) Py[1,x]+Py[x+1,back]))
      M_temp = max(temp); temp = temp - M_temp   # Used to avoid underflow

      total = log(sum(exp(temp)))
      # Equation (4) - Marginalize over all possible placements of the change point.
      temp = exp(temp-total)        # Normalize the vector
      changepoint = sample(1:(back-1), size=1, prob=temp)  # Sample the location of the change point
      chgpt_loc[changepoint] = chgpt_loc[changepoint]+1   # Keep track of change point locations
      samp_holder[i,kk] = changepoint

      #********* (3.3) Sample the Regression Parameters ********
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

    }
    else           # zero change points, so a single homogeneous segment
    {
      #********* (3.3) Sample the Regression Parameters ********
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
  chgpt_loc=chgpt_loc/num_samp  # Posterior probability of a change point


  #*************** Plot the Results **********************

  temp = which(diff(sign(diff(chgpt_loc)))==-2)+1
  chgpts = NULL
  for (i in temp)
  {
    if(chgpt_loc[i] > 0.1)
    {
      chgpts = c(chgpts,i)
    }
  }
  model = unlist(lapply(1:N, function(x) X[x,]%*%BETA[,x]))
  size = abs(max(Y) - min(Y))
  Y = Y+mean.Y
  model = model+mean.Y

  par(mar=c(5, 4, 4, 5))
  plot(t,Y, type = 'l', col ='blue', xlab = NA, ylim = c(min(Y)-0.33*size, max(Y)))
  lines(t,model, type = 'l', col = 'green', xlab = NA)
  par(new=TRUE)
  par(mar=c(5, 4, 4, 5))
  plot(t,chgpt_loc, type = 'l', axes = FALSE, col = 'red', xlab = NA, ylab = NA, ylim = c(0,3))
  axis(side=4, at=seq(0,1,0.2), las = 2, tck = 0.025, cex.axis = 0.75, line = 0, col = 'red', col.axis = 'red')
  axis(side=1,at=chgpts,las=2,cex.axis = 0.75, line = 0,col.axis = 'red',col="red")
  mtext("Probability of a Change Point", side = 4, line = 2, adj = 0)
  title("Change Points in the Data Set", xlab = "Time")


  return( list(k=k, model=model, BETA=BETA, chgpt_loc=chgpt_loc, samp_holder=samp_holder))
}  # end of Bayes_Chgpt function
