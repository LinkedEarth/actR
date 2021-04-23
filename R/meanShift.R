meanShiftCore = function(age, vals, maxDiff, alpha = 0.05, plotOpt = F, gaussOpt = F, figDir = '', datNam = '', varNam = '') {
  ## Written by Hannah Kolus, 09/04/2018
  ## Searches for mean shifts in the specified record, using the changepoint package function cpt.mean().
  ## In addition, tests for significant differences in means from the results of cpt.mean(), with the
  ## potential to reject some changepoints as insignificant.
  if (gaussOpt) {
    vals = gaussianize(vals)
  }

  # interpolate the data to the min resolution of the record
  res = abs(min(diff(age)))
  minSeg = ceiling(maxDiff / res)
  f = approxfun(age,vals)
  X = seq(min(age), max(age), by = res)
  Y = f(X)

  # run the mean shift code
  results = cpt.mean(Y,penalty="SIC",method="PELT",minseglen = minSeg)
  resInds = cpts(results)

  # Delete change points detected at very end of record
  resInds[resInds == length(X)] = NA
  resInds = resInds[!is.na(resInds)]
  breaks = X[resInds]

  # --- Test significance of means delineated by breaks --- #
  segInds = c(1, resInds, length(X))
  segBounds = c(X[1], breaks, X[length(X)])
  pvals = rep(NA, length(breaks))
  reduceBreaks = ''
  for (it in 1:length(breaks)) {

    # Using original data for the t-test
    x = vals[which(age >= segBounds[it] & age < segBounds[it+1])]
    y = vals[which(age >= segBounds[it+1] & age <= segBounds[it+2])]

    # If vector is too short for t-test, record placeholder p-val (1) and skip to next change point
    if (length(x) < 2 || length(y) < 2) {
      print(paste('Insignificant difference in means: CHANGE POINT AGE: ', breaks[it]))
      reduceBreaks = 'REDUCED'
      pvals[it] = 1
      next()
    }

    # Perform t-test
    testRes = try(t.test(x, y, alternative = "two.sided", conf.level = 1 - alpha, var.equal = FALSE), silent = T)

    # record pvals
    if (class(testRes) == "try-error") {
      pvals[it] = 1
      reduceBreaks = 'REDUCED'
      print(paste('Try error'))
    } else {
      pvals[it] = testRes$p.value
    }
  } # end loop thru breaks

  # Store significant breaks in TS_MS
  if (!any(pvals < alpha)) {
    sig_brks = NA
    brk_dirs = NA
  } else {
    sig_brks = breaks[pvals < alpha]

    # Calculate mean over significant segments
    sigBounds = c(X[1], sig_brks, X[length(X)])
    sigBounds = sigBounds[!is.na(sigBounds)]
    origMean = rep(0, length(sigBounds) - 1)
    for (it in 1:length(origMean)) {
      origMean[it] = mean(vals[which(age >= sigBounds[it] & age <= sigBounds[it+1])])
    }

    brk_dirs = ifelse(diff(origMean) < 0, 1, -1)
  }

  ## ------------------------------ PLOTTING ------------------------------ ##
  if (plotOpt) {

    # Calculate mean over significant segments
    sigBounds = c(X[1], sig_brks, X[length(X)])
    sigBounds = sigBounds[!is.na(sigBounds)]
    interpMean = rep(0, length(sigBounds) - 1)
    origMean = rep(0, length(sigBounds) - 1)
    for (it in 1:length(origMean)) {
      origMean[it] = mean(vals[which(age >= sigBounds[it] & age <= sigBounds[it+1])])
      interpMean[it] = mean(Y[which(X >= sigBounds[it] & X <= sigBounds[it+1])])
    }

    # Remove troublesome characters for file saving
    varNam = gsub('/','_',varNam)
    varNam = gsub(',','_',varNam)
    varNam = gsub('%','_',varNam)
    datNam = gsub('/','_',datNam)
    datNam = gsub(',','_',datNam)

    # plot each original and interpolated series with change points and segment means
    pdf(file.path(figDir, paste0(datNam, '_', varNam, '_', reduceBreaks, '.pdf')))
    par(mfrow = c(2,1))
    plot(X, Y, type = 'o', xlab="Age (BP)", ylab = "Values")
    abline(v = sig_brks, col="blue", lwd = 2)
    for (p in 1:(length(sigBounds) - 1)) {
      segments(sigBounds[p], interpMean[p], sigBounds[p+1], interpMean[p], col = 'red', lwd = 2)
    }
    text(sig_brks[1], ifelse(max(vals) < 0, max(vals)*1.1, max(vals)/1.1), ifelse(brk_dirs[1] == 1, 'pos', 'neg'))
    title(paste0(datNam, ':\nInterpolated Data'))

    plot(age, vals, type = 'o', xlab="Age (BP)", ylab = "Values")
    abline(v = sig_brks, col="blue", lwd = 2)
    for (p in 1:(length(sigBounds) - 1)) {
      segments(sigBounds[p], origMean[p], sigBounds[p+1], origMean[p], col = 'red', lwd = 2)
    }
    text(sig_brks[1], ifelse(max(vals) < 0, max(vals)*1.1, max(vals)/1.1), ifelse(brk_dirs[1] == 1, 'pos', 'neg'))
    title('Original Data')
    dev.off()
  } # end plotting

  return(list(sig_brks = sig_brks, brk_dirs = brk_dirs))

} # end function
