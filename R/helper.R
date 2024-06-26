#' Maximum number of consecutive values in a vector
#'
#' @param vec a vector
#' @param val the value to search for consecutive values
#' @param gte return indices for values greater than or equal to this number of consecutive values (default= NA, which is only the max). Note that if gte > the max, the index returned will be NA
#'
#' @concept https://stackoverflow.com/questions/37345853/find-the-maximum-and-mean-length-of-the-consecutive-true-arguments
#' @return the number of consecutive values. 0 means that value is not present in the vector
#' @export
#' @examples
#' maxConsecutive(c(1,2,7,5,7,7,7,3,4),val = 7)
maxConsecutive <- function(vec,val = TRUE,gte = NA){

  consec <- rle(vec)

  # max consecutive values
  wcv <- which(consec$values==val)

  if(length(wcv > 0)){
    mv <- max(consec$lengths[wcv])
  }else{
    mv <- 0
  }
  if(!is.finite(mv)){
    mv <- 0
  }

    if(is.na(gte)){gte <- mv}#get for max only


  if(mv > 0 & gte <= mv){
    io <- c() #initialize
    for(mvi in gte:mv){
      #get indices for max
      fi <- which(consec$lengths == mvi & consec$values == val)
      for(i in fi){
        io <- c(io,seq(cumsum(consec$lengths)[i]-(mvi-1),cumsum(consec$lengths)[i]))
      }
    }
  }else{
    io <- NA
  }

  return(list(max = mv,index = io))
}

createTibbleFromParameterString <- function(params){
  fullstring <- glue::glue("tibble::tibble({params})")
  tib <- eval(parse(text = fullstring))
  return(tib)

}

identicalVectorsList <- function(l){
  #same length?
  t1 <- all(purrr::map2_lgl(l,l[1],~ all(length(.x) == length(.y))))
  if(t1){
    t1 <- all(purrr::map2_lgl(l,l[1],~ all(.x==.y)))
  }
  return(t1)
}

identicalMatrixColumns <- function(l){
  l <- purrr::array_branch(l,2)
  return(all(purrr::map2_lgl(l,l[1],~ all(.x==.y))))
}

list2matrix <- function(l){
  #find max length
  ml <- max(purrr::map_dbl(l,length))

  #setup matrix
  mat <- matrix(NA, nrow = ml,ncol = length(l))

  for(i in 1:length(l)){
    tv <- l[[i]]
    mat[1:length(tv),i] <- tv
  }

  return(mat)

}

hasParameterEnsemble <- function(x){
  np <- purrr::map_dfr(x$event_detection[[1]]$parameters,createTibbleFromParameterString) |>
    purrr::map_dbl(\(x) length(unique(x)))

  return(any(np > 1))
}

whichParametersInEnsemble <- function(x){
  np <- purrr::map_dfr(x$event_detection[[1]]$parameters,createTibbleFromParameterString) |>
    purrr::map_dbl(\(x) length(unique(x)))

  tr <- names(np)[np > 1]
  return(paste(tr,collapse = ", "))
}
