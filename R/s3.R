new_excursion <- function(x = tibble::tibble()) {
  stopifnot(tibble::is_tibble(x))

  structure(x,class = c("excursion",class(tibble::tibble())))
}

new_excursionCore <- function(x = tibble::tibble()) {
  stopifnot(tibble::is_tibble(x))

  structure(x,class = c("excursionCore",class(tibble::tibble())))
}

new_shift <- function(x = list()) {
  stopifnot(tibble::is_tibble(x))

  structure(x,class = c("shift",class(tibble::tibble())))
}

new_shiftCore <- function(x = tibble::tibble()) {
  stopifnot(tibble::is_tibble(x))

  structure(x,class = c("shiftCore",class(tibble::tibble())))
}
