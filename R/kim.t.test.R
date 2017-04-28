
#' Modified t-statistic of Kim et al
#' @description Perform two sample Kim t test on vectors of data.
#' @usage kim.t.test(x, y, alternative = c("two.sided", "less", "greater"))
#' @param x a (non-empty) numeric vector of data values
#' @param y a (non-empty) numeric vector of data values.
#' @param alternative	a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @details
#' --------------Needs to be filled...------------------
#' @return p.value p value of the test
#'
#' @return statistic aliases the value of the z-statistic.
#'
#' @examples
#' vec1 <- c(2,3,NA,4,3,4,NA,5,4,2)
#' vec2 <- c(3,NA,4,4,3,NA,NA,4,NA,2)
#' kim.t.test(vec1, vec2, alternative = "greater")
#'
#'
kim.t.test <- function(x, y, alternative = "two.sided") {

  if(!((all(is.numeric(x)|is.na(x)) & !all(is.na(x))) &
       (all(is.numeric(y)|is.na(y)) & !all(is.na(y))))) {
    stop("Each vector must have at least one numeric values")
  }

  # x is the first sample, y is the second sample
  # d is the difference between the paired samples where there are no NA values
  d <- x-y
  d <- d[!is.na(d)]
  n1 <- length(d)

  # t.rest is the independent sample from x without NAs
  t.rest <- x[(is.na(y)) & (!is.na(x))]
  n2 <- length(t.rest)

  # n.rest is the independent sample from y without NAs
  n.rest <- y[(is.na(x)) & (!is.na(y))]
  n3 <- length(n.rest)

  d.mu <- mean(d)
  t.mu <- mean(t.rest)
  n.mu <- mean(n.rest)
  # get n(h): harmonic mean of t.num and n.num
  n.h <- 2/(1/n2 + 1/n3)

  # as the paper says the null distribution of t is approximated with N(0,1), so we use z-table
  # calculate statistic
  sd <- sqrt(n1*var(d) + n.h^2*(var(n.rest)/n3 + var(t.rest)/n2))
  diff <- n1*mean(d) + n.h*(mean(t.rest) - mean(n.rest))
  t3 <- diff/sd

  p.value <- 0
  if(alternative == "two.sided") {
    p.value <- 2*(1- pnorm(abs(t3)))
  }else if(alternative == "greater") {
    p.value <- 1 - pnorm(t3)
  }else if(alternative == "less") {
    p.value <- pnorm(t3)
  } else {
    stop("arg should be one of \"two.sided\", \"greater\", \"less\"")
  }

  cat("    Modified t test of Kim et al\n      p-value is:", p.value, "\n       statistic:", t3)
}

