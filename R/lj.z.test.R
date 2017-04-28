
#' Corrected Z-test of Looney and Jones
#' @description Perform two sample z test of Looney and Jones on vectors of data.
#' @usage lj.z.test(x, y, alternative = c("two.sided", "less", "greater"))
#' @param x a (non-empty) numeric vector of data values
#' @param y a (non-empty) numeric vector of data values.
#' @param alternative	a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @details
#' --------------Needs to be filled...------------------
#' @return p.value  p value of the test
#' @return statistic aliases the value of the z-statistic.
#'
#' @examples
#' vec1 <- c(2,3,NA,4,3,4,NA,5,4,2)
#' vec2 <- c(3,NA,4,4,3,NA,NA,4,NA,2)
#' lj.z.test(vec1, vec2, alternative = "greater")
#'
lj.z.test <- function(x, y, alternative = "two.sided") {

  # get the paired sample, assign it to t.paired and n.paired
  t.paired <- x[(!is.na(x)) & (!is.na(y))]
  n.paired <- y[(!is.na(x)) & (!is.na(y))]
  n1 <- length(t.paired)  # number of paired sample

  # t.rest is the independent sample from x without NAs
  t.rest <- x[(is.na(y)) & (!is.na(x))]
  n2 <- length(t.rest)

  # n.rest is the independent sample from y without NAs
  n.rest <- y[(is.na(x)) & (!is.na(y))]
  n3 <- length(n.rest)

  #calculate mean.star and var.star
  t.mu.star <- mean(c(t.paired,t.rest))
  n.mu.star <- mean(c(n.paired,n.rest))
  t.var.star <- var(c(t.paired,t.rest))
  n.var.star <- var(c(n.paired,n.rest))
  paired.cov <- cov(t.paired,n.paired)

  sd <- sqrt(t.var.star/(n1+n2) + n.var.star/(n1+n3) - 2*n1*paired.cov/((n1+n2)*(n1+n3)))
  z.corr = (t.mu.star - n.mu.star)/sd

  p.value <- 0
  if(alternative == "two.sided") {
    p.value <- 2*(1- pnorm(abs(z.corr)))
  }else if(alternative == "greater") {
    p.value <- 1 - pnorm(z.corr)
  }else if(alternative == "less") {
    p.value <- pnorm(z.corr)
  } else{
    stop("arg should be one of \"two.sided\", \"greater\", \"less\"")
  }

  cat("    Corrected Z-test of Looney and Jones\n          P value is:", p.value, "\n           statistic:", z.corr)
}
