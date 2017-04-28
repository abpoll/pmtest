# Question 5
# weighted Z-test combination
# One limitation of the p-values pooling approach is in confidence interval estimation

#' Weighted Z-test combination
#' @description Perform two sample weighted z test on vectors of data.
#' @usage weighted.z.test(x, y, alternative = c("two.sided", "less", "greater"))
#' @param x a numeric vector of data values
#' @param y a numeric vector of data values.
#' @param alternative	a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @details
#' The formula interface is only applicable for the 2-sample tests.alternative = "greater" is the alternative that x has a larger mean than y.If paired is TRUE then both x and y must be specified and they must be the same length. Missing values are silently removed (in pairs if paired is TRUE). If var.equal is TRUE then the pooled estimate of the variance is used. By default, if var.equal is FALSE then the variance is estimated separately for both groups and the Welch modification to the degrees of freedom is used.If the input data are effectively constant (compared to the larger of the two means) an error is generated.

#' @return p.valuethe /  p value of the test
#' @return statistic aliases the value of the t-statistic.
#'
#' @examples
#' vec1 <- c(2,3,NA,4,3,4,NA,5,4,2)
#' vec2 <- c(3,NA,4,4,3,NA,NA,4,NA,2)
#' weighted.z.test(vec1, vec2, alternative = "greater")
#'
weighted.z.test <- function(x, y, alternative = "greater") {

  # get the paired sample, assign it to t.paired and n.paired
  t.paired <- x[(!is.na(x)) & (!is.na(y))]
  n.paired <- y[(!is.na(x)) & (!is.na(y))]
  n1 <- length(t.paired)  # number of paired sample

  # get the rest sample which is in T(not NA) and in N(is NA)
  t.rest <- x[(is.na(y)) & (!is.na(x))]
  n2 <- length(t.rest)

  # get the rest sample which is in N(not NA) and in T(is NA)
  n.rest <- y[(is.na(x)) & (!is.na(y))]
  n3 <- length(n.rest)

  p.1i <- t.test(t.paired, n.paired, paired = T, alternative = "greater")$p.value
  p.2i <- t.test(t.rest, n.rest, alternative = "greater")$p.value
  z.1i <- qnorm(1-p.1i)
  z.2i <- qnorm(1-p.2i)
  w1 <- sqrt(2*n1)
  w2 <- sqrt(n2+n3)  # by Zaykin[10]... square root of the sample sizes in practice
  p.ci <- 1 - pnorm((w1*z.1i + w2*z.2i)/sqrt(w1^2 + w2^2))

  p.value <- 0
  if(alternative == "two.sided") {
    if(p.ci < 0.5) {
      p.value <- 2*p.ci
    }else {
      p.value <- 2*(1 - p.ci)
    }
  }else if(alternative == "greater") {
    p.value <- p.ci
  }else if(alternative == "less") {
    p.value <- 1 - p.ci
  }

  cat("    weighted z test\n P value is:", p.value)
}


