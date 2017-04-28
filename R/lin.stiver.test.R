# Question 4
# MLE based test of Lin and Stivers under heteroscedasticity

#' MLE based test of Lin and Stivers
#' @description Perform two sample MLE based t test of Lin and Stivers on vectors of data.
#' @usage lin.stivers.test(x, y, alternative = c("two.sided", "less", "greater"))
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
#' lin.stivers.test(vec1, vec2, alternative = "greater")
#'
lin.stivers.test <- function(x, y, alternative = "two.sided") {
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


  r <- var(t.paired, n.paired)/(sd(t.paired)*sd(n.paired))
  f <- n1*(n1+n3+n2*cov(t.paired, n.paired)/var(t.paired))/((n1+n2)*(n1+n3) - n2*n3*r*r)
  g <- n1*(n1+n2+n3*cov(t.paired, n.paired)/var(n.paired))/((n1+n2)*(n1+n3) - n2*n3*r*r)

  v1.numerator <- (f^2/n1 + (1-f)^2/n2)*var(t.paired)*(n1-1) + (g^2/n1 + (1-g)^2/n3)*var(n.paired)*(n1-1) - 2*f*g*cov(t.paired,n.paired)*(n1-1)/n1
  v1 <- v1.numerator/(n1-1)
  diff <- f*(mean(t.paired) - mean(t.rest)) - g*(mean(n.paired) - mean(n.rest)) + mean(t.rest) - mean(n.rest)

  z.ls <- diff/sqrt(v1)
  df <- n1

  p.value <- 0
  if(alternative == "two.sided") {
    p.value <- 2*(1- pt(abs(z.ls), df))
  }else if(alternative == "greater") {
    p.value <- 1 - pt(z.ls, df)
  }else if(alternative == "less") {
    p.value <- pt(z.ls, df)
  }

  cat("    lin stiver test\n P value is:", p.value, "\n statistic:", z.ls)
}
