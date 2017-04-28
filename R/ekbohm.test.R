# Question 3
# MLE based test of Ekbohm under homoscedasticity
# when the variances of tumor and normal are equal, Ekbohm [5] suggested the following MLE based test statistic

#' MLE based test of Ekbohm
#' @description Perform two sample MLE based t test of Ekbohm on vectors of data.
#' @usage ekbohm.test(x, y, alternative = c("two.sided", "less", "greater"))
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
#' ekbohm.test(vec1, vec2, alternative = "greater")
#'
ekbohm.test <- function(x, y, alternative = "two.sided") {

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

  # Since this Ekbohm test is based on the same variance, we give warning if two samples have different variance at level 0.05
  x.na.rm <- x[!is.na(x)]
  y.na.rm <- y[!is.na(y)]
  f = var.test(x.na.rm, y.na.rm)
  if(f$p.value >= 0.05) {
    warning("Warning: The assumption of Ekbohm test is the same variance of two samples. Since the variance are different between two samples at α = 0.05, Lin Stiver test is recommended", noBreaks. = T)
  }

  r <- var(t.paired, n.paired)/(sd(t.paired)*sd(n.paired))
  f <- n1*(n1+n3+n2*r)/((n1+n2)*(n1+n3) - n2*n3*r^2)
  g <- n1*(n1+n2+n3*r)/((n1+n2)*(n1+n3) -n2*n3*r^2)
  var.hat <- (var(t.paired)*(n1-1) + var(n.paired)*(n1-1) + (1+r^2)*(var(t.rest)*(n2-1) + var(n.rest)*(n3-1)))/(2*(n1-1) + (1+r^2)*(n2+n3-2))
  v1.star <- var.hat*((2*n1*(1-r) + (n2+n3)*(1-r^2))/((n1+n2)*(n1+n3) - n2*n3*r^2))

  diff <- f*(mean(t.paired) - mean(t.rest)) - g*(mean(n.paired) - mean(t.rest)) + mean(t.rest) - mean(n.rest)
  z.e <- diff/sqrt(v1.star)
  df <- n1

  p.value <- 0
  if(alternative == "two.sided") {
    p.value <- 2*(1 - pt(abs(z.e), df))
  }else if(alternative == "greater") {
    p.value <- 1 - pt(z.e, df)
  }else if(alternative == "less") {
    p.value <- pt(z.e, df)
  }

  cat("    Ekbohm test\n P value is:", p.value, "\n statistic:", z.e)
}
