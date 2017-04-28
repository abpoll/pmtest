# This is a testing file...
doubleNA <- read.csv("test_double_NA.csv", header = T)
doubleNA.Typo <- read.csv("test_double_NA&Typo.csv", header = T)

# testing two samples with NA
x <- doubleNA$DiffCol
y <- doubleNA$Black

# testing two samples with NA and Typo(Nonnumerical)
#x <- doubleNA.Typo$DiffCol
#y <- doubleNA.Typo$Black

kim.t.test(x,y, "greater")
kim.t.test(x,y, "less")
kim.t.test(x,y, "two.sided")
kim.t.test(x,y,"typo")

lj.z.test(x,y, "greater")
lj.z.test(x,y, "less")
lj.z.test(x,y, "two.sided")
lj.z.test(x,y,"typo")

ekbohm.test(x,y, "greater")
ekbohm.test(x,y, "less")
ekbohm.test(x,y, "two.sided")
ekbohm.test(x,y,"typo")

lin.stivers.test(x,y, "greater")
lin.stivers.test(x,y, "less")
lin.stivers.test(x,y, "two.sided")
lin.stivers.test(x,y,"typo")

weighted.z.test(x,y, "greater")
weighted.z.test(x,y, "less")
weighted.z.test(x,y, "two.sided")
weighted.z.test(x,y,"typo")
