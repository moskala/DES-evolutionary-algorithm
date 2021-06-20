source("des.R")
source("ackley.R")
source("shubert.R")
source("shekel.R")
source("rastr.R")


des(c(10, 10, 10, 10, 10, 10, 10), ackley, lower=-30, upper=30)

des(c(0,0), shubert, lower=-10, upper=10)

des(c(0,0), shekel, lower=0, upper=10)

des(c(4, 4, 4, 4), rastr, lower=-5.12, upper=5.12)
