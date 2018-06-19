gptest <- newGenopheno(nloci = 3, nalleles = c(2,3,2), alleleNames = list(c("+", "-"), c("A", "B", "0"), c(0,1)), rec = c(0.05,0.5))
gptest <- addLinkageGroup(gptest, nloci = 2, rec = 0.1)
gptest

gptest
