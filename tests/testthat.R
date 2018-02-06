library(testthat)
library(PTMotif)

test_check("PTMotif")


data("motifExampleData")
a <- motifExampleData

v <- detectMotifs(a$background[1:10], a$background[1:100], min.seqs = 3, method = "all", ncores = 4)
v <- detectMotifs(a$background[1:3], a$background[1:100], min.seqs = 5, method = "all", ncores = 4)
v <- detectMotifs(a$background[1:10], a$background[1:11], min.seqs = 5, method = "all", ncores = 4)
