install.packages("../PTMotif", type = "source", repos = NULL)



library(PTMotif)
data("motifExampleData")
a <- motifExampleData

v <- detectMotifs(a$foreground, a$background, method = "all")
vv <- v[grep("\\.", v$motif), ]



load("R/sysdata.rda")


s <- sapply(knwonMotifs$motif, function(x) grep(x, vv$motif))




length(unique(knwonMotifs$rexp))
