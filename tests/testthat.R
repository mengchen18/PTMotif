library(testthat)
library(PTMotif)

data("motifExampleData")
a <- motifExampleData

v <- detectMotifs(a$background[1:10], a$background[1:100], min.seqs = 3, method = "all", ncores = 1)
v <- detectMotifs(a$background[1:3], a$background[1:100], min.seqs = 5, method = "all", ncores = 1)
v <- detectMotifs(a$background[1:10], a$background[1:11], min.seqs = 5, method = "all", ncores = 1)
v <- detectMotifs(a$background[1:5], a$background[1:5], min.seqs = 5, method = "all", ncores = 1)


v1 <- detectMotifs(
  a$background[1:10], a$background[1:100], 
  fg.genes = paste0("g", 1:10), bg.genes = paste0("g", 1:100), 
  min.seqs = 3, method = "all", ncores = 1, max.fdr = 0.5
)

v2 <- detectMotifs(
  a$background[1:10], a$background[1:100], 
  fg.genes = paste0("g", 1:10), bg.genes = paste0("g", 1:100), 
  min.seqs = 3, method = "cvp", ncores = 1, max.fdr = 0.5
)

v3 <- detectMotifs(
  a$background[1:10], a$background[1:100], 
  fg.genes = paste0("g", 1:10), bg.genes = paste0("g", 1:100), 
  min.seqs = 3, method = c("all", "cvp"), ncores = 1, max.fdr = 0.5
)

v4 <- detectMotifs(
  a$background[1:10], a$background[1:100], 
  fg.genes = paste0("g", rep(1, 10)), bg.genes = paste0("g", 1:100), 
  min.seqs = 3, method = c("all", "cvp"), ncores = 1, max.fdr = 0.5
)


seq <- c("ALNQKTSEKMKKRKMSNSFHGIRPPQLEQPE", 
         "REHSYVLSAAKKSTGSPTQETQAPFIAKRVE",
         "SSKKMGSIFDREDQASPRAGSLAALEKRQAE",
         "GELYDKSIIQSAQQDSIKKANMKRENKAYSF",
         "VHRDLKPENILYADDTPGAPVKIIDFGFARL",
         "MSLSAGSSPLHSPKITPHTSPAPRRRSHTPN")

v <- motif_all(seq, min.seqs = 6)


seq <- c("ALNQKTSEKMKKSTMSNSFHGIRPP______", 
         "REHSYVLSAAKKSTGSPTQETQAPF______",
         "SSKKMGSIFDRESTASPRAGSLAAL______",
         "GELYDKSIIQSASTDSIKKANMKRE______",
         "VHRDLKPENILYSTDTPGAPVKIID______",
         "MSLSAGSSPLHSPKITPHTSPAPRR______")

v <- motif_all(seq, min.seqs = 5)
v <- motif_all(seqs = a$background[1:10], min.seqs = 3)

