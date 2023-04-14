library(ape)

###
### Chapter 3
###

x <- paste("AY66", 6597:7060, sep = "")
x <- c(x, "AY724411", "AY724412")
astraptes.seq <- read.GenBank(x)
table(sapply(astraptes.seq, length))
astraptes.seq.ali <- clustal(astraptes.seq)

table(attr(astraptes.seq, "species"))
taxa.astraptes <- attr(astraptes.seq, "species")
names(taxa.astraptes) <- names(astraptes.seq)
save(astraptes.seq.ali, taxa.astraptes,
     file = "astraptes.RData")
###
### Chapter 5
###

M.astraptes.K80 <- dist.dna(astraptes.seq.ali, p = TRUE)
summary(M.astraptes.K80)
hist(M.astraptes.K80)

tr <- nj(M.astraptes.K80)
tr$tip.label <- taxa.astraptes[tr$tip.label]

taxon <- unique(taxa.astraptes)
L <- vector(mode = "list", length = 10)
for (i in 1:10)
    L[[i]] <- grep(taxon[i], tr$tip.label)
zoom(tr, L)
