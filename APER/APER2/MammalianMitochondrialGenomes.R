library(ape)

###
### Chapter 3
###

mtgen.taxa <- read.table("mammal_mtGenome.fasta", skip = 4,
  nrows = 233, sep = "_", comment.char = "(", as.is = TRUE)
mtgen.taxa[1:5, ]
mtgen.taxa$V3 <- NULL
mtgen.taxa$V1 <- gsub("#", "", mtgen.taxa$V1)
mtgen.taxa$V1 <- gsub(" : ", "", mtgen.taxa$V1)
colnames(mtgen.taxa) <- c("code", "species")
mtgen.taxa[1:5, ]

mtgen <- read.dna("mammal_mtGenome.fasta",
                  format = "fasta", skip = 239)
length(mtgen)
names(mtgen)[1:10]

genes <- gsub("^[[:alnum:]]{1,}\\(", "", names(mtgen))
genes <- gsub("\\)$", "", genes)
unique(genes)[15]
i <- grep("Sequence does not exist", names(mtgen))
mtgen <- mtgen[-i]

genes <- gsub("^[[:alnum:]]{1,}\\(", "", names(mtgen))
genes <- gsub("\\)$", "", genes)
table(genes)

BF.sp <- matrix(NA, nrow = 233, ncol = 4)
rownames(BF.sp) <- mtgen.taxa$species
colnames(BF.sp) <- c("A", "C", "G", "T")
for (i in 1:233) {
    x <- grep(mtgen.taxa$code[i], names(mtgen))
    BF.sp[i, ] <- base.freq(mtgen[x])
}
matplot(BF.sp, type = "l", col = 1, xlab = "Species",
        ylab = "Base frequency")
legend(0, 0.22, c("A", "C", "G", "T"), lty=1:4, bty="n")

BF.gene <- matrix(NA, nrow = 39, ncol = 4)
rownames(BF.gene) <- unique(genes)
colnames(BF.gene) <- c("A", "C", "G", "T")
for (i in 1:nrow(BF.gene)) {
  x <- grep(rownames(BF.gene)[i], names(mtgen), fixed=TRUE)
  BF.gene[i, ] <- base.freq(mtgen[x])
}
par(mar = c(8, 3, 3, 2))
barplot(t(BF.gene), las = 2, legend = TRUE,
        args.legend = list(horiz = TRUE, bg = "white"))

cox1 <- mtgen[grep("COX1", names(mtgen))]
table(sapply(cox1, length))
cox1.muscle <- muscle(cox1)
cox1.mafft <- mafft(cox1)
dim(cox1.muscle)
dim(cox1.mafft)
nm <- rownames(cox1.muscle)
identical(cox1.muscle[, 1:1532], cox1.mafft[nm, 1:1532])
cox1.ali <- cox1.muscle[, 1:1532]

BF.cox1 <- matrix(NA, 3, 4)
rownames(BF.cox1) <- paste("codon position", 1:3)
colnames(BF.cox1) <- c("A", "C", "G", "T")
for (i in 1:3) {
    s <- rep(FALSE, 3)
    s[i] <- TRUE
    BF.cox1[i, ] <- base.freq(cox1.ali[, s])
}

par(mar = c(2.5, 4.1, 4.1, 1))
barplot(t(BF.cox1), main = "Cytochrome oxydase I",
        ylab = "Base frequency")
par(cex = 2)
text(0.7, BF.cox1[1, 1]/2, "A", col = "white")
text(0.7, BF.cox1[1,1]+BF.cox1[1,2]/2, "C", col = "white")
text(0.7, sum(BF.cox1[1, 1:2]) + BF.cox1[1, 3]/2, "G")
text(0.7, sum(BF.cox1[1, 1:3]) + BF.cox1[1, 4]/2, "T")

