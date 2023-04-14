##### Chapter 1 #####

ls() # this is an R command

x <- rnorm(1e6)
x


##### Chapter 2 #####

library(ape)
ref <- read.GenBank("U20753")
dir.create("ref/")
setwd("ref/")
write.FASTA(ref, "refU.fas")

library(Rsubread)
buildindex("myref", "refU.fas")

myreads <- list(ref$U20753[1:100], ref$U20753[1001:1100],
                ref$U20753[10001:10100], rDNAbin(100)[[1]])
class(myreads) <- "DNAbin"
write.FASTA(myreads, "x.fas")

align("myref", "x.fas", type = "dna")

propmapped("x.fas.subread.BAM")

bam <- scanBam("x.fas.subread.BAM")
bam

sortBam("x.fas.subread.BAM", "x.sorted.bam")

asSam("x.fas.subread.BAM", "x.sorted.sam")

library(jackalope)
refjack <- read_fasta("refU.fas")
refjack
class(refjack)

pacbio(refjack, "pacbio_reads", n_reads = 100)

library(Biostrings)
readDNAStringSet("pacbio_reads_R1.fq", "fastq")

read.fastq("pacbio_reads_R1.fq")

sublong("myref", "pacbio_reads_R1.fq", "pacbio_reads_R1.bam")

setwd("../")


##### Chapter 3 #####

X <- matrix(1:3, 3)
rownames(X) <- LETTERS[1:3]
d <- dist(X)
str(d)

library(pegas)
x <- data.frame(L1 = c("A/A", "A/a", "a/a"))
x <- as.loci(x)
x
print(x, details = TRUE) # like View(x) but in the console
str(x)

x$population <- factor(c(1, 1, 2))
x
str(x)

library(adegenet) # normally loaded with pegas
y <- loci2genind(x) # function in pegas
y

y@tab
slotNames(y)

pop(y)
ploidy(y)

z <- list(Ind1 = rep(1, 1e6), Ind2 = rep(0, 1e6))
z <- new("genlight", z)
z

str(z)

zs <- new("SnpMatrix", as.raw(0:3))
zs
str(zs)

S <- matrix(c("A", "A", "A", "A", "A", "G"), 3, 2)
rownames(S) <- paste0("Ind", 1:3)
S
S <- as.DNAbin(S)
S

alview(S)

Sb <- as.DNAbin(list(Ind1=c("A", "A"), Ind2=c("A", "A", "G")))
Sb
str(Sb)

library(Biostrings)
DNAString("aaa")
DNAString("aaa") == DNAString("AAA")
DNAString("xaa")
BString("xaa")

library(SNPRelate)
snpgdsVCF2GDS("sampletwo.vcf", "sampletwo.gds")
genotype   { Bit2 2x2, 1B } *

samp <- snpgdsOpen("sampletwo.gds")
samp
   [  ] *

snpgdsSNPRateFreq(samp, with.id = TRUE)

library(pinfsc50)
od <- setwd(system.file("extdata/", package = "pinfsc50"))
dir()
fl <- "pinf_sc50.vcf.gz"
file.size(fl)
info <- VCFloci(fl)

dim(info)
names(info)

VCFlabels(fl)

get(fl, env = pegas:::.cacheVCF)

info

info$INFO
info$FORMAT

getINFO(info)

pegas::getINFO(info, what = "MQ")

snp <- is.snp(info)

table(snp)

range(info$POS)

info$POS[!snp]

X <- read.vcf(fl, which.loci = which(!snp), quiet = TRUE)
X

cat(VCFheader(fl))

length(rangePOS(info, 1, 1e5))

length(selectQUAL(info))
length(selectQUAL(info, threshold = 200))

library(vcfR)
vcf <- read.vcfR(fl, verbose = FALSE)
vcf

str(vcf)

chrom <- create.chromR(vcf)
chromoqc(chrom)

setwd(od)

saveRDS(vcf, "vcf.rds")
vcf2 <- readRDS("vcf.rds")
identical(vcf, vcf2)

num <- c("V00001", "U15717", "KE721553", "WRONG")
SEQ <- read.GenBank(num)
SEQ
attr(SEQ, "species")
attr(SEQ, "description")

names(SEQ)

names(SEQ) <- attr(SEQ, "description")
write.FASTA(SEQ, file = "SEQ.fas")

myfiles <- Xplorefiles()
names(myfiles)


##### Chapter 4 #####

1.2:5.2

-6:-1

-1:2
-(1:2)

args(seq.default)

rep(1:3, times = 3) # or rep(1:3, 3)
rep(1:3, each = 2)
rep(1:3, times = 3, each = 2)
rep(1:3, each = 2, length.out = 3)

args(paste)

paste("x", 1:6)
paste("x", 1:6, sep = "")
paste("x", 1:6, sep = " <- ")

paste(LETTERS[24:26], 1:3, sep = "", collapse = " %*% ")

match("Z", LETTERS)

args(match)

DF1 <- data.frame(x = 1:3, z = 11:13)
DF2 <- data.frame(x = 3:1, z = 13:11)
row.names(DF1) <- paste0("Ind", 1:3)
row.names(DF2) <- paste0("Ind", 3:1)
DF1
DF2

o <- match(row.names(DF2), row.names(DF1))
o
identical(DF1, DF2[o, ])

all(row.names(DF2) %in% row.names(DF1))
all(row.names(DF1) %in% row.names(DF2))

1:3 + 1
1:3 + rep(1, 2)

x <- 1:2
y <- "a"
x[2] <- y
x

s <- x > 0
x[s]

s2 <- is.na(x)

x[!s2] # drop missing values
x[s & !s2] # drop negative values and missing ones

x <- c(TRUE, TRUE, FALSE, FALSE)
y <- c(TRUE, FALSE, TRUE, FALSE)
x & y
x | y
xor(x, y)

system.time(x <- numeric(1e9))
object.size(x)

system.time(y <- x)

system.time(x[1] <- 1)
system.time(x[2] <- 2)

system.time(L <- list(x = x, y = y))
object.size(L)

data(woodmouse)
write.FASTA(woodmouse, "woodmouse.fas")

x <- readDNAStringSet("woodmouse.fas")
y <- read.FASTA("woodmouse.fas")
class(x)
class(y)
identical(as.DNAbin(x), y)

write.FASTA(y, "y.fas")
x2 <- readDNAStringSet("y.fas")
identical(x, x2)

num <- paste0("KX224", 490:529)
library(ape)
catopuma <- read.GenBank(num)
catopuma

catopuma.ali <- muscle(catopuma)
library(ips)
catopuma.ali.mafft <- mafft(catopuma, path = "/usr/bin/mafft")
identical(catopuma.ali, catopuma.ali.mafft)

dim(catopuma.ali)
checkAlignment(catopuma.ali, plot = FALSE)

head(attr(catopuma, "description"))

saveRDS(catopuma.ali, "catopuma.ali.rds")
saveRDS(catopuma, "catopuma.rds")

geo <- read.delim("geo_droso.txt")
str(geo)

samples.info <- read.delim("igsr_samples.tsv")
str(samples.info)

fl <- "ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
labs <- VCFlabels(fl)
all(labs %in% samples.info$Sample.name)

i <- match(labs, samples.info$Sample.name)
vars <- c("Sex", "Population.code", "Superpopulation.code")
DATA <- samples.info[i, vars]
row.names(DATA) <- samples.info$Sample.name[i]
str(DATA)

saveRDS(DATA, "DATA_G1000.rds")

odir <- setwd(system.file("files/", package = "adegenet"))
dir(pattern = "H1N1")

H1N1.HA <- read.dna("pdH1N1-HA.fasta", "fasta")
H1N1.HA
H1N1.NA <- read.dna("pdH1N1-NA.fasta", "fasta")
H1N1.NA
H1N1.DATA <- read.csv("pdH1N1-data.csv", as.is = TRUE)
str(H1N1.DATA)

all(labels(H1N1.HA) == H1N1.DATA$HA.acc.nb)
all(labels(H1N1.NA) == H1N1.DATA$NA.acc.nb)

setwd(odir)

library(pegas)
data(jaguar)
jaguar
names(jaguar)

HP <- read.FASTA("BIGSdb_gene-by-gene_alignment.fasta")
HP

EH.sam <-
"Mesocosm_EH.aligning.unique.sorted.MarkedDuplicates.sam"
SH.sam <-
"Mesocosm_SH.aligning.unique.sorted.MarkedDuplicates.sam"
Mock.sam <- "Mock.aligning.unique.sorted.MarkedDuplicates.sam"
JC.sam <-
"JudayCreek.aligning.unique.sorted.MarkedDuplicates.sam"
ref.fas <-
"mitogenomes_all_corrected_reading_frame.v0301.fasta"

library(Rsubread)
propmapped(EH.sam)

ref <- read.FASTA(ref.fas)
ref


##### Chapter 5 #####

X <- as.loci(data.frame(L1 = c("A/A", "A/A", "G|A"),
                        L2 = c("C/C", "C/T", "T|C")))
s <- summary(X)
s

str(s)

s[["L1"]]

class(s)
plot(s, layout = 4, col = c("grey", "white"))

Xg <- loci2genind(X)
summary(Xg)

summary(dat)

loci2genind(dat)

loci2genind(dat, na.alleles = "")

snpgdsSNPRateFreq(datsnp)

h <- haplotype(S)
h

class(h)

nh <- summary(h)
nh

sort(h)
sort(h, what = "labels")

sort(h, decreasing = FALSE)
sort(h, what = "labels", decreasing = TRUE)

subset(h, minfreq = 2)

str(h)

methods(haplotype)

hX <- haplotype(X)
hX

hap.div(h, variance = TRUE)

hap.div(S, variance = TRUE)

nuc.div(S)
nuc.div(h)

nuc.div(x, variance = TRUE)

dist.hamming(X)
dist.hamming(S)

dist.asd(X)

data(jaguar)
JAG2 <- jaguar[1:4, 2]
print(JAG2, details = TRUE)

library(poppr)
bruvo.dist(loci2genind(JAG2), replen = 4)

dist.asd(JAG2)

X$population <- factor(paste0("Pop", c(1, 1, 2)))
X
print(X, details = TRUE)

by(X)

gr <- gl(2, 2, labels = paste0("pop", 1:2))
gr

o <- outer(gr, gr, "==")
o

o[lower.tri(o)]

o <- outer(gr, gr, function(x, y) x == "pop1" & y == "pop1")
o[lower.tri(o)]

foo <- function(gr, g1, g2, matrix = TRUE) {
    o <- outer(gr, gr, function(x, y)
        x == g1 & y == g2 | x == g2 & y == g1)
    if (matrix) return(o)
    o[lower.tri(o)] | o[upper.tri(o)]
}

foo(gr, "pop1", "pop2")

foo(gr, "pop1", "pop2", matrix = FALSE)

which(foo(gr, "pop1", "pop2"), arr.ind = TRUE)

y <- rnorm(800, rep(0:1, each = 400))
plot(y, type = "l")
lines(runmed(y, k = 201), type = "l", lwd = 3, col = "grey")
legend("topleft", legend = c("y", "runmed(y, k = 201)"),
       lwd = c(1, 3), col = c("black", "grey"))

data(woodmouse)
sw(woodmouse)
sw(woodmouse, rowAverage = TRUE)
sw(woodmouse, 200, 200, rowAverage = TRUE)

sw.wood <- sw(woodmouse, rowAverage = TRUE)
plot(sw.wood, show.ranges = TRUE, col.ranges = "black",
     ylab = "GC content")

snpgdsSlidingWindow(samp, winsize = 1, shift = 1,
                    FUN = function(...) NULL)

Z <- data.frame(L1 = c("A/A", "A/A", "G/G", "G/G"),
                L2 = c("A/A", "G/G", "A/A", "G/G"),
                L3 = c("A/G", "A/A", "A/G", "G/G"),
                L4 = c("A/G", "A/G", "A/G", "A/G"))
Z <- as.loci(Z)
Z.genind <- loci2genind(Z)

Z.genind@tab

pca.Z <- dudi.pca(Z.genind@tab, scannf = FALSE, nf = 2)
pca.Z

biplot(pca.Z)

d.Z <- dist.asd(Z)
d.Z

Zgds <- snpgdsOpen("Z.gds")
pca2.Z <- snpgdsPCA(Zgds)

pca2.Z

plot(pca2.Z)

library(flashpcaR)
args(flashpca)

res.flash <- flashpca(Z.genind@tab, ndim = 1)
str(res.flash)

mds.Z <- cmdscale(d.Z)
hc.Z <- hclust(d.Z)
layout(matrix(1:2, 1))
plot(mds.Z, type = "n")
text(mds.Z, labels = 1:4)
plot(hc.Z, hang = -1)

catopuma.ali <- readRDS("catopuma.ali.rds")

base.freq(catopuma.ali, all = TRUE)

h <- haplotype(catopuma.ali)
nrow(h)
nrow(catopuma.ali)

ss.catopuma <- seg.sites(catopuma.ali)
head(ss.catopuma)

d.K80 <- dist.dna(catopuma.ali)
d.raw <- dist.dna(catopuma.ali, "raw")
plot(d.raw, d.K80)
abline(0, 1)

summary.default(d.K80)
summary.default(d.raw)

hist(d.raw)
rug(d.raw)

fl <- "global.pop.GATK.SNP.hard.filters.V3.phased_all.pop.maf.05.recode.vcf.gz"
info.droso <- VCFloci(fl)
info.droso

table(info.droso$CHROM)

SNP <- is.snp(info.droso)
table(SNP)

geo <- read.delim("geo_droso.txt")
head(geo)

labs <- VCFlabels(fl)
all(geo$ID == labs)

table(geo$Region)

nonsnp <- which(!SNP)
droso.nonsnp <- read.vcf(fl, which.loci = nonsnp)
a <- getAlleles(droso.nonsnp)
a
all(unlist(lapply(a, nchar)) == 1)

table(lengths(a))

res <- by(droso.nonsnp, geo$Region)

unique(names(res))

chr.nonsnp <- info.droso$CHROM[nonsnp]
pos.nonsnp <- info.droso$POS[nonsnp]
names(res) <- paste(chr.nonsnp, pos.nonsnp, sep = ".")
res[1:3] # print the first three loci

table(chr.nonsnp)

library(SNPRelate)
snpgdsVCF2GDS(fl, "drosoSNP.gds", method = "biallelic.only")
genotype   { Bit2 121x1047913, 30.2M } *

snp.gds <- snpgdsOpen("drosoSNP.gds")
pca.snp <- snpgdsPCA(snp.gds, autosome.only = FALSE)

barplot(pca.snp$eigenval[1:10])

pca.snp$varprop

X <- loci2genind(droso.nonsnp)@tab
dim(X)
dim(X)

7878 * 3 + 27 * 4

pca.nonsnp <- dudi.pca(X, scannf = FALSE, nf = 2)

pca.nonsnp$eig/sum(pca.nonsnp$eig)

chromosomes <- sort(unique(info.droso$CHROM))
pca <- vector("list", 5)
names(pca) <- chromosomes

for (chr in chromosomes) {
    subx <- which(info.droso$CHROM == chr)
    j <- subx[seq(1, length(subx), length.out = 1e4)]
    dat <- loci2genind(read.vcf(fl, which.loci = j))
    pca[[chr]] <- dudi.pca(dat@tab, scannf = FALSE, nf = 2)
}

layout(matrix(1:6, 3, 2, byrow = TRUE))
for (i in 1:5)
    screeplot(pca[[i]], npcs = 10,
          main = paste("Chromosome", chromosomes[i]), las = 1)

DATA <- readRDS("DATA_G1000.rds")
lapply(DATA, table)

fl <- "ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz"
labs <- VCFlabels(fl)
length(labs)

all(labs %in% samples.info$Sample.name)

info <- VCFloci(fl)
info

plot(sw(rep(1, nrow(info)), 200, 200, POS=info$POS, FUN=sum))

table(is.snp(info))

MITO <- read.vcf(fl, to = nrow(info))
i <- match(row.names(MITO), samples.info$Sample.name)
MITO$population <- samples.info$Population.code[i]
MITO$Continent <- samples.info$Superpopulation.code[i]

g <- loci2genind(MITO, ploidy = 1)
pca <- dudi.pca(g@tab, scannf = FALSE, nf = 2)
pcasvd <- prcomp(g@tab, scale. = TRUE)

head(pca$eig/sum(pca$eig))
head(pcasvd$sdev^2/sum(pcasvd$sdev^2))

layout(matrix(1:2, 1))
plot(pca$li, cex = 0.5, main = "PCA with ade4")
plot(pcasvd$x[, 1:2], cex = 0.5, main = "PCA by SVD")

hist(pca$co[, 1], 100, main = "ade4")
hist(pcasvd$rotation[, 1], 100, main = "SVD")

checkAlignment(H1N1.HA, plot = FALSE)

checkAlignment(H1N1.NA, plot = FALSE)

summary.default(dist.dna(H1N1.HA, "N"))
summary.default(dist.dna(H1N1.NA, "N"))

X <- cbind(H1N1.HA, H1N1.NA, check.names = FALSE)
rownames(X) <- H1N1.DATA$X
X

base.freq(X, freq = TRUE, all = TRUE)
h <- haplotype(X)
h

h.NA <- sort(haplotype(H1N1.NA))
h.HA <- sort(haplotype(H1N1.HA))
layout(matrix(c(1, 3, 2, 3), 2, 2))
plot(h.NA, las = 2, axisnames = FALSE, main = "NA")
plot(h.HA, las = 2, axisnames = FALSE, main = "HA")
plot(sort(h), las = 2, axisnames = FALSE, main = "Combined")

d <- dist.dna(X)
tr <- nj(d)
dates <- as.Date(H1N1.DATA$date)
str(dates)
all(H1N1.DATA$X == tr$tip.label)
plotTreeTime(tr, dates, color = FALSE, edge.width = 0.5)

lengths(getAlleles(jaguar))

lengths(getGenotypes(jaguar))

s <- summary(jaguar)
plot(s, what = "alleles", layout = 16, col = "grey", las = 2)

table(jaguar$population)

bypop <- by(jaguar)
bypop[c(1, 13)]

X <- loci2genind(na.omit(jaguar))
dim(X@tab)

acp.jaguar <- prcomp(X@tab, scaled. = TRUE)

screeplot(acp.jaguar, npcs = 40)

vasr <- acp.jaguar$sdev^2
vars/sum(vars)

pop <- jaguar$population
plot(acp.jaguar.svd$x[, 1:2], asp = 1, pch = (1:4)[pop])
legend("bottomleft", legend = levels(pop), pch = 1:4)

round(allelicrichness(jaguar), 1)
allelicrichness(jaguar, method = "raw")

rhost(jaguar)
rhost(jaguar, method = "rarefaction")

round(base.freq(HP, all = TRUE), 4)

HP <- as.matrix(HP)
foo <- function(x) base.freq(x, all = TRUE)["-"]
o4 <- sw(HP, width = 1e4, step = 1e4, FUN = foo)
dim(o4)

image(o4, col=grey((10:0)/10), axes=FALSE, xlab="Position (Mb)")
at <- seq(0.2, 1.6, 0.2)
axis(1, at = 1e6 * at/ncol(HP), labels = at)
box()
mtext("Sequences", 2, 1)

d <- dist.dna(HP, "N")
summary(as.vector(d))

d2 <- dist.dna(HP, "N", p = TRUE)
cor(d, d2)

mds <- cmdscale(d, 10, eig = TRUE)
barplot(mds$eig[1:10])

pco <- mds$points
layout(matrix(1:3, 3))
plot(pco[, 1:2], cex = 0.5, xlab = "Axis 1", ylab = "Axis 2")
plot(pco[, c(1, 3)], cex = 0.5, xlab = "Axis 1", ylab = "Axis 3")
plot(pco[, c(1, 4)], cex = 0.5, xlab = "Axis 1", ylab = "Axis 4")

i <- which(pco[, 2] > 1e4)
names(i)

library(Rsamtools)
asBam(EH.sam, "EH.bam")
EHreads <- scanBam("EH.bam")

EHtab <- table(EHreads[[1]]$rname)
EHtab

names(EHtab) <- gsub("_", " ", stripLabel(names(EHtab)))

barplot(EHtab, horiz = TRUE, las = 1, main = "EH")


##### Chapter 6 #####

G <- data.frame(L1 = c("A|A", "A|A", "G|G", "G|G"),
                L2 = c("C|C", "C|C", "T|T", "T|T"),
                L3 = c("C|C", "T|T", "C|C", "T|T"))
G <- as.loci(G)

LD(G)

LD(G, locus = c(1, 3))

H <- unphase(G)
LD2(H)
LD2(H, c(1, 3))

Y <- data.frame(L1=sample(c("A/A", "A/a", "a/a"), 5, rep=TRUE),
                L2=sample(c("B/B", "B/b", "b/b"), 5, rep=TRUE),
                L3=sample(c("C/C", "C/c", "c/c"), 5, rep=TRUE))
Y
Y <- as.loci(Y)

library(haplo.stats)
dat <- setupGeno(loci2alleles(Y))
dat

hapem <- haplo.em(dat)
hapem

fr <- hapem$hap.prob
names(fr) <- t(apply(hapem$haplotype, 1, paste, collapse="-"))
barplot(fr)

Gx <- loci2SnpMatrix(G) # function in pegas
library(snpStats)
rules <- snp.imputation(Gx)
rules

impute.snps(rules, Gx)

haplotype(G, locus = 1:3)

haplotype(G, locus = 1:3, compress = FALSE)

ldG <- LDscan(G, quiet = TRUE)
ldG

LDscan(G, depth = 1, quiet = TRUE)
LDscan(G, depth = 2, quiet = TRUE)

LDmap(ldG, col=grey(10:1/10), border=TRUE, scale.legend=0.2)

snpgdsVCF2GDS("G.vcf", "G.gds", verbose = FALSE)
Ggds <- snpgdsOpen("G.gds")

library(SNPRelate)
snpgdsLDMat(Ggds, slide = 0, method = "corr")
ld.rel

as.dist(ld.rel$LD)

library(snpStats)
ld.stat <- ld(loci2SnpMatrix(G), depth=2, stats="R.squared")
ld.stat

as.dist(t(as.matrix(ld.stat)))

X <- read.vcf(fl, which.loci = 1:nrow(info.droso))
any(!is.phased(X))

isnpx <- which(info.droso$CHROM == "X" & SNP)
length(isnpx)

layout(matrix(1:6, 3, 2, byrow = TRUE))
for (shift in 0:5 * 2e4) {
    sel <- isnpx[1:100 + shift]
    x <- read.vcf(fl, which.loci = sel, quiet = TRUE)
    s <- LDscan(x, quiet = TRUE)
    LDmap(s, info.droso$POS[sel], scale.legend = 5,
          col = grey(10:1/10))
}

x <- read.vcf(fl, which.loci = isnpx)
ldlong <- LDscan(x, depth = 1e5)
round(summary(ldlong[[1]]), 4)
hist(ldlong[[1]], main = "")

n <- 2 * 121
PrA <- 0.5
PrB <- 0.5

nrep <- 10000
r <- numeric(nrep)
for (i in 1:nrep) {
    SA <- sample(1:2, n, replace=TRUE, prob=c(PrA, 1 - PrA))
    SB <- sample(1:2, n, replace=TRUE, prob=c(PrB, 1 - PrB))
    PA <- sum(SA == 1)/n
    PB <- sum(SB == 1)/n
    PAB <- sum(SA == 1 & SB == 1)/n
    r[i] <- abs((PAB - PA*PB)/sqrt(PA*(1 - PA) * PB*(1 - PB)))
}

quantile(r, 0.95)

table(ldlong[[1]] > 0.125)

fl <- "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
library(SNPRelate)
snpgdsVCF2GDS(fl, "chr22.gds")

x <- snpgdsOpen("chr22.gds")
fx <- snpgdsSNPRateFreq(x)
str(fx)

table(cut(fx$MinorFreq, 0:5/10))

s <- which(fx$MinorFreq > 0.495)
sel <- read.gdsn(index.gdsn(x, "snp.id"))[s]
length(sel)

ld <- snpgdsLDMat(x, method = "r", snp.id = sel, slide = -1)

pos <- read.gdsn(index.gdsn(x, "snp.position"))[s]

LDmap(as.dist(ld$LD^2), pos/1e6, col = grey(9:0/9))

locnms <- names(jaguar)[1:13]

res <- list()
nms <- character()
for (i in 1:12) {
    for (j in (i + 1):13) {
        res <- c(res, list(LD2(jaguar, c(i, j))))
        nms <- c(nms, paste(locnms[i], locnms[j], sep = "-"))
    }
}
names(res) <- nms

Pvals <- sapply(res, function(x) x$T2["P-val"])
names(Pvals) <- nms
head(Pvals)

hist(Pvals, 20)
rug(Pvals)

Pmat <- matrix(NA, 13, 13, dimnames = list(locnms, locnms))
Pmat[lower.tri(Pmat)] <- Pvals >= 0.05

image(1:13, 1:13, t(Pmat), col=grey(1:2/3), xaxt="n", yaxt="n",
      xlab = "", ylab = "")
mtext(locnms, at = 1:13)
mtext(locnms, 2, at = 1:13, las = 1, adj = 1)

text(1:13, 1:13, lengths(getAlleles(jaguar)), font = 2)

layout(matrix(1:81, 9, 9, byrow = TRUE))
par(mar = c(0.1, 0.1, 1.3, 0.1))
for (i in 1:78) {
    image(res[[i]]$Delta, col=grey(0:4/4), axes=FALSE,
          xlab="", ylab="")
    mtext(names(res)[i], cex = .9,
          font = ifelse(Pvals[i] < 0.05, 2, 3))
}

n.alleles <- sapply(res, function(x) sum(dim(x$Delta)))
plot(n.alleles, Pvals)
abline(h = 0.05, lty = 2)


##### Chapter 7 #####

proba.genotype()
proba.genotype(alleles=c("A", "G"), p=c(0.9, 0.1), ploidy=4)

length(proba.genotype(1:20))
length(proba.genotype(1:20, ploidy = 4))

hw.test(unphase(X))

hw.test(as.loci(c("A/A/A/A", "T/T/T/T")))

Z$population <- gl(2, 2)

H(Z, variance = TRUE, observed = TRUE)

by(Z, FUN = H)

Fst(Z)

library(mmod)
Z.genind <- loci2genind(Z)
diff_stats(Z.genind)

D_Jost(Z.genind)

Z.bs <- chao_bootstrap(Z.genind)
Z.bs

test.Z.D <- summarise_bootstrap(Z.bs, D_Jost)
test.Z.D

library(snpStats)

Zsnp <- loci2SnpMatrix(Z)
Fst(Zsnp, Z$population)

Zgds <- snpgdsOpen("Z.gds")
snpgdsFst(Zgds, Z$population, "W&H02", verbose = FALSE)

X <- diag(3)
rownames(X) <- LETTERS[1:3]
X
d <- dist(X, "manhattan")
d

nt <- mst(d)
nt

str(nt)

mn <- msn(d)
mn

str(mn)

rnt <- rmst(d)
rnt

all.equal(rnt, mn)
all.equal(nt, mn)

hapX <- haplotype(X, labels = rownames(X))
ntcs <- haploNet(hapX, d)
all.equal(rnt, ntcs)

jn <- mjn(X)
jn

str(jn)
attr(jn, "data")

Y <- sample(c("A", "G"), size = 45, replace = TRUE)
Y <- matrix(Y, 9, 5)
Y
Yloc <- as.loci(as.data.frame(Y))
Yloc$population <- gl(3, 3, labels = paste0("pop", 1:3))
Yg <- loci2genind(Yloc, ploidy = 1)

res.dapc <- dapc(Yg, n.pc = 5, n.da = 2)
res.dapc

col <- c("grey30", "grey60", "grey90")
layout(matrix(1:2, 1))
scatter(res.dapc, col = col)
compoplot(res.dapc, col = col)

n <- 90; p <- 100; K <- 3
X <- sample(0:1, n * p, replace = TRUE)
dim(X) <- c(n, p)
X <- genind(X, ploidy = 1)
pop(X) <- gl(K, n/K, labels = paste0("Pop", 1:K))

dapc5 <- dapc(X, n.pca = 5, n.da = 2)
dapc80 <- dapc(X, n.pca = 80, n.da = 2)

layout(matrix(1:2, 1))
scatter(dapc5, col = col)
scatter(dapc80, col = col)

a.score(dapc5)$mean
a.score(dapc80)$mean

optim.a.score(dapc80)

snapclust(Yg, k = 2)

snapclust.choose.k(5, Yg)

snapclust.choose.k(5, Yg, IC = "AICc")

library(rhierbaps)
dat.snp <- load_fasta(as.DNAbin(Y))
str(dat.snp)

res.baps <- hierBAPS(dat.snp, max.depth = 2, n.pops = 5,
             assignment.probs = TRUE, quiet = TRUE)
res.baps

apply(res.baps$cluster.assignment.prob[[1]], 1, which.max)
res.snap$group

layout(matrix(1:2, 1))
compoplot(res.baps$cluster.assignment.prob[[1]], col = col[c(1, 3)])
compoplot(res.snap, col = col[c(1, 3)])

library(fastbaps)
Ybaps <- import_fasta_sparse_nt(as.DNAbin(Y))

Ybaps

Ybaps.opt <- optimise_prior(Ybaps)
Ybaps.opt$prior

res.fast <- fast_baps(Ybaps.opt, 2, quiet = TRUE)

plot(res.fast)

G <- sapply(Yloc[, -6], as.integer) - 1
write.table(t(G), "G.geno", sep = "", row.names = FALSE,
            col.names = FALSE)

K <- 1:3
res.snmf <- snmf("G.geno", K=K, repetitions=10, entropy=TRUE)

names(K) <- paste("K", K, sep = "=")
sapply(K, function(k) cross.entropy(res.snmf, k))

colMeans(sapply(K, function(k) cross.entropy(res.snmf, k)))

layout(matrix(1:2, 1))
for (i in 2:3) {
    o <- barchart(res.snmf, K = i, run = 10, sort.by.Q = FALSE,
                  space = 0, col = cols[1:i], paste("K =", i))
    mtext(rownames(Yloc), 1, at = 1:nrow(Yloc) - 0.5, las = 3)
}

remove.snmfProject("G.snmfProject")

Q <- list(snap = as.data.frame(res.snap$proba))

plotQ(Q)

library(admixturegraph)
leaves <- c("W", "Y", "X")
inner_nodes <- c("WY", "WYX")
edges <- parent_edges(c(edge("W", "WY"), edge("Y", "WY"),
         edge("WY", "WYX"), edge("X", "WYX")))
graph <- agraph(leaves, inner_nodes, edges)
graph

inner_nodes2 <- c("w", "y", "x", "XWY")
edges2 <- parent_edges(c(edge("W", "w"), edge("w", "XWY"),
                    edge("X", "x"), edge("x", "XWY"),
                    edge("Y", "y"),
                    admixture_edge("y", "w", "x", "alpha")))
graph2 <- agraph(leaves, inner_nodes2, edges2)

layout(matrix(1:2, 1))
plot(graph, col = "grey")
plot(graph2, col = "grey")

d <- dist.dna(catopuma.ali, "N")
nt <- rmst(d)
nt

all.equal(mst(d), nt, use.steps = FALSE)

msn(d)

plotNetMDS(nt, d, col = "black", font = 1)

xy <- replot()

library(network)
xy <- plot(as.network(nt))

xy <- list(x = xy[, 1], y = xy[, 2])

plot(nt, labels=FALSE, show.mutation=3, threshold=c(1, Inf))
replot(xy)

is <- which(SNP & info.droso$CHROM == "2L")
x <- read.vcf(fl, which.loci = sample(is, size = 1e4))
x$population  <- geo$Region
z <- loci2genind(x)
res <- dapc(z, n.pca = 20, n.da = 3)

scatter(res, col = "black")

xtabs(~ res$assign + x$population)

layout(matrix(1:6, 3, 2, byrow = TRUE))
for (i in 1:5)
    optim.a.score(res[[i]], main = names(res)[i])

x <- read.vcf(fl, which.loci = which(!SNP))
z <- loci2genind(x)

o <- snapclust.choose.k(10, z)
which.min(o)

barplot(o - o[6], xlab = "K", ylab = expression(delta*"AIC"))

snap.droso <- snapclust(z, k = 6)

xtabs(~ snap.droso$group + geo$Region)

summary(as.vector(snap.droso$proba))

vcf2geno(fl, "droso.geno")

K <- 1:10
droso.snmf <- snmf("droso.geno", K=K, repetitions=10,
                   entropy=TRUE)

droso.snmf <- load.snmfProject("droso.snmfProject")

names(K) <- paste("K", K, sep = "=")
foo <- function(k) cross.entropy(droso.snmf, k)[1:5]
sapply(K, foo)
round(colMeans(sapply(K, foo)), 5)

o <- barchart(droso.snmf, K = 5, run = 5, sort.by.Q = FALSE,
              space = 0, col = grey(0:4/4), "K = 5",
              ylab = "Membership probability")
mtext(geo$Region, 1, at = 1:121 - 0.5, las = 3, cex = 0.85)

sel <- seq(1, sum(SNP), length.out = 1e4)
x <- read.vcf(fl, which.loci = which(SNP)[sel])

x$population <- geo$Region
fbypop <- by(xsel)

POPS <- levels(x$population)
K <- length(POPS)
F2.droso <- matrix(0, K, K)

for (i in 1:(K - 1))
  for (j in (i + 1):K)
    F2.droso[i, j] <- F2.droso[j, i] <-
      F2(allele.freq=fbypop, jack=0, B=0, pops=POPS[c(i, j)])

dimnames(F2.droso) <- list(POPS, POPS)
F2.d <- as.dist(F2.droso)
F2.d

nt <- rmst(F2.d)
tr <- nj(F2.d)

d.HA <- dist.dna(h.HA, "N")
d.NA <- dist.dna(h.NA, "N")

nt.HA <- rmst(d.HA)
nt.NA <- rmst(d.NA)
nt.HA
nt.NA

table(rbind(nt.HA, attr(nt.HA, "alter.links"))[, 3])
table(rbind(nt.NA, attr(nt.NA, "alter.links"))[, 3])

layout(matrix(1:2, 2))
freq.HA <- summary(h.HA)
co.HA <- grey((max(freq.HA) - freq.HA)/max(freq.HA))
plot(nt.HA, labels = FALSE, bg = co.HA)
freq.NA <- summary(h.NA)
co.NA <- grey((max(freq.NA) - freq.NA)/max(freq.NA))
plot(nt.NA, labels = FALSE, bg = co.NA)

hw.test(jaguar)

Fst(jaguar)
Rst(jaguar)

Rst(jaguar)

library(mmod)
jaguar.genind <- loci2genind(jaguar)
diff_stats(jaguar.genind)

jaguar.bs <- chao_bootstrap(jaguar.genind, 1e4)
jaguar.test.Z.D <- summarise_bootstrap(jaguar.bs, D_Jost)
jaguar.test.Z.D


##### Chapter 8 #####

lon <- c(4, 7); lat <- c(43.6, 43.6)
geod(lon, lat)

geod(lon, c(0, 0))

dy <- dist.hamming(Yloc)
dy

res.amova <- amova(dy ~ population, data = Yloc)
res.amova

mds <- cmdscale(dy)
mds

o <- order(mds[, 1])
o

pop <- NULL # or: pop <- integer(nrow(Yloc))
pop[o] <- rep(1:3, each = 3)
pop <- factor(pop)
pop

amova(dy ~ pop)

vignette("MoranI")

cny <- chooseCN(xy, ask = FALSE, type = 3)
cny

spca.Y <- spca(Yg, cn=cny, scannf=FALSE, nfposi=2, nfnega=2)
spca.Y

plot(spca.Y, axis = 1)

geno <- loci2alleles(Yloc)
geno <- ifelse(geno == "A", 1, 2) - 1L
geno
xy <- as.matrix(xy)

library(tess3r)
res.tess <- tess3(X = geno, coord = xy, K = 1:4,
                  ploidy = 1, rep = 10)

plot(res.tess, xlab = "K", ylab = "Cross-validation score")

qmat <- qmatrix(res.tess, K = 2)
n <- 9
o <- barplot(qmat, border = NA, space = 0, xlab = "Individuals",
         ylab = "Ancestry proportions", palette.length = 2,
         col.palette = CreatePalette(c("grey", "lightgrey"), 2))
axis(1, at = 1:n - 0.5, labels = o$order, las = 3)

dir.create("Geneland_run/")
setwd("Geneland_run/")
library(Geneland)

geno <- geno + 1L

MCMC(xy, geno, path.mcmc = ".", nit = 1e4, npopmax = 4)

PostProcessChain(xy, "./", 50, 50, 0)

PlotTessellation(xy, "./")

Plotnpop("./", burnin = 9000)

setwd("../")

i <- which(SNP & info.droso$CHROM == "3R")
x <- read.vcf(fl, which.loci=sample(i, size=1e4), quiet=TRUE)
d <- dist.asd(x)
amova(d ~ Region/Locality, geo)

x <- MITO[, is.snp(MITO)]
x <- as.DNAbin(sapply(x, as.character))

checkAlignment(x, plot = FALSE)

dx <- dist.dna(x, "N")

am <- amova(dx ~ Continent/population, MITO, nperm = 100)
am

ic <- outer(MITO$Continent, MITO$Continent, "==")
ip <- outer(MITO$population, MITO$population, "==")
ic <- ic[lower.tri(ic)]
ip <- ip[lower.tri(ip)]

layout(matrix(1:4, 2, 2, byrow = TRUE))
hist(dx[ic], main = "Within continents")
hist(dx[!ic], main = "Between continents")
hist(dx[ip], main = "Within populations")
hist(dx[!ip], main = "Between populations")

conti <- levels(MITO$Continent)

layout(matrix(1:6, 3, 2, byrow = TRUE))
for (i in 1:5) {
  j <- foo(MITO$Continent, conti[i], conti[i], FALSE)
  hist(dx[ip & j],
       main = paste("Within populations in", conti[i]))


##### Chapter 9 #####

summary(rgeom(1e6, 1e-3))
summary(rexp(1e6, 1e-3))

tr <- rcoal(50)
plot(tr, type = "c", show.tip.label = FALSE)
axisPhylo()

tr$edge.length <- 2 * tr$edge.length

tr <- rcoal(3)
x <- simSeq(tr, 5, type = "USER", levels = 0:1)
x
as.character(x)

library(phyclust)
x <- ms(10, opts = "-T")
x

tr <- read.tree(text = x[3])
tr

theta.tree(tr)

res <- ms(5, opts = "-T -t 1")
res
str(res)

ms()

library(scrm)
scrm("5 2 -T")

scrm("5 2 -T -t 1.5")

scrm("3 2 -T -r 1 30")

subst <- sub_JC69(lambda = 1e-3)

tr <- rcoal(5)
vars <- create_variants(refjack, vars_phylo(tr), subst)
vars

illumina(vars, out_prefix = "illumina", n_reads = 5e4,
   read_length = 100, paired = FALSE, sep_file = TRUE)

data(jaguar)
theta.msat(jaguar)

tr <- rcoal(50)
res <- theta.tree(tr)
res

THETA <- seq(0.5, 2, 0.01)
log.lik <- theta.tree(tr, THETA, fixed = TRUE)
log.lik

plot(THETA, log.lik, type = "l")
abline(v = res$theta, lty = 3) # estimated THETA
abline(v = res$theta + c(-1.96, 1.96) * res$se, lty = 2)
abline(h = res$logLik - 1.95, lty = 4)
legend("bottomright", legend = expression("log-likelihood",
    hat(theta) * " (MLE)", "95%\ conf. interv.", "ML - 1.96"),
    lty = c(1, 3, 2, 4))

plot(skyline(rcoal(100)))

tr <- rcoal(100)
res <- mcmc.popsize(tr, nstep = 1e4, progress.bar = FALSE)
names(res)

N  <- extract.popsize(res, burn.in = 1e3, thinning = 10)
str(N)

plot(skyline(tr), lty = 3)
lines(N)
legend("bottomleft", legend = c("Skyline", "MCMC",
       "95% cred. int."), lty = c(3, 1, 1), lwd = c(1, 3, 1))

tr <- rcoal(100)
res.bnpr <- BNPR(tr)

plot_BNPR(res.bnpr)

th <- rtree(100)
res.bnpr.h <- BNPR(th)
plot_BNPR(res.bnpr.h)

site.spectrum(Yloc)

site.spectrum(Yloc, folded = FALSE, ancestral = rep("A", 5))

data(woodmouse)
sp <- site.spectrum(woodmouse)
sp
plot(sp, col = "lightgrey")

stw.wood <- stairway(sp, epoch = 1:14)
stw.wood

plot(stw.wood)

library(CubSFS)
res <- estimateCubSFS(sp, attr(sp, "sample.size"), n.knots=5,
             t_m = 1, alpha = 0.25, is.folded = TRUE)
plot(res$CoalRate[, 1:2], type = "l")
text(0.2, 1.5, expression(alpha==0.25))

x <- t(replicate(1000, branching.times(rcoal(50))))

td <- seq(0, 1.5, 1e-1)
N <- Popsicle(x, td)
N

plot(td, c(N, N[length(N)]), type="s", xlab="Time", ylab="Ne")
abline(h = 1, lty = 2)
legend("topleft", legend = "True Ne", lty = 2)

theta.s(catopuma.ali, var = TRUE)

nuc.div(catopuma.ali, var = TRUE)

library(coalescentMCMC)
o <- coalescentMCMC(catopuma.ali, 1e6, moves = c(1, 3))

effectiveSize(o)

autocorr.plot(o)

o.sub <- subset(o, 1e5, 1e3)
dim(o.sub)

autocorr.plot(o.sub)
effectiveSize(o.sub)

plot(o.sub)

summary(o.sub)
HPDinterval(o.sub)

sp <- site.spectrum(catopuma.ali)
plot(sp, col = "lightgrey", main = "")

tr <- nj(d.K80)
library(phangorn)
rtr <- midpoint(tr)
chr <- chronos(rtr)

res <- mcmc.popsize(chr, 1e4, progress.bar = FALSE)
N <- extract.popsize(res)
plot(skyline(chr), lty = 2)
lines(N)

s <- which(SNP & info.droso$CHROM == "2L")
droso <- read.vcf(fl, which.loci = s, quiet = TRUE)
xf <- site.spectrum(droso)
xf

reg <- geo$Region
SFS <- vector("list", nlevels(reg))
for (i in seq_along(res)) {
    ind <- reg == levels(reg)[i]
    SFS[[i]] <- site.spectrum(droso[ind, ])
}

layout(matrix(1:6, 2, 3, byrow = TRUE))
for (i in seq_along(SFS)) {
    plot(SFS[[i]], col = "lightgrey", main = "")
    title(levels(reg)[i], cex.main = 1.5)
}

url <- "ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.41_FB2011_09/fasta/dmel-all-chromosome-r5.41.fasta.gz"
reffile <- "dmel-all-chromosome-r5.41.fasta.gz"
download.file(url, reffile)

ref <- read.FASTA(reffile)

names(ref)

names(ref) <- gsub(" .*", "", names(ref))
s <- match(unique(info.droso$CHROM), names(ref))
ref2 <- ref[s]
ref2

library(psmcr)
x <- VCF2DNAbin(fl, ref2)
x <- seqBinning(x)

x

o <- psmc(x, "1+1+1+1+1+1", B = 100, trunksize=1e5)

plot(o, scaled = TRUE)

x22 <- VCF2DNAbin(fl, ref2, individual = 22)

library(maps)
map()
points(H1N1.DATA[, 4:5], pch = 20)

theta.s(H1N1.HA, variance = TRUE)
theta.s(H1N1.NA, variance = TRUE)
nuc.div(H1N1.HA, variance = TRUE)
nuc.div(H1N1.NA, variance = TRUE)

dgen <- dist.dna(H1N1.HA, "N", p = TRUE, as.matrix = TRUE)
dspat <- 20000 - geod(H1N1.DATA[, 4:5])
date <- as.POSIXct(H1N1.DATA$date)
res <- seqTrack(dgen, H1N1.DATA$X, date, prox.mat=dspat)
str(res)

map()
plotSeqTrack(res, H1N1.DATA[, 4:5], add = TRUE)

fgaps <- del.colgapsonly(HP, freq.only = TRUE)
sum(fgaps == 0)

x <- HP[, which(fgaps == 0)]

checkAlignment(x, plot = FALSE)

spx <- site.spectrum(x)
plot(spx, col = "lightgrey")


##### Chapter 10 #####

data(woodmouse)
dnds.woodm <- dnds(woodmouse, code = 2)
str(dnds.woodm)

summary(dnds.woodm)

library(phangorn)
data(woodmouse)
X <- dna2codon(phyDat(woodmouse))
X

trw <- nj(dist.dna(woodmouse))
m0 <- pml(trw, X)
m0
m1 <- optim.pml(m0, control = pml.control(trace = 0))
m1
anova(m0, m1)

Xall <- read.vcf(fl, which.loci = 1:nrow(info))

sel <- apply(is.phased(Xall), 2, all) & snp
hap <- haplotype(Xall[, sel], locus=1:sum(sel), compress=FALSE)
str(hap)

write.table(t(hap), "tmp.hap", quote=FALSE, col.names=FALSE)

write.table(info[sel, c(1, 2, 4, 5)], "tmp.map", quote=FALSE,
            col.names = FALSE)

library(rehh)
dat <- data2haplohh("tmp.hap", "tmp.map",
                    allele_coding = "none")

sc <- scan_hh(dat)
res <- ihh2ihs(sc, freqbin = 0.01)
str(res)

manhattanplot(ihs, pval = TRUE)

i <- which.max(res$ihs$IHS)
i

pos <- res$ihs$POSITION[i]
pos

foc <- which(dat@positions == pos)
foc

furc <- calc_furcation(dat, foc)
plot(furc, col = c("grey", "black"))

res.ehhs <- calc_ehhs(dat, mrk = foc)
str(res.ehhs)
plot(res.ehhs)

n <- 200
p <- 100
X <- sample(c("A/A", "A/G", "G/G"), n * p, replace = TRUE)
dim(X) <- c(n, p)
X <- as.loci(as.data.frame(X))
X$population <- gl(2, n/2)

Xg <- sapply(X[, 1:p], as.integer) - 1

FST <- MakeDiploidFSTMat(Xg, locusNames = 1:p,
                         popNames = X$population)
str(FST)

outf <- OutFLANK(FST, NumberOfSamples = n)
str(outf)

P <- pOutlierFinderChiSqNoCorr(FST, Fstbar = outf$FSTNoCorrbar,
                               dfInferred = outf$dfInferred)
str(P)

plot(-log10(P$pvalues), type = "h")

library(MINOTAUR)
MINOTAUR()

x <- matrix(runif(50), 5, 10)
library(poolSeq)
estimateSH(x, t = 0:9, Ne = 100, haploid = TRUE, h = 0.5)

xs <- t(apply(x, 1, sort))
estimateSH(xs, t = 0:9, Ne = 100, haploid = TRUE, h = 0.5)

layout(matrix(1:2, 1))
matplot(t(x), type = "l", lty = 1, col = 1)
matplot(t(xs), type = "l", lty = 1, col = 1)

tajima.test(catopuma.ali)

f <- function(x) tajima.test(x)$Pval.beta
sw(catopuma.ali, 1e3, 1e3, FUN = f, rowAverage = TRUE)

R2.test(catopuma.ali, plot = FALSE)

x <- read.delim("mtGenome_Catopuma_KX224490.txt")
str(x)

i <- grep("COX1", x$seq)
s <- x$start[i]
e <- x$end[i]
summary.default(dist.aa(trans(catopuma.ali[, s:e], 2)))

pval <- snmf.pvalues(droso.snmf, TRUE, entropy = TRUE,
                      K = 5, ploidy = 2)
str(pval)

y <- -log10(pval$pvalues)

CHR <- info.droso$CHROM[SNP]
POS <- info.droso$POS[SNP]
chr <- unique(CHR)

layout(matrix(1:6, 2, byrow = TRUE))
for (i in chr) {
   s <- CHR == i
   plot(POS[s]/1e6, y[s], type="h", ylim=c(0, 300), main=i,
        xlab="Position (Mb)", ylab=expression(-log[10](P)))

x <- read.pcadapt(fl, "vcf", "matrix", allele.sep = "|")

res.pcaa <- pcadapt(x, K = 5)
str(res.pcaa)

plot(res.pcaa)

plot(res.pcaa, option = "scores")

s <- which(info.droso == "X" & SNP)
X <- read.vcf(fl, which.loci = s)

any(!is.phased(X))

write.table(t(hap), "X_droso.hap", quote = FALSE,
            col.names = FALSE)
write.table(info.droso[s, c(1, 2, 4, 5)], "X_droso.map",
            quote = FALSE, col.names = FALSE)

dat <- data2haplohh("X_droso.hap", "X_droso.map",
                    allele_coding = "none")

sc <- scan_hh(dat, threads = 2)
ihs <- ihh2ihs(sc, freqbin = 0.01)

manhattanplot(ihs, pval = TRUE)

i <- which.max(ihs$ihs$IHS)
pos <- ihs$ihs$POSITION[i]
pos
foc <- which(dat@positions == pos)

furc <- calc_furcation(dat, foc)
plot(furc, col = c("grey", "black"))

plot(calc_ehhs(dat, foc))

G <- LEA::read.geno("droso.geno")
dim(G)

tajima.test(H1N1.HA)
tajima.test(H1N1.NA)

R2.test(H1N1.HA, B = 100, plot = FALSE)
R2.test(H1N1.NA, B = 100, plot = FALSE)

for (i in 1:3)
alview(trans(h.NA, codonstart = i)[1, 1:60], showpos = FALSE)

dnds.ha <- dnds(h.HA, codonstart = 3, quiet = TRUE)
dnds.na <- dnds(h.NA, codonstart = 3, quiet = TRUE)
summary.default(dnds.ha)
summary.default(dnds.na)

layout(matrix(1:2, 1))
hist(dnds.ha[dnds.ha >= 0], 50, main = "HA")
hist(dnds.na[dnds.na >= 0], 50, main = "NA")

aa.NA <- trans(h.NA, codonstart = 3)
aa.HA <- trans(h.HA, codonstart = 3)

summary.default(dist.aa(aa.NA))
summary.default(dist.aa(aa.HA))

sum(!duplicated(as.character(aa.NA)))
sum(!duplicated(as.character(aa.HA)))


##### Appendix A #####

mirrors <- getCRANmirrors()
dim(mirrors)
names(mirrors)

url.cran.paris <- mirrors$URL[grep("Paris", mirrors$City)]
url.cran.paris

repos <- c(CRAN = url.cran.paris)
repos

repos <- c(repos,
  BioCsoft = "https://bioconductor.org/packages/release/bioc")

pkgs <- c("adegenet", "ape", "Biostrings", "pegas",
          "SNPRelate", "snpStats")
install.packages(pkgs, dependencies = TRUE, repos = repos)

pkgs <- c("blwaltoft/CubSFS", "gtonkinhill/fastbaps",
          "gabraham/flashpca/tree/master/flashpcaR",
          "whitlock/OutFLANK", "mdkarcher/phylodyn",
          "ThomasTaus/poolSeq", "emmanuelparadis/psmcr"
          "rwdavies/STITCH", "bcm-uga/TESS3_encho_sen")
library(remotes)
install_github(pkgs)

install.packages("Geneland_4.0.7.tar.gz", repos = NULL)


##### Appendix B #####

library(ape)
n <- 1000 # number of sequences
s <- 10000 # sequence length

x <- rDNAbin(nrow = n, ncol = s)

system.time(saveRDS(x, "x.rds"))
system.time(write.FASTA(x, "x.fas"))

system.time(system("gzip x.fas"))
system.time(system("bzip2 x.fas"))

system.time(a <- readRDS("x.rds"))
system.time(b <- read.FASTA(gzfile("x.fas.gz")))
system.time(c <- read.FASTA(bzfile("x.fas.bz2")))
system.time(d <- read.FASTA("x.fas"))

library(Biostrings)
system.time(e <- readDNAStringSet("x.fas"))
system.time(f <- readDNAStringSet("x.fas.gz"))

system.time(g <- as.DNAbin(e))

system.time(writeXStringSet(e, "x.bios.fas"))
system.time(writeXStringSet(e, "x.bios.fas.gz", compress = TRUE))

fs0 <- file.size("x.fas")
file.size("x.rds") / fs0
file.size("x.fas.gz") / fs0
file.size("x.fas.bz2") / fs0
