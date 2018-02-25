#Notes March 30
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("GenomicRanges")
library(IRanges)
rng <- IRanges(start=4, end=13)
rng
rng2 <- IRanges(start=4, width=3)
rng2
x <- IRanges(start=c(4, 7, 2, 20), end=c(13, 7, 5, 23))
x
names(x) <- letters[1:4]
x
class(x)
width(x)
start(x)
end(x)
end(x) <- end(x)+4
end(x)
x
range(x)
x[2:3]
start(x) < 5
x[start(x) < 5]
x[width(x) > 8]
a <- IRanges(start=7, width=4)
a
b <- IRanges(start=2, end=5)
b
c <- c(a, b)
c
class(c)
x <- IRanges(start=c(40, 80), end=c(67, 114))
x + 4L #actually increases total range by 8
x - 10L
y <- IRanges(start=c(4, 6, 10, 12), width=13)
y
restrict(y, 5, 10) #takes out all parts of ranges that don't fall within the restrict range of 5 to 10
flank(x, width=7) #creates upstream sequence to investigate things like promoters
flank(x, width=7, start=FALSE) #can specify forward versus revers
promoters <- flank(x, width=20)
promoters

set.seed(0)
alns <- IRanges(start=sample(seq_len(50), 20), width=5)
head(alns, 10)

reduce(alns) #collapse into super reads
gaps(alns) #pull out the gaps in your alignment

a <- IRanges(start=4, end=13)
b <- IRanges(start=12, end=17)
union(a,b)
a
b
intersect(a,b)
setdiff(a,b)

qry <- IRanges(start=c(1, 26, 19, 11, 21, 7), end=c(16, 30, 19, 15, 24, 8),
               names=letters[1:6])
sbj <- IRanges(start=c(1, 19, 10), end=c(5, 29, 16), names=letters[24:26])
qry
sbj
hts <- findOverlaps(qry,sbj) #query seq must always come first, subject seq second
class(hts)
hts
names(qry)[queryHits(hts)] #give names of which ranges map to which
names(sbj)[subjectHits(hts)]
#findOverlaps() only gives whether any overlap is present, not whether all of our range resides in the entirity of the other
hts_within <- findOverlaps(qry, sbj, type="within") #specify that you want entire overlap
hts_within

findOverlaps(qry, sbj, select="first")

sbj_it <- IntervalTree(sbj) #look for matches but limit out the ones in the query that aren't in same region as the subject
class(sbj_it)
findOverlaps(qry, sbjit)

as.matrix(hts)
countQueryHits(hts) #gives number of subject hits for each query range
setNames(countQueryHits(hts), names(qry))
countSubjectHits(hts)
ranges(hts, qry, sbj) #shows ranges of overlap between query and subject
countOverlaps(qry, sbj)
subsetByOverlaps(qry, sbj)

qry <- IRanges(start=6, end=13, name='query')
sbj <- IRanges(start=c(2, 4, 18, 19), end=c(4, 7, 21, 24), names=1:4)
nearest(qry, sbj)
sbj
precede(qry, sbj)
follow(qry, sbj)

qry <- IRanges(sample(seq_len(1000), 5), width=10)
sbj <- IRanges(sample(seq_len(1000), 5), width=10)
distanceToNearest(qry, sbj)
distance(qry, sbj) #compares them in order
qry
sbj
library(GenomicRanges)
gr <- GRanges(seqname=c("chr1", "chr1", "chr2", "chr3"),
              ranges=IRanges(start=5:8, width=10),
              strand=c("+", "-", "-", "+"))
gr
gr <- GRanges(seqname=c("chr1", "chr1", "chr2", "chr3"), ranges=IRanges(start=5:8, width=10),
              strand=c("+", "-", "-", "+"), gc=round(runif(4), 3))
gr
seqlens <- c(chr1=152, chr2=432, chr3=903) #specify sequence lengths to calculate gaps in downstream analysis
gr <- GRanges(seqname=c("chr1", "chr1", "chr2", "chr3"),
              ranges=IRanges(start=5:8, width=10),
              strand=c("+", "-", "-", "+"),
              gc=round(runif(4), 3),
              seqlengths=seqlens)
gr
start(gr)
end(gr)
seqnames(gr)
strand(gr)
names(gr) <- letters[1:length(gr)]
gr
table(seqnames(gr))
mcols(gr) #access metadata columns, returns a list
mcols(gr)$gc #access a specific metadata column, returns a vector
mean(mcols(gr[seqnames(gr) == "chr1"])$gc) #takes the mean for the gc content metadata of chromosome 1
mcols(gr[seqnames(gr) == "chr1"])$gc #pulls out metadata for the first chromosome in the gc content metadata column

gr1 <- GRanges(c("chr1", "chr2"), IRanges(start=c(32, 95), width=c(24, 123)))
gr2 <- GRanges(c("chr8", "chr2"), IRanges(start=c(27, 12), width=c(42, 34)))
grl <- GRangesList(gr1, gr2) #combines the granges objects into a list
grl

chrs <- c("chr3", "chr1", "chr2", "chr2", "chr3", "chr1")
gr <- GRanges(chrs, IRanges(sample(1:100, 6, replace=TRUE),
                            width=sample(3:30, 6, replace=TRUE)))
gr_split <- split(gr, seqnames(gr)) #splits the object into a list containing ranges for each chromosome
gr_split[[1]] #can access parts of the list individually
gr_split[[2]]

lapply(gr_split, function(x) order(width(x))) #can operate across an entire list
sapply(gr_split, function(x) min(start(x)))
sapply(gr_split, length)

library(BiocInstaller)
biocLite("GenomicFeatures")
biocLite("TxDb.Mmusculus.UCSC.mm10.ensGene") #download seq. data
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene
genes(txdb)
transcripts(txdb)
exons(txdb)
promoters(txdb)
mm_exons_by_gn <- exonsBy(txdb, by="gene") #separate into separate lists by gene, then order the exons in that gene
mm_exons_by_gn
