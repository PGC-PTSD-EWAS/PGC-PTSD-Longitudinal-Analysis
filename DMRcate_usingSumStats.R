########################################################################################
# DMR Analysis
########################################################################################

library(DMRcate)
library(GenomicRanges)

# Load summary statistics from meta-analysis
results <- read.csv(" ", as.is = T)

# Make data suitable for DMRcate 
annotated <- matrix(ncol=7,nrow=nrow(results))
colnames(annotated) <- c("ID", "stat", "CHR", "pos", "diff", "ind.fdr", "is.sig")
annotated <- as.data.frame(annotated)
annotated$ID <- results$IlmnID
annotated$stat <- results$z.combined
annotated$CHR <- results$CHR
annotated$pos <- results$MAPINFO
annotated$diff <- results$betaC
annotated$ind.fdr <- results$p
annotated$is.sig[annotated$ind.fdr < 0.0001] <- TRUE #change cutoff here 
annotated$is.sig[annotated$ind.fdr > 0.0001] <- FALSE
annotated$CHR <- as.factor(paste0("chr",annotated$CHR))
rownames(annotated) <- annotated$ID
dat <- makeGRangesFromDataFrame(annotated,
                                keep.extra.columns=TRUE,
                                ignore.strand=TRUE,
                                seqinfo=NULL,
                                seqnames.field=c("CHR"),
                                start.field="pos",
                                end.field=c("pos"),
                                starts.in.df.are.0based=FALSE)

setClass("ranges",slots = "ranges")
dat2 <- new("ranges",ranges = dat)
class(dat2) <- "CpGannotated"

test <- dmrcate(dat2, lambda = 1000, C = 2)
results.ranges <- extractRanges(test, genome = "hg19")
results.ranges.tb <- as.data.frame(results.ranges)
write.csv(results.ranges.tb, file = " ", row.names = F)

