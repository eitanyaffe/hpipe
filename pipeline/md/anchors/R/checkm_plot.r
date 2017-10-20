plot.analysis=function(ifn.anchors, ifn.qa, ifn.ca, ifn.contigs, ifn.info, type, fdir)
{
    df = load.table(ifn.anchors)
    qa = load.table(ifn.qa)
    info = load.table(ifn.info)

    df = df[is.element(df$set, qa$Bin.Id),]

    ca = load.table(ifn.ca)
    fc = field.count(ca, field="contig")
    ca$is.multi = fc$count[match(ca$contig,fc$contig)]>1

    contigs = load.table(ifn.contigs)
    ca$length = contigs$length[match(ca$contig, contigs$contig)]
    ids = df$id
    N = length(ids)

    if (type == "A")
        ca = ca[ca$anchor == ca$contig_anchor,]
    if (type == "X")
        ca = ca[ca$contig_anchor == 0,]
    if (type == "nS")
        ca = ca[!ca$is.multi,]
    if (type == "S")
        ca = ca[ca$is.multi,]

    qa$anchor.id = df$id[match(qa$Bin.Id, df$set)]
    qa$coverage = info$coverage[match(qa$Bin.Id, info$anchor)]

    sx = sapply(split(ca$length, ca$anchor), median)
    qa$median.contig.length = sx[match(qa$Bin.Id, names(sx))] / 1000

    sx = sapply(split(ca$length, ca$anchor), sum)
    qa$genome.size = sx[match(qa$Bin.Id, names(sx))] / 1000

    qa = qa[order(match(qa$anchor.id, ids)),]

    width = 1+N*0.15

    #########################################################################################################

    fig.start(fdir=fdir, ofn=paste(fdir, "/complete.pdf", sep=""), type="pdf", width=width, height=3)
    barplot(qa$Completeness, names.arg=ids, border=NA, col="darkgreen", ylim=c(0,100), las=2, cex.names=0.8)
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/complete_sorted.pdf", sep=""), type="pdf", width=width, height=3)
    ix = order(qa$Completeness)
    barplot(qa$Completeness[ix], names.arg=ids[ix], border=NA, col="darkgreen", ylim=c(0,100), las=2, cex.names=0.8)
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/contamination.pdf", sep=""), type="pdf", width=width, height=3)
    barplot(qa$Contamination, names.arg=ids, border=NA, col="red", ylim=c(0,100), las=2, cex.names=0.8)
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/hetero.pdf", sep=""), type="pdf", width=width, height=3)
    barplot(qa$Strain.heterogeneity, names.arg=ids, border=NA, col="blue", ylim=c(0,100), las=2, cex.names=0.8)
    fig.end()

    #########################################################################################################
    # coverage length

    xlim = range(qa$coverage)
    ylim = range(qa$Completeness)
    cc = round(cor(qa$coverage, qa$Completeness, method="spearman"), 2)
    fig.start(fdir=fdir, ofn=paste(fdir, "/complete_vs_coverage.pdf", sep=""), type="pdf", width=5, height=5)
    plot.init(xlim=xlim, ylim=ylim, xlab="coverage", ylab="complete", log="x", main=paste("spearman=", cc, sep=""))
    points(qa$coverage, qa$Completeness, pch="+")
    fig.end()

    #########################################################################################################
    # median contig length

    xlim = c(min(qa$median.contig.length), 2*max(qa$median.contig.length))
    ylim = range(qa$Completeness)
    cc = round(cor(qa$median.contig.length, qa$Completeness, method="spearman"), 2)
    fig.start(fdir=fdir, ofn=paste(fdir, "/complete_vs_contig_length.pdf", sep=""), type="pdf", width=5, height=5)
    plot.init(xlim=xlim, ylim=ylim, xlab="kb", ylab="complete", log="x", main=paste("spearman=", cc, sep=""))
    points(qa$median.contig.length, qa$Completeness, pch="+")
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/complete_vs_contig_length_label.pdf", sep=""), type="pdf", width=5, height=5)
    plot.init(xlim=xlim, ylim=ylim, xlab="kb", ylab="complete", log="x", main=paste("spearman=", cc, sep=""))
    points(qa$median.contig.length, qa$Completeness, pch="+")
    text(qa$median.contig.length, qa$Completeness, pos=4, qa$anchor.id, cex=0.5)
    fig.end()

    #########################################################################################################
    # genome size

    xlim = c(min(qa$genome.size), max(qa$genome.size))
    ylim = range(qa$Completeness)
    cc = round(cor(qa$genome.size, qa$Completeness, method="spearman"), 2)
    fig.start(fdir=fdir, ofn=paste(fdir, "/complete_vs_genome_size.pdf", sep=""), type="pdf", width=5, height=5)
    plot.init(xlim=xlim, ylim=ylim, xlab="kb", ylab="complete", main=paste("spearman=", cc, sep=""))
    points(qa$genome.size, qa$Completeness, pch="+")
    fig.end()

    xlim = c(min(qa$genome.size), max(qa$genome.size))
    ylim = range(qa$Completeness)
    cc = round(cor(qa$genome.size, qa$Completeness, method="spearman"), 2)
    fig.start(fdir=fdir, ofn=paste(fdir, "/complete_vs_genome_size.pdf", sep=""), type="pdf", width=5, height=5)
    plot.init(xlim=xlim, ylim=ylim, xlab="kb", ylab="complete", main=paste("spearman=", cc, sep=""))
    points(qa$genome.size, qa$Completeness, pch="+")
    text(qa$genome.size, qa$Completeness, pos=match(qa$anchor.id,ids) %% 4 + 1, qa$anchor.id, cex=0.5)
    fig.end()

    #########################################################################################################
    # contamination vs contig size

    xlim = c(min(qa$median.contig.length), 2*max(qa$median.contig.length))
    ylim = range(qa$Contamination)
    cc = round(cor(qa$median.contig.length, qa$Contamination, method="spearman"), 2)
    fig.start(fdir=fdir, ofn=paste(fdir, "/contamination_vs_contig_length.pdf", sep=""), type="pdf", width=5, height=5)
    plot.init(xlim=xlim, ylim=ylim, xlab="kb", ylab="contamination", log="x", main=paste("spearman=", cc, sep=""))
    points(qa$median.contig.length, qa$Contamination, pch="+")
    fig.end()

    #########################################################################################################
    # complete hist

    fig.start(fdir=fdir, ofn=paste(fdir, "/complete_hist.pdf", sep=""), type="pdf", width=5, height=5)
    hist(qa$Completeness, col="gray", main="completeness", breaks=10, xlim=c(0,100))
    fig.end()

    #########################################################################################################
    # stats

    fc = file(paste(fdir, "/stats.txt", sep=""))
    lines = c(
        sprintf("median completeness: %.1f", median(qa$Completeness)),
        sprintf("mean completeness: %.1f", mean(qa$Completeness)),
        sprintf("number of unions over 20%% completeness: %d (%.1f%%)", sum(qa$Completeness>20), 100*sum(qa$Completeness>20)/dim(qa)[1])
        )
    writeLines(lines, fc)
    close(fc)

}
