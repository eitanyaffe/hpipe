plot.abundance=function(ifn.coverage, ifn.seeds, ifn.ca, ifn.gc, fdir)
{
    table = load.table(ifn.coverage)
    gc = load.table(ifn.gc)
    ca = load.table(ifn.ca)

    df = data.frame(contig=table$contig, a=table$abundance.enrichment, l=table$length, gc=gc$gc[match(table$contig,gc$contig)])

    # compute relative mass of anchored
    df$is.anchored = is.element(df$contig, ca$contig)
    df$mass = df$l * 10^df$a
    df$mass = 100 * df$mass / sum(df$mass)
    anchored.mass = sum(df$mass[df$is.anchored])

    seeds = load.table(ifn.seeds)
    seeds = seeds[seeds$cluster != -1,]

    binsize = 10*10^6
    N = 200 + 1
    breaks = seq(min(df$a), max(df$a), length.out=N)
    delta = breaks[2] - breaks[1]

    dense=function(df) {
        s = sapply(split(df$l, cut(df$a, breaks=breaks, include.lowest=T)), sum)
        result = data.frame(ab=breaks[-N], bp=s, cs=cumsum(s))
        result$dd = (result$bp / sum(result$bp)) / delta
        result
    }
    result.all = dense(df)
    indices = c(which(c(T, diff(round(result.all$cs / binsize)) != 0)), dim(result.all)[1])

    fig.start(fdir=fdir, ofn=paste(fdir, "/all.png", sep=""), width=800)
    par(yaxs="i")
    plot.init(xlim=range(result.all$a), ylim=range(result.all$dd), xlab="abundance", ylab="dd")
    bcol = c("white", "gray")
    for (i in 2:length(indices)) {
        ab = result.all$ab[indices[i-1]:indices[i]]
        dd = result.all$dd[indices[i-1]:indices[i]]
        polygon(x=c(ab,rev(ab)), y=c(dd,rep(0,length(dd))), col=bcol[1+i%%2], border=NA)
    }
    lines(x=result.all$a, y=result.all$dd)
    fig.end()

    # residual
    # df.resid = df[!is.element(df$contig, ca$contig),]
    # result.resid = dense(df.resid)
    df.anchored = df[is.element(df$contig, ca$contig),]
    result.anchored = dense(df.anchored)
    fig.start(fdir=fdir, ofn=paste(fdir, "/anchored.png", sep=""), width=800)
    par(yaxs="i")
    plot.init(xlim=range(result.all$a), ylim=range(result.all$bp), xlab="abundance", ylab="bp", main=paste("anchored mass: ", round(anchored.mass,1), "%", sep=""))
    bcol = c("white", "gray")
    for (i in 2:length(indices)) {
        ab = result.all$ab[indices[i-1]:indices[i]]
        bp = result.all$bp[indices[i-1]:indices[i]]
        polygon(x=c(ab,rev(ab)), y=c(bp,rep(0,length(bp))), col=bcol[1+i%%2], border=NA)
    }
    lines(x=result.all$a, y=result.all$bp)
    lines(x=result.anchored$a, y=result.anchored$bp, col="red", lwd=2)
    fig.end()

    # add anchor means
    ca$abundance.enrichment = table$abundance.enrichment[match(ca$contig, table$contig)]
    ss = sapply(split(ca$abundance.enrichment, ca$anchor), mean)
    fig.start(fdir=fdir, ofn=paste(fdir, "/anchors.png", sep=""), width=800)
    plot(result.all$a, result.all$bp, col="gray", type="l", lwd=2, xlab="abundance")
    abline(v=ss, col="red", lwd=1)
    fig.end()

    # GC vs abundance

    width = 800
    height = 880

    # clean
    fig.start(fdir=fdir, ofn=paste(fdir, "/gc_abundance_all.png", sep=""), width=width, height=height)
    plot(df$a, df$gc, xlab="abundance", ylab="GC", pch=".", col=1, cex=1)
    fig.end()

    # color anchors
    ix = is.element(df$contig, ca$contig)
    fig.start(fdir=fdir, ofn=paste(fdir, "/gc_abundance_anchors.png", sep=""), width=width, height=height)
    plot(df$a, df$gc, xlab="abundance", ylab="GC", pch=".", col="gray", cex=1)
    points(df$a[ix], df$gc[ix], xlab="abundance", ylab="GC", pch=".", col="darkred", cex=2)
    fig.end()

    # color seed anchors
    ix = is.element(df$contig, seeds$contig)
    fig.start(fdir=fdir, ofn=paste(fdir, "/gc_abundance_anchor_seeds.png", sep=""), width=width, height=height)
    plot(df$a, df$gc, xlab="abundance", ylab="GC", pch=".", col="gray", cex=1)
    points(df$a[ix], df$gc[ix], xlab="abundance", ylab="GC", pch=".", col="darkred", cex=2)
    fig.end()
}
