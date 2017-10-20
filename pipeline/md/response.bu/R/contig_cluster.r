norm.conditions=function(ifn, ofn)
{
    tt = load.table(ifn)
    contigs = tt$contig
    m = as.matrix(tt[,-1])
    # m = t(t(m)/colSums(m))
    m = log10((m[,-1]+1) / (m[,1]+1))
    r = data.frame(contig=contigs, as.data.frame(m))
    save.table(r, ofn)
}

fold.change=function(ifn, max.change, min.coverage, ofn)
{
    tt = load.table(ifn)
    max.change.log = log10(max.change)

    # remove contigs under threshold
    total = rowSums(as.matrix(tt[,-1]))
    tt = tt[total > min.coverage,]
    cat(sprintf("contigs with enough supporting reads: %d\n", dim(tt)[1]))

    contigs = tt$contig
    m = as.matrix(tt[,-1])

    # normalize total reads per condition
    m = t(t(m) / colSums(m))

    result = NULL
    for (i in 2:dim(m)[2]) {
        x = ifelse(m[,i] == 0 & m[,i-1] == 0, 0, log10(m[,i] / m[,i-1]))
        x[x > max.change.log] = max.change.log
        x[x< -max.change.log] = -max.change.log
        result = cbind(result, round(x,5))
    }
    colnames(result) = names(tt)[-(1:2)]

    result = data.frame(contig=contigs, result)
    save.table(result, ofn)
}

log.norm=function(ifn, min.coverage, ofn)
{
    tt = load.table(ifn)

    # remove contigs under threshold
    total = rowSums(as.matrix(tt[,-1]))
    tt = tt[total > min.coverage,]
    cat(sprintf("contigs with enough supporting reads: %d\n", dim(tt)[1]))

    contigs = tt$contig
    m = as.matrix(tt[,-1])

    # normalize total reads per condition
    m = t(t(m) / colSums(m))

    result = NULL
    for (i in 2:dim(m)[2]) {
        x = ifelse(m[,i] == 0 & m[,i-1] == 0, 0, log10(m[,i] / m[,i-1]))
        x[x > max.change.log] = max.change.log
        x[x< -max.change.log] = -max.change.log
        result = cbind(result, round(x,5))
    }
    colnames(result) = names(tt)[-(1:2)]

    result = data.frame(contig=contigs, result)
    save.table(result, ofn)
}

cluster.order=function(ifn, ofn)
{
    table = load.table(ifn)
    m = as.matrix(table[,-(1:2)])
    names = table$centroid
    hh = hclust(dist(m))

    df = data.frame(centroid=names[hh$order])
    save.table(df, ofn)
}

# plots using lines and confidence interval, does not scale up
plot.centroids.lines=function(ifn.mean, ifn.sd, ifn.order, fdir)
{
    oo = load.table(ifn.order)
    centroids = oo$centroid
    N = length(centroids)

    means = load.table(ifn.mean)
    sds = load.table(ifn.sd)

    M = dim(means)[2] - 2
    xx = 1:M
    ylim = 1.1 * range(means[,-(1:2)] + sds[,-(1:2)])
    xlim = range(xx)

    size.max = max(means$size)

    fig.start(fdir=fdir, ofn=paste(fdir, "/patterns.png", sep=""), width=800, height=(200 + N*80))

    layout(matrix(1:(2*(N+1)), N+1, 2, byrow=T), widths=c(2,1), heights=c(rep(1, N), 0.4))
    for (centroid in centroids) {
        par(mai=c(0.05, 1, 0.01, 0.5))

        vv = means[means$centroid == centroid,-(1:2)]
        vv.sd = sds[sds$centroid == centroid,-(1:2)]

        # pattern
        plot.init(xlim=xlim, ylim=ylim, xlab="", ylab=centroid, add.box=T, x.axis=F, y.axis=T, add.grid=T)
        abline(h=0, lty=2)
        polygon(c(xx, rev(xx)), c(vv+vv.sd, rev(vv-vv.sd)), col="gray", border=NA)
        lines(xx, vv, lwd=2)

        # size
        par(mai=c(0.05, 0.1, 0.01, 0.5))
        vv = means[means$centroid == centroid,"size"]
        plot.init(xlim=c(0, 1.1 * size.max), ylim=c(0, 1), xlab="", ylab="", add.box=F, x.axis=F, y.axis=F, add.grid=T)
        rect(xleft=0, ybottom=0.1, xright=vv, ytop=0.9, col="gray", border=NA)
    }

    par(mai=c(0.3, 1, 0.01, 0.5))
    plot.init(xlim=xlim, ylim=ylim, xlab="", ylab="", add.box=F, x.axis=F, y.axis=F, add.grid=F)
    # axis(1)
    title(xlab="time")

    par(mai=c(0.3, 0.1, 0.01, 0.5))
    plot.init(xlim=c(0, 1.1 * size.max), ylim=c(0, 1), xlab="", ylab="", add.box=F, x.axis=F, y.axis=F, add.grid=F)
    axis(1)
    title(xlab="#contigs")

    fig.end()
}

# plots using a color scale
plot.centroids.colors=function(ifn.mean, ifn.order, ifn.contigs, ifn.contig.table, fdir)
{
    contigs = load.table(ifn.contigs)
    contig.table = load.table(ifn.contig.table)

    oo = load.table(ifn.order)
    centroids = oo$centroid
    N = length(centroids)

    means = load.table(ifn.mean)

    M = dim(means)[2] - 2
    xx = 1:M
    ylim = c(0,1)
    xlim = c(0, M)

    colors = c("black", "blue", "white", "red", "orange")
    breaks = c(-2, -1, 0, 1, 2)
    colors = c("blue", "white", "red")
    breaks = c(-2, 0, 2)
    panel = make.color.panel(colors)

    # generate total length table
    sizes = NULL
    for (centroid in centroids) {
        cen.contigs = contigs$contig[contigs$cluster == centroid]
        sizes = rbind(sizes, data.frame(centroid=centroid, length=sum(contig.table$length[is.element(contig.table$contig, cen.contigs)])))
    }
    size.range = c(10000, max(expand.range(sizes$length, 0.8, 1.2)))

    fig.start(fdir=fdir, ofn=paste(fdir, "/patterns.png", sep=""), width=600, height=(200 + N*10))

    layout(matrix(1:(2*(N+1)), N+1, 2, byrow=T), widths=c(2,1), heights=c(rep(1,N), 3))
    for (centroid in centroids) {
        cen.contigs = contigs$contig[contigs$cluster == centroid]
        par(mai=c(0.01, 1, 0.01, 0.5))

        vv = means[means$centroid == centroid,-(1:2)]
        cen.contigs = contigs$contig[contigs$cluster == centroid]

        # pattern
        plot.init(xlim=xlim, ylim=ylim, xlab="", ylab="", add.box=F, x.axis=F, y.axis=F, add.grid=F)
        mtext(centroid, 2, line=3, las=2)
        vcols = panel[vals.to.cols(t(vv), breaks, ncols=256)]
        rect(xleft=xx-1, ybottom=0, xright=xx, ytop=1, col=vcols, border=NA)

        # size
        par(mai=c(0.01, 0.1, 0.01, 0.5))
        vv = sizes$length[sizes$centroid == centroid]
        plot.init(xlim=log10(size.range), ylim=c(0, 1), xlab="", ylab="", add.box=F, x.axis=F, y.axis=F,
                  add.grid=T, grid.ny=NA)
        rect(xleft=0, ybottom=0.1, xright=log10(vv), ytop=0.9, col="gray", border=NA)
    }

    plot.new()

    par(mai=c(0.3, 0.1, 0.01, 0.5))
    plot.init(xlim=log10(size.range), ylim=c(0, 1), xlab="", ylab="", add.box=F, x.axis=F, y.axis=F, add.grid=T, grid.ny=NA)
    axis(1)

    fig.end()
}
