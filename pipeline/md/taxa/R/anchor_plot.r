plot.matrix=function(smat, field, anchors, ids, add.text.x=F, add.text.y=T, default.value, mai, breaks, colors, rotated, matrix.border)
{
    N = length(anchors)
    dim = c(N, N)

    ncols = 256
    lim = c(0,N)

    smat$set1.i = match(smat$set1, anchors)
    smat$set2.i = match(smat$set2, anchors)
    smat = smat[!is.na(smat$set1.i) & !is.na(smat$set2.i),]

    m = smatrix2matrix(smat, dim=dim, i.field="set1.i", j.field="set2.i", value.field=field, default.value=default.value)
    m = (m + t(m)) / 2
    sm = matrix2smatrix(m)
    sm[sm$i == sm$j,"value"] = 100

    # colors
    if (length(breaks) != length(colors))
        stop("#breaks != #colors")
    panel = make.color.panel(colors, ncols)
    sm$col = panel[vals.to.cols(sm$value, breaks, ncols)]

    mai[1] = 0

    par(mai=mai)
    plot.new()

    rotate=function(x,y) { list(x=(x+y)/2, y=(y-x)/2) }
    if (!rotated) {
        plot.window(xlim=lim, ylim=lim)
        rect(sm$j-1, sm$i-1, sm$j, sm$i, col=sm$col, bordersm$col)
        if (add.text.y)
            mtext(text=ids, side=2, at=1:N-0.5, las=2, line=-1, cex=0.8)
        if (add.text.x)
            mtext(text=ids, side=1, at=1:N-0.5, las=2, line=0, cex=0.8)
    } else {
        xlim = lim
        ylim = xlim/2
        p1 = rotate(x=sm$j-1, y=sm$i-1)
        p2 = rotate(x=sm$j, y=sm$i-1)
        p3 = rotate(x=sm$j, y=sm$i)
        p4 = rotate(x=sm$j-1, y=sm$i)
        plot.window(xlim=xlim, ylim=ylim, yaxs="i")
        x = as.vector(t(cbind(p1$x, p2$x, p3$x, p4$x, NA)))
        y = as.vector(t(cbind(p1$y, p2$y, p3$y, p4$y, NA)))
        col = sm$col
        polygon(x=x, y=y, col=col, border=col)
        polygon(x=x, y=y, col=NA, border=matrix.border)
    }
}

plot.box=function(x, anchors, ids, col="blue", ylim, add.axis=T, do.log=F, mai, add.grid=F, label, factor=1, plot.top95.line=T, axis.side=2)
{
    x$coord = match(x$anchor, anchors) - 0.5
    N = length(anchors)
    coords = 0:N

    fields = c("bottom05", "bottom25", "median", "top75", "top95")
    for (field in fields)
        x[field] = x[,field] * factor

    xlim = c(0,N)
    ylim = range(c(x$bottom05, 1.1*x$top95))

    log = ifelse(do.log, "y", "")

    bg.cols = c("white", colors()[246])
    par(mai=mai)
    plot.new()
    plot.window(xlim=xlim, ylim=ylim, log=log)
    # title(ylab=label)
    mtext(label, side=axis.side, outer=F, col=1, line=4, las=2)
    rect(xleft=coords, ybottom=rep(ylim[1],N+1), xright=coords+1, ytop=rep(ylim[2], N+1), col=bg.cols, border=bg.cols)
    if (add.grid) grid(col="darkgray", lty=1)
    if (plot.top95.line)  segments(x0=x$coord, x1=x$coord, y0=x$bottom05, y1=x$top95)
    rect(xleft=x$coord-0.25, ybottom=x$bottom25, xright=x$coord+0.25, ytop=x$top75, col=col, border=col)
    segments(x0=x$coord-0.25, x1=x$coord+0.25, y0=x$median, y1=x$median, lwd=2)
    segments(x0=0, x1=N, y0=ylim[1], y1=ylim[1])
    if (add.axis)
        axis(axis.side, las=2, cex.axis=0.75)
}

plot.bar=function(x, field, anchors, ids, col="blue", ylim, add.axis=T, axis.at=NULL, axis.labels=NULL, mai, label, ybottom=0,
    ytop=NULL, axis.side=2)
{
    x$coord = match(x$anchor, anchors) - 0.5
    N = length(anchors)
    coords = 0:N

    xlim = c(0,N)
    ytop = ifelse(is.null(ytop), 1.1 * max(x[,field]), ytop)
    ylim = c(ybottom, ytop)
    bg.cols = c("white", colors()[246])
    par(mai=mai)
    plot.new()
    plot.window(xlim=xlim, ylim=ylim)
    mtext(label, side=axis.side, outer=F, col=1, line=4, las=2)
    rect(xleft=coords, ybottom=rep(ylim[1],N+1), xright=coords+1, ytop=rep(ylim[2], N+1), col=bg.cols, border=bg.cols)
    rect(xleft=x$coord-0.25, ybottom=ybottom, xright=x$coord+0.25, ytop=x[,field], col=col, border=NA)
    segments(x0=0, x1=N, y0=ybottom, y1=ybottom)
    if (add.axis) {
        if (!is.null(axis.at))
            axis(axis.side, at=axis.at, labels=axis.labels, las=2, cex.axis=0.75)
        else
            axis(axis.side, las=2, cex.axis=0.75)
    }
}

plot.bar.append=function(x, field, anchors, col, ybottom=0)
{
    x$coord = match(x$anchor, anchors) - 0.5
    N = length(anchors)
    coords = 0:N
    rect(xleft=x$coord-0.25, ybottom=ybottom, xright=x$coord+0.25, ytop=x[,field], col=col, border=NA)
}

plot.names=function(anchors, ids, mai, labels=ids, srt=90, adj=1, cex=1.1, ycoord=1)
{
    N = length(anchors)
    coords = 1:N

    xlim = c(0, N)
    ylim = c(0, 1)

    mai[1] = 0
    mai[3] = 0

    par(mai=mai)
    plot.new()
    plot.window(xlim=xlim, ylim=ylim)
    text(x=coords - 0.5,  y=ycoord, labels=labels, srt=srt, cex=cex, adj=adj)
}

plot.taxa=function(taxa, anchors, ids, mai, axis.side)
{
    ix = match(anchors, taxa$anchor)
    col = taxa$color[ix]
    is.cag = taxa$is.cag[ix]
    N = length(anchors)
    coords = 1:N - 1

    xlim = c(0, N)
    ylim = c(0, 1)

    mai[1] = 0.05
    mai[3] = 0

    par(mai=mai)
    plot.new()
    plot.window(xlim=xlim, ylim=ylim)
    mtext("taxa", side=axis.side, outer=F, col=1, line=4, las=2)

    n.cag = sum(is.cag)
    cag.line.height = 0.83

    par(xpd=F)
    # rect(xleft=coords[is.cag], xright=coords[is.cag]+1, ybottom=0, ytop=1, border=NA, col=col[is.cag], density=25, angle=45, lwd=1)
    # rect(xleft=coords[!is.cag], xright=coords[!is.cag]+1, ybottom=0, ytop=1, border=NA, col=col[!is.cag])
    # rect(xleft=coords, xright=coords+1, ybottom=0, ytop=1, border=1, col=NA)
    rect(xleft=coords, xright=coords+1, ybottom=0, ytop=1, border=col, col=col)
    text(x=coords+0.5, y=0.5, taxa$letter[ix], cex=1.2, adj=c(0.5, 0.5))
    segments(x0=coords[is.cag]+0.25, x1=coords[is.cag]+0.75, y0=rep(cag.line.height, n.cag), y1=rep(cag.line.height, n.cag), lwd=1)
}

plot.group=function(groups, anchors, ids, mai)
{
    N = length(anchors)

    ix = match(anchors, groups$anchor)
    colors = rainbow(length(unique(groups$group)))
    groups$coord = 1:N - 1
    s = sapply(split(groups$coord, groups$group), range)

    # col = groups$color[ix]

    # N = length(anchors)
    # coords = 1:N - 1

    xlim = c(0, N)
    ylim = c(0, 1)

    mai[1] = 0.05
    mai[3] = 0

    par(mai=mai)
    plot.new()
    plot.window(xlim=xlim, ylim=ylim)
    mtext("group", side=2, outer=F, col=1, line=4, las=2)

    for (i in 1:dim(s)[2]) {
        group = as.numeric(colnames(s)[i])
        col = colors[group]
        start = s[1,i]
        end = s[2,i]+1
        xcoords = c(start, end)
        middle = (start + end) / 2
        # rect(xleft=start, xright=end, ybottom=0, ytop=1, border=1, col=col)
        segments(x0=xcoords, x1=xcoords, y0=0, y1=1)
        text(x=middle, y=0.5, labels=group, cex=0.85, srt=90)
    }
    # rect(xleft=coords, xright=coords+1, ybottom=0, ytop=1, border=col, col=col)
    # text(x=coords+0.5, y=0.5, groups$group[ix], cex=0.85, srt=90)
}

compute.quantiles=function(x, field, anchors)
{
    result = NULL
    for (anchor in anchors) {
        xx = x[x$anchor == anchor,]
        values = xx[,field]
        df = data.frame(
            anchor=anchor,
            count=length(values),
            bottom05=quantile(values, 0.05),
            bottom25=quantile(values, 0.25),
            median=median(values),
            top75=quantile(values, 0.75),
            top95=quantile(values, 0.95))
        result = rbind(result, df)
    }
    result
}

plot.complete.basic=function(contig.ifn, ca.ifn, coverage.ifn, gc.ifn, mat.ifn, order.ifn, info.ifn, group.ifn, taxa.ifn, checkm.ifn, min.identity, fdir)
{
    contigs = load.table(contig.ifn)
    ca = load.table(ca.ifn)
    coverage = load.table(coverage.ifn)
    gc = load.table(gc.ifn)
    smat = load.table(mat.ifn)
    taxa = load.table(taxa.ifn)
    groups = load.table(group.ifn)
    checkm = load.table(checkm.ifn)
    checkm$anchor = checkm$Bin.Id

    # reps$gene.divergence = log10(1+100-reps$anchor.coverage)
    # reps$local.divergence = log10(1+100-reps$anchor.identity)
    div.axis.labels = c(0, 2, 10, 50)
    div.axis.at = log10(1 + div.axis.labels)

    ca$contig.length.k = contigs$length[match(ca$contig, contigs$contig)] / 1000
    s = sapply(split(ca$contig.length.k, ca$anchor), sum)
    genome.size = data.frame(anchor=names(s), mb=s/1000)

    info = load.table(info.ifn)
    info$anchor.size.m = info$anchor.size/10^6
    info$union.size.m = info$union.size/10^6

    anchor.table = load.table(order.ifn)
    anchors = anchor.table$set
    ids = anchor.table$id

    # enrichment summary
    enrichment = compute.quantiles(x=ca, field="enrichment", anchors=anchors)

    # contig length
    contig.length = compute.quantiles(x=ca, field="contig.length.k", anchors=anchors)

    # identity
    identity.breaks = c(30, 50, 60, 80, 100)
    identity.colors = c("darkblue", "blue", "red", "orange", "yellow")
    identity.panel = make.color.panel(identity.colors)
    wlegend2(fdir=fdir, panel=identity.panel, breaks=identity.breaks, title="identity")

    # shared genes
    shared.breaks = c(0, 1, 10, 100, 1000)
    shared.colors = c("white", "blue", "blue", "red", "orange")
    shared.panel = make.color.panel(shared.colors)
    wlegend2(fdir=fdir, panel=shared.panel, breaks=shared.breaks, title="shared")

    ll = list(
        list(field="identity", default=min.identity, colors=identity.colors, breaks=identity.breaks, matrix.border=NA),
        list(field="shared", default=0, colors=shared.colors, breaks=shared.breaks, matrix.border="lightgray"))

    color = colors()[139]
    axis.side = 4

    for (i in 1:2) {
        llx = ll[[i]]
        field = llx$field

        fig.start(fdir=fdir, ofn=paste(fdir, "/complete_", field, ".pdf", sep=""), type="pdf", height=7, width=10)

        mai = c(0.1,0.1,0.01,2)
        par(xaxs="i")

        N.small = 2
        N.big = 5
        N = 1 + N.big + N.small
        layout(matrix(1:N, N, 1), heights=c(6, 0.4, rep(1,N.big), 0.4))

        # matrix
        field = "identity"
        default.value = min.identity

        plot.matrix(smat=smat, field=llx$field, anchors=anchors, ids=ids, add.text.y=T, mai=mai,
                    breaks=llx$breaks, colors=llx$colors, default.value=default.value, rotated=T, matrix.border=llx$matrix.border)

        # small
        # plot.group(groups=groups, anchors=anchors, ids=ids, mai=mai)

        # small
        plot.names(anchors=anchors, ids=ids, mai=mai)

        # big
        # plot.box(x=gc, anchors=anchors, ids=ids, do.log=F, mai=mai, col=colors()[640], label="GC")
        plot.box(x=gc, anchors=anchors, ids=ids, do.log=F, mai=mai, col=color, label="GC", axis.side=axis.side)

        # big
        # plot.box(x=coverage, anchors=anchors, ids=ids, do.log=F, mai=mai, col=colors()[121], label="abundance")
        plot.box(x=coverage, anchors=anchors, ids=ids, do.log=F, mai=mai, col=color, label="abundance", axis.side=axis.side)

        # big
        # plot.bar(x=checkm, field="Completeness", anchors=anchors, ids=ids, mai=mai, col=colors()[139], label="completeness (%)")
        plot.bar(x=checkm, field="Completeness", anchors=anchors, ids=ids, mai=mai, col=color, label="completeness (%)", axis.side=axis.side)
        plot.bar(x=checkm, field="Contamination", anchors=anchors, ids=ids, mai=mai, col=color, label="multi-strain (%)", axis.side=axis.side)

        # big
        # plot.bar(x=info, field="length.m", anchors=anchors, ids=ids, mai=mai, col=colors()[640], label="size (Mb)")
        plot.bar(x=info, field="union.size.m", anchors=anchors, ids=ids, mai=mai, col="red", label="size (Mb)", axis.side=axis.side)
        plot.bar.append(x=info, field="anchor.size.m", anchors=anchors, col=color, ybottom=0)

        # small
        plot.names(anchors=anchors, ids=ids, mai=mai)

        fig.end()
    }
}

plot.complete.taxa=function(order.ifn, taxa.ifn, res.ifn, reps.ifn, fdir)
{
    taxa = load.table(taxa.ifn)
    res = load.table(res.ifn)
    reps = load.table(reps.ifn)
    anchor.table = load.table(order.ifn)
    anchors = anchor.table$set
    ids = anchor.table$id
    resolve = res$short.level[match(ids, res$anchor.id)]
    color = colors()[139]

    fig.start(fdir=fdir, ofn=paste(fdir, "/complete_taxa.pdf", sep=""), type="pdf", height=2, width=10)

    mai = c(0.1,0.1,0.01,2)
    par(xaxs="i")

    N.small = 3
    N.big = 2
    N = N.big + N.small
    layout(matrix(1:N, N, 1), heights=c( c(0.2, 0.2), rep(0.4,N.big), 0.2))

    color = colors()[139]
    axis.side = 4

    # small
    plot.names(anchors=anchors, ids=ids, labels=resolve, mai=mai, srt=90, cex=0.9, ycoord=0.1, adj=0)

    # small
    plot.taxa(taxa=taxa, anchors=anchors, ids=ids, mai=mai, axis.side=axis.side)

    # big
    plot.bar(x=res, field="coverage", anchors=anchors, ids=ids, mai=mai, col=color, label="Coverage", axis.side=axis.side)

    # big
    plot.bar(x=res, field="identity", anchors=anchors, ids=ids, mai=mai, col=color, label="AAI",
             ybottom=40, ytop=100, axis.side=axis.side)

    # small
    plot.names(anchors=anchors, ids=ids, mai=mai)

    fig.end()
}

