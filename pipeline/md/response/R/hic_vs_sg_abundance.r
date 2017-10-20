plot.scatter=function(ifn.ca, ifn.anchors, ifn.hic.coverage, ifn.sg.coverage, ifn.taxa.legend, fdir)
{
    ca = load.table(ifn.ca)
    anchor.table = load.table(ifn.anchors)
    hic = load.table(ifn.hic.coverage)
    sg = load.table(ifn.sg.coverage)
    tlegend = load.table(ifn.taxa.legend)

    anchors = anchor.table$set

    m = merge(hic, sg, by="contig")
    df = data.frame(contig=m$contig, hic=m$abundance.enrichment.x, sg=m$abundance.enrichment.y)
    df$anchor = ca$contig_anchor[match(df$contig, ca$contig)]
    df = df[!is.na(df$anchor) & df$anchor != 0,]
    get.anchor.abundance=function(field) {
        s = sapply(split(df[,field], df$anchor), median)
        s[match(anchors,names(s))]
    }
    anchor.table$abundance.hic = get.anchor.abundance("hic")
    anchor.table$abundance.sg = get.anchor.abundance("sg")
    anchor.table$color = tlegend$color[match(anchor.table$id, tlegend$anchor.id)]

    height = 6
    width = 6
    plot.scatter.anchors=function(add.labels) {
        lim = range(c(anchor.table$abundance.hic, anchor.table$abundance.sg))
        cc = round(cor(anchor.table$abundance.hic, anchor.table$abundance.sg),2)
        plot.init(xlim=lim, ylim=lim, main=paste("anchors, rho:", cc), xlab="Hi-C", ylab="Shotgun")
        grid()
        abline(a=0, b=1, col="gray")
        points(anchor.table$abundance.hic, anchor.table$abundance.sg, pch=19, col=anchor.table$color)
        if (add.labels)
            text(anchor.table$abundance.hic, anchor.table$abundance.sg, anchor.table$id, pos=4)
    }
    fig.start(fdir=fdir, ofn=paste(fdir, "/anchor_scatter_clean.pdf", sep=""), type="pdf", height=height, width=width)
    plot.scatter.anchors(add.labels=F)
    fig.end()
    fig.start(fdir=fdir, ofn=paste(fdir, "/anchor_scatter_labels.pdf", sep=""), type="pdf", height=height, width=width)
    plot.scatter.anchors(add.labels=T)
    fig.end()
}
