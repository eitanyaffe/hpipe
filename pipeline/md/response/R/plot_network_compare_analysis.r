
plot.abundance.summary=function(ifn.elements, ifn.anchors, ifn.clusters, ifn.map, min.contacts, fdir)
{
    element.table = load.table(ifn.elements)
    map = load.table(ifn.map)
    anchor.table = load.table(ifn.anchors)
    cluster.table = load.table(ifn.clusters)

    element.table$abundance1 = cluster.table$abundance1[match(element.table$cluster, cluster.table$cluster)]
    element.table$abundance2 = cluster.table$abundance2[match(element.table$cluster, cluster.table$cluster)]

    height = 6
    width = 6
    plot.scatter.anchors=function(add.labels) {
        lim = range(c(anchor.table$abundance1, anchor.table$abundance2))
        cc = round(cor(anchor.table$abundance1, anchor.table$abundance2),2)
        plot.init(xlim=lim, ylim=lim, main=paste("anchors, rho:", cc), xlab="cond1", ylab="cond2")
        grid()
        abline(a=0, b=1, col="gray")
        points(anchor.table$abundance1, anchor.table$abundance2, pch=19, col="orange")
        if (add.labels)
            text(anchor.table$abundance1, anchor.table$abundance2, anchor.table$anchor.id, pos=4)
    }
    fig.start(fdir=fdir, ofn=paste(fdir, "/anchor_scatter_clean.pdf", sep=""), type="pdf", height=height, width=width)
    plot.scatter.anchors(add.labels=F)
    fig.end()
    fig.start(fdir=fdir, ofn=paste(fdir, "/anchor_scatter_labels.pdf", sep=""), type="pdf", height=height, width=width)
    plot.scatter.anchors(add.labels=T)
    fig.end()

    plot.scatter.elements=function(add.labels) {
        lim = range(c(element.table$abundance1, element.table$abundance2))
        cc = round(cor(element.table$abundance1, element.table$abundance2),2)
        plot.init(xlim=lim, ylim=lim, main=paste("elements, rho:", cc), xlab="cond1", ylab="cond2")
        grid()
        abline(a=0, b=1, col="gray")
        points(element.table$abundance1, element.table$abundance2, pch=19, col="orange")
        if (add.labels)
            text(element.table$abundance1, element.table$abundance2, element.table$id, pos=4)
    }
    fig.start(fdir=fdir, ofn=paste(fdir, "/element_scatter_clean.pdf", sep=""), type="pdf", height=height, width=width)
    plot.scatter.elements(add.labels=F)
    fig.end()
    fig.start(fdir=fdir, ofn=paste(fdir, "/element_scatter_labels.pdf", sep=""), type="pdf", height=height, width=width)
    plot.scatter.elements(add.labels=T)
    fig.end()

}

plot.scatter.summary=function(ifn.elements, ifn.anchors, ifn.clusters, ifn.map, min.contacts, fdir)
{
    element.table = load.table(ifn.elements)
    map = load.table(ifn.map)
    anchor.table = load.table(ifn.anchors)
    cluster.table = load.table(ifn.clusters)

    # post abundant elements
    clusters = cluster.table$cluster

    # post abundant anchors
    anchors = anchor.table$anchor
    map = map[is.element(map$cluster, clusters) & is.element(map$anchor, anchors),]
    map$s1 = map$score1
    map$s2 = map$score2
    map$no.observed = map$observed1 == 0 | map$observed2 == 0
    map.high = map[map$observed1 > min.contacts | map$observed2 > min.contacts,]
    height = 6
    width = 5

    # observed scatter
    fig.start(fdir=fdir, ofn=paste(fdir, "/observed.pdf", sep=""), type="pdf", height=height, width=width)
    mm = median(map.high$observed2 / map.high$observed1)
    cc = cor(map.high$observed2, map.high$observed1)
    olim = range(c(0, log10(1+map$observed1), log10(1+map$observed2)))
    plot.init(xlim=olim, ylim=olim, xlab="log10(o1)", ylab="log10(o2)", main=sprintf("median(o2/o1)=%.2f, rho=%.2f", mm, cc))
    grid()
    points(log10(1+map$observed1), log10(1+map$observed2), pch=".", cex=2, col=ifelse(!map$no.observed, 1, 2))
    fig.end()

    # enrichment scatter
    fig.start(fdir=fdir, ofn=paste(fdir, "/enrichment.pdf", sep=""), type="pdf", height=height, width=width)
    mm = median(map.high$s2 - map.high$s1)
    cc = cor(map.high$s2, map.high$s1)
    clim = range(c(0, map.high$s1, map.high$s2))
    plot.init(xlim=clim, ylim=clim, xlab="log10(o1/e1)", ylab="log10(o2/e2)", main=sprintf("median(c2-c1)=%.2f, rho=%.2f", mm, cc))
    grid()
    points(map.high$s1, map.high$s2, pch=".", cex=2, col=ifelse(!map.high$no.observed, 1, 2))
    fig.end()
}


