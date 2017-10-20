
plot.anchor.map=function(ifn.majors, ifn.map, ifn.elements, ifn.order, fdir)
{
    anchors = load.table(ifn.order)$anchor
    majors = load.table(ifn.majors)
    map = load.table(ifn.map)

    map$en = log10(map$observed/map$expected)
    map$low = map$observed == 0
    map$lexp = log10(map$expected)
    elements = load.table(ifn.elements)
    majors = majors[match(anchors, majors$anchor),]

    breaks.pos = c(0,3)
    colors.pos = c("white", "red")
    panel.pos = make.color.panel(colors=colors.pos)
    wlegend2(fdir=fdir, panel=panel.pos, breaks=breaks.pos, title="element_map_enrichment")

    breaks.expected = c(0,1)
    colors.expected = c("white", "blue")
    panel.expected = make.color.panel(colors=colors.expected)
    wlegend2(fdir=fdir, panel=panel.expected, breaks=breaks.expected, title="element_map_enrichment")

    map$col.en = panel.pos[vals.to.cols(map$en, breaks.pos)]
    map$col.exp = panel.expected[vals.to.cols(map$expected, breaks.expected)]

    N = length(anchors)
    fig.start(fdir=fdir, ofn=paste(fdir, "/element_hic_map.png", sep=""), height=1200, width=1200)
    layout(matrix(1:(ceiling(N/8)*8), ceiling(N/8), 8, byrow=T))
    par(mai=c(0.4, 0.1, 0.3, 0.1))
    for (anchor in anchors) {
        major = majors$max.cluster[anchor == majors$anchor]
        minors = elements$cluster[elements$anchor == anchor & !elements$major]
        aclusters = c(major, minors)
        M = length(aclusters)
        if (M == 1)
            next
        amap = map[is.element(map$cluster1, aclusters) & is.element(map$cluster2, aclusters),]
        amap$dist = ifelse(amap$en > 0, max(amap$en) - amap$en, max(amap$en))
        amap$i1 = match(amap$cluster1, aclusters)
        amap$i2 = match(amap$cluster2, aclusters)
        hh = hclust(as.dist(smatrix2matrix(smatrix=amap, dim=c(M,M), i.field="i1", j.field="i2", value.field="dist", default.value=0)))
        oclusters = aclusters[hh$order]
        aclusters = c(major, setdiff(oclusters, major))
        amap$i1 = match(amap$cluster1, aclusters)
        amap$i2 = match(amap$cluster2, aclusters)
        amap$col = ifelse(amap$i1 > amap$i2, amap$col.en, amap$col.exp)

        plot.new()
        plot.window(xlim=c(0,M), ylim=c(0,M), xaxs="i", yaxs="i")
        title(main=anchor)

        rect(1:M-1, 1:M-1, 1:M, 1:M, col="gray", border=NA)
        rect(amap$i1-1, amap$i2-1, amap$i1, amap$i2, col=amap$col, border=NA)
        axis(1, at=1:M-0.5, labels=aclusters, las=2)
        abline(h=1, v=1)
        box()
    }
    fig.end()
}

# map should be complete with cluster1,cluster2
plot.hic.map=function(ifn.map, ifn.network, ifn.elements, ifn.anchors, min.contacts, lid, fdir)
{
    element.table = load.table(ifn.elements)
    anchor.table = load.table(ifn.anchors)

    map = load.table(ifn.map)
    map$en = map$score
    map$en[map$sum.observed < min.contacts] = 0

    plot.matrix=function(m, Nx, Ny, labels.x, labels.y, title, col.field, ofn, plot.network=F) {
        fig.start(fdir=fdir, ofn=paste(fdir, "/", ofn, ".png", sep=""), height=200+Ny*10, width=200+Nx*10)
        plot.new()
        plot.window(xlim=c(0,Nx), ylim=c(0,Ny), xaxs="i", yaxs="i")
        title(main=paste(lid, title))

        rect(m$x-1, m$y-1, m$x, m$y, col=m[,col.field], border=NA)
        axis(1, at=1:Nx-0.5, labels=labels.x, las=2)
        axis(2, at=1:Ny-0.5, labels=labels.y, las=2)
        if (plot.network) points(network$x-0.5, network$y-0.5, pch=19, cex=0.75)
        abline(v=0:Nx, col="darkgray")
        abline(h=0:Ny, col="darkgray")
        fig.end()
    }
    make.matrix=function(df, x, y) {
        m = merge(df, map)
        m$x = match(m$anchor, x)
        m$y = match(m$cluster, y)
        m$col.en = panel.en[vals.to.cols(m$en, breaks.en)]
        m$col.obs = panel.obs[vals.to.cols(m$observed, breaks.obs)]
        m
    }

#    breaks.en = c(0,1,2,3)
#    colors.en = c("white", "blue", "red", "orange")
    breaks.en = c(0,4)
    colors.en = c("white","red")
    panel.en = make.color.panel(colors=colors.en)
    wlegend2(fdir=fdir, panel=panel.en, breaks=breaks.en, title="hic_enrichment")

    breaks.obs = c(0,10,100,1000)
    colors.obs = c("white", "blue", "red", "orange")
    panel.obs = make.color.panel(colors=colors.obs)
    wlegend2(fdir=fdir, panel=panel.obs, breaks=breaks.obs, title="hic_obs")

    network = load.table(ifn.network)
    anchors = anchor.table$set
    anchor.ids = anchor.table$id
    element.ids = element.table$id
    elements = element.table$cluster
    # elements = sort(unique(network$cluster))

    N = length(anchors)
    M = length(elements)
    network$x = match(network$anchor, anchors)
    network$y = match(network$cluster, elements)

    mixed.mat = expand.grid(anchors, elements)
    mixed.df = data.frame(anchor=mixed.mat[,1], cluster=mixed.mat[,2])
    m.mixed = make.matrix(mixed.df, x=anchors, y=elements)

    # mixed
    plot.matrix(m=m.mixed, Nx=N, Ny=M, labels.x=anchor.ids, labels.y=element.ids, title="mixed en", col.field="col.en", ofn="mixed_en", plot.network=T)
    plot.matrix(m=m.mixed, Nx=N, Ny=M, labels.x=anchor.ids, labels.y=element.ids, title="mixed obs", col.field="col.obs", ofn="mixed_obs", plot.network=T)
}
