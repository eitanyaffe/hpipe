
plot.hic.map.compare.matrix=function(ifn.map1, ifn.map2, ifn.network1, ifn.network2,
    min.contacts, min.element.abundance, min.anchor.abundance, ifn.elements, ifn.anchors,
    ifn.anchors.compare, ifn.clusters.compare, fdir)
{
    element.table = load.table(ifn.elements)
    anchor.table = load.table(ifn.anchors)

    anchors.compare = load.table(ifn.anchors.compare)
    clusters.compare = load.table(ifn.clusters.compare)

    anchor.table$abundance2 = anchors.compare$anchor.abundance2[match(anchor.table$set, anchors.compare$anchor)]
    element.table$abundance2 = clusters.compare$abundance2[match(element.table$cluster, clusters.compare$cluster)]

    anchor.table = anchor.table[anchor.table$abundance2 > min.anchor.abundance,]
    element.table = element.table[element.table$abundance2 > min.element.abundance,]

    load.map=function(ifn) {
        map = load.table(ifn)
        map$en = ifelse(map$observed>min.contacts, log10(map$observed/map$expected), -Inf)
        map
    }
    map1 = load.map(ifn.map1)
    map2 = load.map(ifn.map2)

    make.triangle.coords=function(m, top) {
        mx = as.matrix(data.frame(x1=m$x-1, x2=m$x-1*top, x3=m$x, dummy=NA))
        my = as.matrix(data.frame(y1=m$y-1, y2=m$y-1*!top, y3=m$y, dummy=NA))
        list(x=as.vector(t(mx)), y=as.vector(t(my)))
    }
    plot.matrix=function(m1, m2, labels.x, labels.y, title, col.field, ofn, label.field, add.labels=F) {
        height = 200+Ny*10
        width = 200+Nx*10
        if (add.labels) {
            height = height * 2
            width = width * 2
        }
        fig.start(fdir=fdir, ofn=paste(fdir, "/", ofn, ".png", sep=""), height=height, width=width)
        plot.new()
        par(mai=c(1,1,1,1))
        plot.window(xlim=c(0,Nx), ylim=c(0,Ny), xaxs="i", yaxs="i")
        title(main=paste(title))

        coords1 = make.triangle.coords(m1, top=T)
        coords2 = make.triangle.coords(m2, top=F)
        polygon(x=coords1$x, y=coords1$y, col=m1[,col.field], border=NA)
        polygon(x=coords2$x, y=coords2$y, col=m2[,col.field], border=NA)
        abline(h=0:Ny, v=0:Nx, col="gray")

        axis(1, at=1:Nx-0.5, labels=labels.x, las=2)
        axis(3, at=1:Nx-0.5, labels=labels.x, las=2)
        axis(2, at=1:Ny-0.5, labels=labels.y, las=2)
        axis(4, at=1:Ny-0.5, labels=labels.y, las=2)
        # points(network1$x-0.75, network1$y-0.25, pch=19, cex=0.5)
        # points(network2$x-0.25, network2$y-0.75, pch=19, cex=0.5)

        if (add.labels) {
            text(m1$x-0.75, m1$y-0.25, m1[,label.field], cex=0.75)
            text(m2$x-0.25, m2$y-0.75, m2[,label.field], cex=0.75)
        }
        fig.end()
    }
    make.matrix=function(df, x, y, map) {
        m = merge(df, map)
        m$x = match(m$anchor, x)
        m$y = match(m$cluster, y)
        m$col.en = panel.en[vals.to.cols(m$en, breaks.en)]
        m$col.en[!is.finite(m$en)] = "lightgray"
        m$col.obs = panel.obs[vals.to.cols(m$observed, breaks.obs)]
        m
    }

    breaks.en = c(0,1,2,3)
    colors.en = c("white", "blue", "red", "orange")
    panel.en = make.color.panel(colors=colors.en)
    wlegend2(fdir=fdir, panel=panel.en, breaks=breaks.en, title="hic_enrichment")

    breaks.obs = c(0,10,100,1000)
    colors.obs = c("white", "blue", "red", "orange")
    panel.obs = make.color.panel(colors=colors.obs)
    wlegend2(fdir=fdir, panel=panel.obs, breaks=breaks.obs, title="hic_obs")

    anchors = anchor.table$set
    anchor.ids = anchor.table$id
    major.cluster = anchors

    network1 = load.table(ifn.network1)
    network2 = load.table(ifn.network2)

    element.ids = element.table$id
    elements = element.table$cluster
    # elements = sort(unique(c(network1$cluster, network2$cluster)))

    Nx = length(anchors)
    Ny = length(elements)
    network1$x = match(network1$anchor, anchors)
    network1$y = match(network1$cluster, elements)
    network2$x = match(network2$anchor, anchors)
    network2$y = match(network2$cluster, elements)

    mat = expand.grid(major.cluster, elements)
    df = data.frame(cluster1=mat[,1], cluster2=mat[,2])

    m1 = make.matrix(df=df, x=major.cluster, y=elements, map=map1)
    m2 = make.matrix(df=df, x=major.cluster, y=elements, map=map2)

    # mixed
    plot.matrix(m1=m1, m2=m2, labels.x=anchor.ids, labels.y=element.ids, title="", col.field="col.en", ofn="enrichment")
    plot.matrix(m1=m1, m2=m2, labels.x=anchor.ids, labels.y=element.ids, title="", col.field="col.en", ofn="enrichment_labels",
                label.field="observed", add.labels=T)
    plot.matrix(m1=m1, m2=m2, labels.x=anchor.ids, labels.y=element.ids, title="", col.field="col.obs", ofn="observed",
                label.field="observed", add.labels=T)
}
