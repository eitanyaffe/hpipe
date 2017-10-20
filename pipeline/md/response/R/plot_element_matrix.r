plot.internal=function(network, anchor.table, element.table, anchor.ids, element.ids, cluster.class, tlegend)
{
    par(xpd=NA)

    N = length(anchor.ids)
    M = length(element.ids)
    K = dim(cluster.class)[2] - 1

    network$anchor.index = match(network$anchor.id, anchor.ids)
    network$element.index = match(network$element.id, element.ids)

    index2id.element=function(index) { network$element.id[match(index, network$element.index)] }
    index2id.anchor=function(index) { network$anchor.id[match(index, network$anchor.index)] }

    ss = split(network$anchor.index, network$element.index)
    smin = sapply(ss, min)
    smax = sapply(ss, max)
    bb.elements = data.frame(index=as.numeric(names(ss)), min=smin, max=smax)

    ss = split(network$element.index, network$anchor.index)
    smin = sapply(ss, min)
    smax = sapply(ss, max)
    bb.anchors = data.frame(index=as.numeric(names(ss)), min=smin, max=smax)

    anchor.ybottom = bb.anchors$min - 0.7

    plot.new()
    par(mai=c(0.5, 0.5, 0.1, 0.1))
    plot.window(xlim=c(0,N), ylim=c(0,M))
    segments(x0=bb.elements$min, x1=bb.elements$max, y0=bb.elements$index, y1=bb.elements$index, col="darkgray", lwd=1)
    segments(y0=bb.anchors$max, y1=anchor.ybottom, x0=bb.anchors$index, x1=bb.anchors$index, col="darkgray", lwd=1)
    text(x=bb.elements$min, y=bb.elements$index, labels=index2id.element(bb.elements$index), cex=0.8, adj=1.5)

    # anchors
    alabels = index2id.anchor(bb.anchors$index)
    sw = strwidth(alabels)
    awidth = max(strheight(alabels))-0.1
    acolors = tlegend$color[match(alabels, tlegend$anchor.id)]
    rect(ybottom=anchor.ybottom-sw-0.2, ytop=anchor.ybottom,
         xleft=bb.anchors$index-awidth/2, xright=bb.anchors$index+awidth/2, border=NA, col=acolors)
    text(y=anchor.ybottom, x=bb.anchors$index, srt=90, labels=index2id.anchor(bb.anchors$index), cex=0.8, adj=1.2)

    # connecting nodes
    points(x=network$anchor.index, y=network$element.index, pch=19, col="white")
    points(x=network$anchor.index, y=network$element.index, pch=1, col=1)

    class.size = 0.5
    class.offset.1 = 1.3
    class.offset.2 = 2.1
    class.width = 2
    cluster.class = cluster.class[is.element(cluster.class$cluster, element.table$cluster),]
    df = cluster.class[match(element.table$cluster, cluster.class$cluster),-1]
    types = names(df)
    K = length(types)
    right = ifelse(bb.elements$index<10, bb.elements$min - class.offset.1, bb.elements$min - class.offset.2)
    left = right - class.width
    for (i in 1:length(types)) {
        type = types[i]
        type.on = df[,i] == 1
        xleft = left + class.width*(i-1)/K
        xright = left + class.width*(i)/K
        rect(xleft=xleft, xright=xright, ybottom=bb.elements$index-class.size/2, ytop=bb.elements$index+class.size/2, col="gray", border="darkgray")
        rect(xleft=xleft, xright=xright, ybottom=bb.elements$index-class.size/2, ytop=bb.elements$index+class.size/2, col=ifelse(type.on, "blue", NA),
             border="darkgray")
    }
}

plot.element.matrix=function(ifn.anchors, ifn.elements, ifn.network, ifn.cluster.class, ifn.taxa.legend, fdir)
{
    anchor.table = load.table(ifn.anchors)
    element.table = load.table(ifn.elements)
    network = load.table(ifn.network)
    cluster.class = load.table(ifn.cluster.class)
    tlegend = load.table(ifn.taxa.legend)
    types = colnames(cluster.class)[-1]

    anchor.ids = anchor.table$id
    element.ids = element.table$id

    network$anchor.id = anchor.table$id[match(network$anchor, anchor.table$set)]
    network$element.id = element.table$id[match(network$cluster, element.table$cluster)]

    N = length(anchor.ids)
    M = length(element.ids)
    K = dim(cluster.class)[2] - 1

    fig.start(fdir=fdir, ofn=paste(fdir, "/element_matrix.pdf", sep=""), type="pdf", width=1+0.13*N, height=1+0.13*M)
    plot.internal(network=network, anchor.ids=anchor.ids, element.ids=element.ids,
                  anchor.table=anchor.table, element.table=element.table, cluster.class=cluster.class, tlegend=tlegend)
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/element_matrix_legend.pdf", sep=""), type="pdf", width=0.5+0.5*K, height=2.2)
    plot.new()
    par(mai=c(1.5, 0, 0, 0))
    plot.window(xlim=c(0,K), ylim=c(0,1))
    for (i in 1:length(types))
        rect(xleft=i-1, xright=i, ybottom=-1, ytop=1, col="blue", border="darkgray")
    axis(1, at=1:K-0.5, labels=types, las=2, tick=F, line=-0.5)
    fig.end()
}

plot.element.matrix.compare=function(
    ifn.network, ifn.elements, ifn.anchors, ifn.clusters, ifn.map, ifn.map.select,
    ifn.cluster.class, ifn.taxa.legend, fdir)
{
    network = load.table(ifn.network)
    element.table = load.table(ifn.elements)
    map = load.table(ifn.map)
    map.select = load.table(ifn.map.select)
    anchor.table = load.table(ifn.anchors)
    cluster.table = load.table(ifn.clusters)
    cluster.class = load.table(ifn.cluster.class)
    tlegend = load.table(ifn.taxa.legend)
    types = colnames(cluster.class)[-1]

    anchors = sort(unique(map.select$anchor))
    clusters = sort(unique(map.select$cluster))

    network = network[is.element(network$anchor, anchors) | is.element(network$cluster, clusters),]
    network$anchor.id = anchor.table$anchor.id[match(network$anchor, anchor.table$anchor)]
    network$element.id = element.table$id[match(network$cluster, element.table$cluster)]

    anchor.ids = anchor.table$anchor.id[is.element(anchor.table$anchor, network$anchor)]
    element.ids = element.table$id[is.element(element.table$cluster, network$cluster)]

    N = length(anchors)
    M = length(clusters)
    K = dim(cluster.class)[2] - 1

    fig.start(fdir=fdir, ofn=paste(fdir, "/element_matrix.pdf", sep=""), type="pdf", width=2+0.13*N, height=2+0.13*M)
    plot.internal(network=network, anchor.ids=anchor.ids, element.ids=element.ids,
                  anchor.table=anchor.table, element.table=element.table, cluster.class=cluster.class, tlegend=tlegend)
    fig.end()
}
