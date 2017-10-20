plot.network.compare.details=function(ifn.elements, ifn.anchors, ifn.clusters, ifn.map, min.contacts, ifn.map.select, legend1, legend2, fdir)
{
    system(paste("rm -rf", fdir))
    element.table = load.table(ifn.elements)
    map = load.table(ifn.map)
    anchor.table = load.table(ifn.anchors)
    cluster.table = load.table(ifn.clusters)
    map.select = load.table(ifn.map.select)

    anchor.table = anchor.table[!anchor.table$lost,]
    cluster.table = cluster.table[!cluster.table$lost,]

    # limit to relevant anchors
    map = map[is.element(map$anchor, anchor.table$anchor),]

    element.table = element.table[is.element(element.table$cluster, cluster.table$cluster),]
    map$is.element = is.element(map$cluster, element.table$cluster)

    # map = map[map$is.element,]
    map$element.id = ifelse(map$is.element, element.table$id[match(map$cluster, element.table$cluster)], NA)

    map$anchor.id = anchor.table$anchor.id[match(map$anchor, anchor.table$anchor)]

    ix = match(map$mid, map.select$mid)
    map$type = ifelse(!is.na(ix), map.select$type[ix], "unknown")
    map$select = map$type != "unknown"
    map$col =
        ifelse(map$type == "gained", "red",
               ifelse(map$type == "lost", "blue",
                      ifelse(map$type == "stable", "green", "black")))

    map = map[map$observed1 >= min.contacts | map$observed2 >= min.contacts,]

    # anchor.table = anchor.table[is.element(anchor.table$anchor, sort(unique(map$anchor))),]
    # cluster.table = cluster.table[is.element(cluster.table$cluster, sort(unique(map$cluster))),]
    map = map[is.element(map$cluster, cluster.table$cluster) & is.element(map$anchor, anchor.table$anchor),]
    element.table = element.table[is.element(element.table$cluster, map$cluster),]

    lim = range(c(map$score1, map$score2))

    plot.map=function(imap, add.labels=F, main, multi=F) {
        ilim = ifelse(rep(multi,2), lim, 1.1*range(c(0, imap$score1, imap$score2)))
        if (ilim[2] < 1) ilim[2] = 1
        plot.init(xlim=ilim, ylim=ilim, xlab=legend1, ylab=legend2, main=main, x.axis=!multi, y.axis=!multi)
        if (multi) {
            axis(1, labels=F)
            axis(2, labels=F)
        }
        abline(a=0, b=1, col=1, lty=3)
        abline(h=0, col=1, lty=3)
        abline(v=0, col=1, lty=3)

        left = ifelse(imap$score1 != 0, imap$score1 - imap$sd.score1, 0)
        right = imap$score1 + imap$sd.score1
        bottom = ifelse(imap$score2 != 0, imap$score2 - imap$sd.score2, 0)
        top = imap$score2 + imap$sd.score2
        segments(x0=imap$score1, x1=imap$score1, y0=bottom, y1=top, col="gray", lwd=2)
        segments(x0=left, x1=right, y0=imap$score2, y1=imap$score2, col="gray", lwd=2)

        ix = imap$select
        points(x=imap$score1[!ix], y=imap$score2[!ix], col="black", pch=19, cex=0.5)
        points(x=imap$score1[ix], y=imap$score2[ix], col=imap$col[ix], pch=19)
        # points(x=imap$score1, y=imap$score2, col=ifelse(imap$select,1,NA), pch=1)

        if (add.labels && dim(imap)[1] > 0 && is.element("label", colnames(imap)))
            text(x=imap$score1, y=imap$score2, pos=4, labels=imap$label, cex=0.75)
    }

    height = 6
    width = 5

    height.clean = height * 0.75
    width.clean = width * 0.75

    fdir.anchors.clean = paste(fdir, "/anchors/clean", sep="")
    fdir.anchors.labels = paste(fdir, "/anchors/labels", sep="")
    fdir.elements.clean = paste(fdir, "/elements/clean", sep="")
    fdir.elements.labels = paste(fdir, "/elements/labels", sep="")
    system(sprintf("mkdir -p %s %s %s %s", fdir.anchors.clean, fdir.anchors.labels, fdir.elements.clean, fdir.elements.labels))

    for (i in 1:dim(anchor.table)[1]) {
        anchor = anchor.table$anchor[i]
        anchor.id = anchor.table$anchor.id[i]
        a1 = round(anchor.table$abundance1[i],2)
        a2 = round(anchor.table$abundance2[i],2)
        main = paste(anchor.id, ", a1=", a1, ", a2=", a2, sep="")
        imap = map[map$anchor == anchor,]
        imap$label = imap$element.id

        fig.start(fdir=fdir.anchors.clean, ofn=paste(fdir.anchors.clean, "/", anchor.id, ".pdf", sep=""),
                  type="pdf", height=height.clean, width=width.clean)
        plot.map(imap=imap, add.labels=F, main=main)
        fig.end()

        fig.start(fdir=fdir.anchors.labels, ofn=paste(fdir.anchors.labels, "/", anchor.id, ".pdf", sep=""),
                  type="pdf", height=height, width=width)
        plot.map(imap=imap, add.labels=T, main=main)
        fig.end()
    }

    for (i in 1:dim(element.table)[1]) {
        cluster = element.table$cluster[i]
        element.id = element.table$id[i]
        a1 = round(cluster.table$abundance1[match(cluster, cluster.table$cluster)],2)
        a2 = round(cluster.table$abundance2[match(cluster, cluster.table$cluster)],2)
        main = paste(element.id, ", a1=", a1, ", a2=", a2, sep="")
        imap = map[map$cluster == cluster,]
        imap$label = imap$anchor.id

        fig.start(fdir=fdir.elements.clean, ofn=paste(fdir.elements.clean, "/", element.id, ".pdf", sep=""),
                  type="pdf", height=height.clean, width=width.clean)
        plot.map(imap=imap, add.labels=F, main=main)
        fig.end()

        fig.start(fdir=fdir.elements.labels, ofn=paste(fdir.elements.labels, "/", element.id, ".pdf", sep=""),
                  type="pdf", height=height, width=width)
        plot.map(imap=imap, add.labels=T, main=main)
        fig.end()
    }

    # one anchor large plot
    N = dim(anchor.table)[1]
    Ny = ceiling(sqrt(N))
    Nx = ceiling(N/Ny)
    NN = c(1:N,rep(N+1,Nx*Ny-N))
    fig.start(fdir=fdir, ofn=paste(fdir, "/anchor_summary.pdf", sep=""), type="pdf", height=Ny*1, width=Nx*1)
    layout(matrix(NN, Nx, Ny, byrow=T))
    par(mai=c(0.1, 0.1, 0.25, 0.05))
    for (i in 1:dim(anchor.table)[1]) {
        anchor = anchor.table$anchor[i]
        anchor.id = anchor.table$anchor.id[i]
        imap = map[map$anchor == anchor,]
        plot.map(imap=imap, add.labels=F, main=anchor.id, multi=T)
    }
    fig.end()


    # one element large plot
    N = dim(element.table)[1]
    Ny = ceiling(sqrt(N))
    Nx = ceiling(N/Ny)
    NN = c(1:N,rep(N+1,Nx*Ny-N))
    fig.start(fdir=fdir, ofn=paste(fdir, "/element_summary.pdf", sep=""), type="pdf", height=Ny*1, width=Nx*1)
    layout(matrix(NN, Nx, Ny, byrow=T))
    par(mai=c(0.1, 0.1, 0.25, 0.05))
    for (i in 1:dim(element.table)[1]) {
        cluster = element.table$cluster[i]
        element.id = element.table$id[i]
        imap = map[map$cluster == cluster,]
        plot.map(imap=imap, add.labels=F, main=element.id, multi=T)
    }
    fig.end()
}
