plot.network.compare.matrix=function(ifn.elements, ifn.anchors, ifn.clusters, ifn.map, min.contacts, ifn.map.select, zoom, fdir)
{
    make.triangle.coords=function(m, top) {
        mx = as.matrix(data.frame(x1=m$xcoord-1, x2=m$xcoord-1*top, x3=m$xcoord, dummy=NA))
        my = as.matrix(data.frame(y1=m$ycoord-1, y2=m$ycoord-1*!top, y3=m$ycoord, dummy=NA))
        list(x=as.vector(t(mx)), y=as.vector(t(my)))
    }

    element.table = load.table(ifn.elements)
    map = load.table(ifn.map)
    map.select = load.table(ifn.map.select)
    anchor.table = load.table(ifn.anchors)
    cluster.table = load.table(ifn.clusters)

    map$omit1 = map$observed1 < min.contacts
    map$omit2 = map$observed2 < min.contacts
    map$omit = map$omit1 & map$omit2

    # focus on anchors that have hi-c results (low abundance anchors are missing)
    valid.anchors = sort(unique(map$anchor))
    anchor.table = anchor.table[is.element(anchor.table$anchor, valid.anchors) & !anchor.table$lost,]

    cluster.table = cluster.table[!cluster.table$lost,]
    element.table = element.table[is.element(element.table$cluster, cluster.table$cluster),]

    map.select = map.select[map.select$type != "stable" & map.select$type != "unknown",]

    if (zoom) {
        sclusters = map.select$cluster
        xmap = map[is.element(map$cluster, sclusters),]
        sanchors = unique(c(map.select$anchor, xmap$anchor[xmap$type == "stable"]))
        anchor.table = anchor.table[is.element(anchor.table$anchor, sanchors),]
        element.table = element.table[is.element(element.table$cluster, sclusters),]
    }

    anchor.ids = anchor.table$anchor.id
    element.ids = element.table$id

    map = map[is.element(map$cluster, element.table$cluster) & is.element(map$anchor, anchor.table$anchor),]
    map$element.id = element.table$id[match(map$cluster, element.table$cluster)]
    map$anchor.id = anchor.table$anchor.id[match(map$anchor, anchor.table$anchor)]

    # add coords
    map$xcoord = match(map$anchor.id, anchor.ids)
    map$ycoord = match(map$element.id, element.ids)

    Nx = length(anchor.ids)
    Ny = length(element.ids)
    labels.x = anchor.ids
    labels.y = element.ids

    height = 2+Ny*0.15
    width = 1+Nx*0.15

    if (zoom) {
        height = height + 1.5
        width = width + 1.5
    }

    colors = c("blue", "white", "red")
    panel = make.color.panel(colors=colors)
    breaks = c(-3, 0, 3)
    wlegend2(fdir=fdir, panel=panel, breaks=breaks, title="delta_score")

    breaks.en = c(0,4)
    colors.en = c("white","red")
    panel.en = make.color.panel(colors=colors.en)
    wlegend2(fdir=fdir, panel=panel.en, breaks=breaks.en, title="score")

    map$col = panel[vals.to.cols(map$delta.score, breaks)]
    map$col[map$omit] = "lightgray"

    map$col1 = panel.en[vals.to.cols(map$score1, breaks.en)]
    map$col2 = panel.en[vals.to.cols(map$score2, breaks.en)]
    map$col1[map$omit1] = "lightgray"
    map$col2[map$omit2] = "lightgray"
    map$col =
        ifelse(map$type == "gained", "red",
               ifelse(map$type == "lost", "blue",
                      ifelse(map$type == "stable", "green", "white")))

    names = c("gained", "increase", "stable", "decrease", "lost")
    cols = c("red", "pink", "green", "lightblue", "blue")
    wlegend(fdir=fdir, names=names, cols, title="delta", border=1)

    map$select.sign =
        ifelse(map$type == "gained", "+",
               ifelse(map$type == "lost", "-",
                      ifelse(map$type == "stable", "=", "")))

    plot.f=function(trig) {
        trig.tag = ifelse(trig, "cmp", "delta")
        zoom.tag = ifelse(zoom, "_zoom", "")
        ofn = sprintf("%s/%s%s.pdf", fdir, trig.tag, zoom.tag)
        fig.start(fdir=fdir, ofn=ofn, type="pdf", height=height, width=width)
        plot.new()
        par(mai=c(1,1,1,1))
        plot.window(xlim=c(0,Nx), ylim=c(0,Ny), xaxs="i", yaxs="i")
        if (trig) {
            top = make.triangle.coords(map, top=T)
            bottom = make.triangle.coords(map, top=F)
            polygon(x=top$x, y=top$y, col=map$col1, border=NA)
            polygon(x=bottom$x, y=bottom$y, col=map$col2, border=NA)
            if (zoom) text(x=map$xcoord-0.5, y=map$ycoord-0.5, labels=map$select.sign)
        } else {
            rect(xleft=map$xcoord-1, xright=map$xcoord, ybottom=map$ycoord-1, ytop=map$ycoord, border=NA, col=map$col)
            # points(x=map$xcoord-0.5, y=map$ycoord-0.5, col=ifelse(map$select, 1, NA), pch=19, cex=0.75)
        }
        abline(v=0:Nx, col="darkgray")
        abline(h=0:Ny, col="darkgray")
        axis(1, at=1:Nx-0.5, labels=labels.x, las=2)
        axis(2, at=1:Ny-0.5, labels=labels.y, las=2)
        if (!zoom) {
            axis(3, at=1:Nx-0.5, labels=labels.x, las=2)
            axis(4, at=1:Ny-0.5, labels=labels.y, las=2)
        }
        fig.end()
    }
    plot.f(T)
    plot.f(F)
}
