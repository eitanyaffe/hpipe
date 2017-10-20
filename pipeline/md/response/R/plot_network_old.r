plot.network.elements=function(lid, ifn.network, ifn.obs, ifn.exp, labels, reg, fdir)
{
    obs = load.table(ifn.obs)
    exp = load.table(ifn.exp)
    mat = load.table(ifn.network)
    elements = sort(unique(mat$cluster))
    N = dim(obs)[2] - 1
    coords = 1:N

    ylim = c(-3.5, 1.5)

    # separate files
    for (element in elements) {
        xx = mat[mat$cluster == element,]
        M = dim(xx)[2]
        majors = xx$major.cluster
        anchors = xx$anchor

        major.obs = obs[match(majors, obs$cluster), -1]
        major.exp = exp[match(majors, exp$cluster), -1]
        major.norm = log10((major.obs+reg)/(major.exp+reg))
        major.norm.sum = log10(colSums(10^major.norm))

        element.obs = unlist(obs[element == obs$cluster, -1])
        element.exp = unlist(exp[element == exp$cluster, -1])
        element.norm = log10((element.obs+reg)/(element.exp+reg))

        # confidence interval
        element.sqrt = sqrt(element.obs)
        element.high = element.obs+element.sqrt
        element.low = element.obs-element.sqrt
        element.low = ifelse(element.low < 0, 0, element.low)
        high = log10((element.high+reg)/(element.exp+reg))
        low = log10((element.low+reg)/(element.exp+reg))

        fig.start(fdir=fdir, ofn=paste(fdir, "/", element, ".png", sep=""), height=400, width=300)
        plot.init(xlim=c(1,N), ylim=ylim,
                  main=paste(lid, ", eid=", element, ", anchors=", paste(anchors, sep=",", collapse=","), sep=""), x.axis=F)
        abline(h=0)
        for (i in 1:M)
            lines(x=coords, y=major.norm[i,], col=i+1)
        if (M>1)
            lines(x=coords, y=major.norm.sum, lty=2)

        polygon(x=c(coords,rev(coords)), y=c(high, rev(low)), col="gray", border=NA)
        lines(x=coords, y=element.norm, lwd=2)
        at = c(3,7)+0.5
        abline(v=at, lwd=2, lty=2)
        axis(side=1, labels=labels, at=1:N, las=2)
        fig.end()

        wlegend(fdir, names=anchors, cols=1:M+1, title="anchor_colors", ofn.prefix=element, width=600, make.file=T, border=NA)

    }

}

plot.network=function(ifn.majors, ifn.network, ifn.coords, ifn.taxa.legend, ifn.element.info, ifn.classify, lid, fdir)
{
    majors = load.table(ifn.majors)
    mat = load.table(ifn.network)
    taxa.legend = load.table(ifn.taxa.legend)
    info = load.table(ifn.element.info)
    coords = load.table(ifn.coords)
    class.table = load.table(ifn.classify)

    a.ids = unique(mat$anchor)
    e.ids = unique(mat$cluster)

    n.anchors = length(a.ids)
    n.elements = length(e.ids)

    df.anchors = data.frame(id=paste("g", a.ids, sep=""), type="anchor", class="anchor",
        color=taxa.legend$group.color[match(a.ids, taxa.legend$anchor)], shape="square", label=a.ids, size=6)

    df.elements = data.frame(id=e.ids, type="element", class=info$class[match(e.ids, info$cluster)], shape="circle", label=e.ids, size=2)

    df.elements$color =
        ifelse(df.elements$class == "none", "lightgray",
               ifelse(df.elements$class == "phage", "red",
                      ifelse(df.elements$class == "plamid", "blue",
                             ifelse(df.elements$class == "transposase", "green", "orange"))))

    vs = rbind(df.anchors, df.elements)
    vs$x = coords$x[match(vs$id, coords$id)]
    vs$y = coords$y[match(vs$id, coords$id)]

    es = data.frame(from=paste("g", mat$anchor, sep=""), to=mat$cluster)
    es = es[is.element(es$from, vs$id) & is.element(es$to, vs$id),]

    es$from.x = vs$x[match(es$from, vs$id)]
    es$from.y = vs$y[match(es$from, vs$id)]
    es$to.x = vs$x[match(es$to, vs$id)]
    es$to.y = vs$y[match(es$to, vs$id)]

    vs.anchors = vs[vs$type == "anchor",]
    vs.elements = vs[vs$type == "element",]

    # determine position for element labels
    rad2deg = function(rad) {(rad * 180) / (pi)}
    vs.elements$pos = 0
    for (i in 1:dim(vs.elements)[1])
    {
        id = vs.elements$id[i]
        oids = es$from[es$to == id]
        mm = data.frame(
            x=vs$x[match(oids, vs$id)] - vs$x[match(id, vs$id)],
            y=vs$y[match(oids, vs$id)] - vs$y[match(id, vs$id)])
        mm$length = sqrt(mm$x^2 + mm$y^2)
        mm$x = mm$x / mm$length
        mm$y = mm$y / mm$length

        x = -sum(mm$x)
        y = -sum(mm$y)
        tetha = rad2deg(atan2(y=y, x=x))
        if (tetha < -45 && tetha > -135)
            pos = 1
        if (tetha < -135 || tetha > 135)
            pos = 2
        if (tetha < 135 && tetha > 45)
            pos = 3
        if (tetha < 45 && tetha > -45)
            pos = 4
        vs.elements$pos[i] = pos
    }

    pplot=function(ofn, labels) {
        fig.start(fdir=fdir, ofn=ofn, height=800, width=800)
        plot.new()
        plot.window(xlim=expand.range(vs$x, 0.1), ylim=expand.range(vs$y, 0.1))

        # plot edges
        segments(x0=es$from.x, x1=es$to.x, y0=es$from.y, y1=es$to.y)

        # anchor vertices
        points(x=vs.anchors$x, y=vs.anchors$y, pch=19, col=vs.anchors$color, cex=4)
        # points(x=vs.anchors$x, y=vs.anchors$y, cex=4)

        # element vertices
        points(x=vs.elements$x, y=vs.elements$y, pch=15, col=vs.elements$color, cex=2)

        if (labels) {
            text(x=vs.anchors$x, y=vs.anchors$y, labels=vs.anchors$label)
            text(x=vs.elements$x, y=vs.elements$y, labels=vs.elements$label, pos=vs.elements$pos, offset=0.7)
        }

        title(main=paste(lid, ": n.anchors=", dim(vs.anchors)[1], ", n.elements=", dim(vs.elements)[1], sep=""))
        fig.end()
    }

    pplot(ofn=paste(fdir, "/element_network.png", sep=""), labels=T)
    pplot(ofn=paste(fdir, "/element_network_clean.png", sep=""), labels=F)

}
