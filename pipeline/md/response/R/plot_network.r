plot.network.element.support=function(ifn.anchors, ifn.clusters, ifn.map, ifn.elements, ifn.network, ifn.ca, ifn.ca.matrix, fdir)
{
    anchor.table = load.table(ifn.anchors)
    element.table = load.table(ifn.elements)
    emap = load.table(ifn.map)
    cluster.table = load.table(ifn.clusters)
    ca = load.table(ifn.ca)
    ca.matrix = load.table(ifn.ca.matrix)
    mat = load.table(ifn.network)
    clusters = sort(unique(mat$cluster))

    ca.matrix$o = ca.matrix$contig_total_count
    ca.matrix$e = ca.matrix$contig_expected

    emap$o = emap$sum.observed
    emap$e = emap$sum.expected

    xlim = range(ca.matrix$e)
    ylim = c(1,max(ca.matrix$o))

    for (cluster in clusters) {
        emat = mat[mat$cluster == cluster,]
        anchors = emat$anchor
        contigs = cluster.table$contig[cluster.table$cluster == cluster]
        element.id = element.table$id[match(cluster, element.table$cluster)]
        for (anchor in anchors) {
            anchor.id = anchor.table$id[match(anchor, anchor.table$set)]
            emap.x = emap[emap$cluster == cluster & emap$anchor == anchor,]
            anchor.contigs = ca$contig[ca$anchor == anchor & ca$contig_anchor != 0]
            cmat = ca.matrix[ca.matrix$anchor == anchor,]
            cmat.a = cmat[is.element(cmat$contig, anchor.contigs),]

            cmat.i = ca.matrix[ca.matrix$anchor == anchor & ca.matrix$contig_anchor != 0 & ca.matrix$contig_anchor != anchor,]

            # limit to element
            cmat = cmat[is.element(cmat$contig,contigs),]
            zeros = sum(!is.element(contigs,cmat$contig))

            # are contigs associated with anchor?
            selected.contigs = ca$contig[ca$anchor == anchor & is.element(ca$contig, contigs)]
            cmat$color = ifelse(is.element(cmat$contig, selected.contigs), "red", "blue")
            zeros = zeros + sum(cmat$contig_total_count == 0)
            main = paste("element=", element.id, ", anchor=", anchor.id, "\nassociated=", sum(is.element(cmat$contig, selected.contigs)), "/", length(contigs), "\nno obs=", zeros, sep="")
            fig.start(fdir=fdir, ofn=paste(fdir, "/", element.id, "_", anchor.id, ".pdf", sep=""), type="pdf", height=4, width=3)

            plot.init(xlim=xlim, ylim=ylim, log="xy", main=main)
            abline(a=0,b=1, lty=2)
            points(cmat.i$e, cmat.i$o, col="gray", pch=19, cex=0.5)
            points(cmat.a$e, cmat.a$o, col="black", pch=19, cex=0.5)
            points(cmat$e, cmat$o, col=cmat$color, pch=19)
            points(emap.x$e, emap.x$o, col="green", pch=19)

            fig.end()
        }
    }
}

plot.network.element.summary=function(lid, ifn.reps, ifn.anchors, ifn.elements, ifn.network, ifn.obs, ifn.exp, ifn.detection, labels, fdir)
{
    rep.table = load.table(ifn.reps)
    element.table = load.table(ifn.elements)
    anchor.table = load.table(ifn.anchors)
    obs = load.table(ifn.obs)
    exp = load.table(ifn.exp)
    mat = load.table(ifn.network)
    min.score = log10(load.table(ifn.detection)[1,1])
    elements = sort(unique(mat$cluster))

    result = NULL
    for (element in elements) {
        element.id = element.table$id[match(element, element.table$cluster)]
        xx = mat[mat$cluster == element,]
        M = dim(xx)[1]
        anchors = xx$anchor
        majors = anchors

        major.obs = obs[match(majors, obs$cluster), -1]
        major.exp = exp[match(majors, exp$cluster), -1]
        major.norm = log10(major.obs/major.exp)
        major.norm[major.norm<min.score] = min.score
        major.norm.sum = log10(colSums(ifelse(major.norm>min.score,10^major.norm,0)))
        major.norm.sum[major.norm.sum<min.score] = min.score

        element.obs = unlist(obs[element == obs$cluster, -1])
        element.exp = unlist(exp[element == exp$cluster, -1])
        element.norm = log10(element.obs/element.exp)
        element.norm[element.norm<min.score] = min.score

        delta = element.norm - major.norm.sum
        delta = delta - delta[1]
        delta = log2(10^delta)
        result = rbind(result, delta)
    }

    mm = sapply(as.data.frame(result), mean)
    sd = sapply(as.data.frame(result), sd)
    top = mm + sd
    bottom = mm - sd
    middle = mm

    #qnt = sapply(as.data.frame(result), function(x) { quantile(x, c(0.25, 0.5, 0.75))})
    #bottom = qnt[1,]
    #middle = qnt[2,]
    #top = qnt[3,]

    N = length(labels)
    ylim = range(c(top, bottom))
    coord = 1:N-0.5
    width = 0.5
    xleft = coord-width/2
    xright = coord+width/2
    fig.start(fdir=fdir, ofn=paste(fdir, "/element_trend.pdf", sep=""), type="pdf", height=4, width=4)
    plot.init(xlim=c(0,N), ylim=ylim, x.axis=F, y.axis=F)
    grid()
    #rect(xleft=xleft, xright=xright, ybottom=bottom, ytop=top, col="gray", border=NA)
    #segments(x0=xleft, x1=xright, y0=middle, y1=middle, col=1, lwd=2)

    segments(x0=xleft, x1=xright, y0=top, y1=top, col=1, lwd=1)
    segments(x0=xleft, x1=xright, y0=bottom, y1=bottom, col=1, lwd=1)
    segments(x0=coord, x1=coord, y0=bottom, y1=top, col=1, lwd=1)
    points(x=coord, y=middle, pch=15, col=1, cex=1)
    axis(2, las=2)
    axis(side=1, labels=labels, at=1:N, las=2)
    fig.end()
}

plot.network=function(ifn.reps, ifn.majors, ifn.network, ifn.coords, ifn.taxa.legend, ifn.elements, ifn.classify, lid, fdir)
{
    rfactor = 60
    size = 12

    majors = load.table(ifn.majors)
    network = load.table(ifn.network)

    taxa.legend = load.table(ifn.taxa.legend)
    element.table = load.table(ifn.elements)
    coords = load.table(ifn.coords)

    # classes
    class.table = load.table(ifn.classify)
    n.classes = dim(class.table)[2] - 1
    classes = names(class.table)[-1]
    class.colors = rainbow(n.classes)

    network = network[,c("anchor", "cluster")]

    anchors = unique(network$anchor)
    elements = unique(network$cluster)

    n.anchors = length(anchors)
    n.elements = length(elements)

    ix = match(anchors, taxa.legend$anchor)
    anchor.table = data.frame(
        id=anchors,
        anchor.id=taxa.legend$anchor.id[match(anchors, taxa.legend$anchor)],
        color=taxa.legend$color[ix],
        sign=paste(taxa.legend$letter[ix], sep=""))
    anchor.table$label = anchor.table$anchor.id

    anchor.table$x = coords$x[match(paste("g", anchor.table$id, sep=""), coords$id)]
    anchor.table$y = coords$y[match(paste("g", anchor.table$id, sep=""), coords$id)]

    etable = data.frame(id=elements)
    etable$x = coords$x[match(etable$id, coords$id)]
    etable$y = coords$y[match(etable$id, coords$id)]

    elements.classes = class.table[match(elements, class.table$cluster),]

    network$anchor.x = anchor.table$x[match(network$anchor, anchor.table$id)]
    network$anchor.y = anchor.table$y[match(network$anchor, anchor.table$id)]
    network$element.x = etable$x[match(network$cluster, etable$id)]
    network$element.y = etable$y[match(network$cluster, etable$id)]

    # determine position for element labels
    rad2deg = function(rad) {(rad * 180) / (pi)}
    etable$pos = 0
    for (i in 1:dim(etable)[1])
    {
        id = etable$id[i]
        oids = network$anchor[network$cluster == id]
        mm = data.frame(
            x=anchor.table$x[match(oids, anchor.table$id)] - etable$x[match(id, etable$id)],
            y=anchor.table$y[match(oids, anchor.table$id)] - etable$y[match(id, etable$id)])
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
        etable$pos[i] = pos
    }

    xlim = expand.range(c(anchor.table$x, etable$x), 0.1)
    ylim = expand.range(c(anchor.table$y, etable$y), 0.1)
    radius = min(diff(xlim), diff(ylim))/rfactor

    circle.coords=function(df, N, offset=0, radius=1) {
        for (i in 1:N) {
            rot = 2*((i-1+0.25+offset)/N)
            df[,paste("x",i,sep="_")] = df$x + cospi(rot)*radius
            df[,paste("y",i,sep="_")] = df$y + sinpi(rot)*radius
        }
        df
    }
    triangle.coords=function(df, N, prefix="x") {
        df$null = NA
        fields = NULL
        for (i in 1:(N-1))
            fields = c(fields, c(prefix, paste(prefix, i, sep="_"), paste(prefix, i+1, sep="_"), "null"))
        fields = c(fields, c(prefix, paste(prefix, 1, sep="_"), paste(prefix, N, sep="_"), "null"))
        m = as.matrix(df[,fields])
        as.vector(t(m))
    }
    triangle.colors=function(df, colors) {
        m = as.matrix(df)
        v = as.vector(t(m))
        xcolors = rep(colors, length(v))
        ifelse(v, xcolors, "white")
    }
    get.triangles=function(df, N, radius, colors, class.matrix) {
        df = circle.coords(df=df, N=N, radius=radius)
        tri.coords.x = triangle.coords(df=df, N=N, prefix="x")
        tri.coords.y = triangle.coords(df=df, N=N, prefix="y")
        tri.colors = triangle.colors(df=class.matrix, class.colors)
        list(x=tri.coords.x, y=tri.coords.y, color=tri.colors)
    }
    element.triangles = get.triangles(df=etable, N=n.classes, radius=radius, colors=class.colors, class.matrix=elements.classes[,-1] == 1)

    # prepare legend
    legend.triangles = get.triangles(df=data.frame(x=0,y=0), N=n.classes, radius=0.5, colors=class.colors, class.matrix=matrix(T,1,n.classes))
    legend.text.points = circle.coords(df=data.frame(x=0,y=0), N=n.classes, radius=0.75, offset=0.5)

    fig.type = "pdf"

    plegend=function(clean) {
        fig.start(fdir=fdir, type=fig.type, ofn=paste(fdir, "/element_legend", ifelse(clean,"_clean", ""), ".", fig.type, sep=""), height=4, width=4)
        plot.new()
        par(mai=c(0,0,0,0))
        plot.window(xlim=c(-1,1), ylim=c(-1,1))
        polygon(x=legend.triangles$x, y=legend.triangles$y, col=legend.triangles$color)
        if (!clean) {
            for (i in 1:n.classes)
                text(x=legend.text.points[1,paste("x", i, sep="_")], y=legend.text.points[1,paste("y", i, sep="_")], labels=classes[i])
        }
        fig.end()
    }
    plegend(T)
    plegend(F)

    aspect.ratio = diff(xlim) / diff(ylim)
    pplot=function(ofn, labels) {
        fig.start(fdir=fdir, type=fig.type, ofn=ofn, height=size, width=size*aspect.ratio)
        plot.new()
        plot.window(xlim=xlim, ylim=ylim)

        # plot edges
        segments(x0=network$anchor.x, x1=network$element.x, y0=network$anchor.y, y1=network$element.y, lwd=1, col="darkgray")

        # anchor vertices
        points(x=anchor.table$x, y=anchor.table$y, pch=19, col=anchor.table$color, cex=3)
        # points(x=vs.anchors$x, y=vs.anchors$y, cex=4)

        # polygon(x=element.triangles$x, y=element.triangles$y, col=element.triangles$color)

        text(x=anchor.table$x, y=anchor.table$y, labels=anchor.table$label, cex=0.75)
        if (labels) {
            text(x=etable$x, y=etable$y, labels=element.table$id[match(etable$id, element.table$cluster)], pos=etable$pos, offset=0.5, cex=0.75)
        }

        # element vertices
        sign =
            ifelse(elements.classes$phage | elements.classes$phage, "p",
                   ifelse(elements.classes$defense, "d",
                          ifelse(elements.classes$recombination, "r", "")))
        points(x=etable$x, y=etable$y, pch=19, col="white", cex=2)
        points(x=etable$x, y=etable$y, pch=1, cex=2)
        text(x=etable$x, y=etable$y, sign, cex=0.8)


        title(main=paste(lid, ": n.anchors=", dim(anchor.table)[1], ", n.elements=", dim(etable)[1], sep=""))
        fig.end()
    }

    pplot(ofn=paste(fdir, "/element_network.", fig.type, sep=""), labels=T)
    pplot(ofn=paste(fdir, "/element_network_clean.", fig.type, sep=""), labels=F)

}
