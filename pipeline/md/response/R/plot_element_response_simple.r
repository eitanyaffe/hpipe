plot.element.response.simple=function(
    ifn.reps, ifn.anchors, ifn.elements, ifn.contigs, ifn.norm, ifn.bottom5, ifn.bottom25, ifn.median, ifn.top75, ifn.top95,
    ifn.detection, labels, disturb.ids, fdir)
{
    rep.table = load.table(ifn.reps)
    element.table = load.table(ifn.elements)
    anchor.table = load.table(ifn.anchors)
    q5 = load.table(ifn.bottom5)
    q25 = load.table(ifn.bottom25)
    q50 = load.table(ifn.median)
    q75 = load.table(ifn.top75)
    q95 = load.table(ifn.top95)
    min.score = log10(load.table(ifn.detection)[1,1])
    norm = load.table(ifn.norm)
    contigs = load.table(ifn.contigs)

    elements = element.table$id
    N = length(elements)

    ids = colnames(q50)[-(1:2)]
    M = length(ids)
    drange = range(which(is.element(ids, disturb.ids)))

    ylim = c(min.score, 1.1*log10(max(norm[,-1])))
    plot.element=function(element, multi=F) {
        cluster = element.table$cluster[match(element, element.table$id)]
        size = q50$size[match(cluster, q50$cluster)]
        ccontigs = contigs$contig[contigs$cluster == cluster]

        x.median = q50[match(cluster, q50$cluster),-(1:2)]
        x.top = q95[match(cluster, q95$cluster),-(1:2)]
        x.top75 = q75[match(cluster, q75$cluster),-(1:2)]
        x.bottom25 = q25[match(cluster, q25$cluster),-(1:2)]
        x.bottom = q5[match(cluster, q5$cluster),-(1:2)]

        high = log10(x.top)
        low = log10(x.bottom)
        high.25 = log10(x.top75)
        low.25 = log10(x.bottom25)

        x = 1:M

        main=paste(element, " n=", size, sep="")

        plot.init(xlim=c(1,M), ylim=ylim,
                  main=main,
                  x.axis=F, y.axis=F, xaxs="i", yaxs="i")

        if (!multi)
            axis(side=2, las=2)

        abline(h=0, lty=3)

        color.5 = colors()[121]
        color.25 = colors()[563]
        #polygon(x=c(x,rev(x)), y=c(high, rev(low)), col=color.5, border=NA)
        #polygon(x=c(x,rev(x)), y=c(high.25, rev(low.25)), col=color.25, border=NA)
        for (contig in ccontigs) {
            coverage = log10(norm[norm$contig == contig,-1])
            coverage[coverage<min.score] = min.score
            lines(x=x, y=coverage, lwd=1, col="gray")
        }
        lines(x=x, y=log10(x.median), lwd=2)

        at = drange-0.5
        abline(v=at, lwd=2, lty=2)
        if (!multi) {
            segments(x0=x, y0=t(low), x1=x, y1=t(high))
            axis(side=1, labels=labels, at=1:M, las=2)
        }
    }

    Ny = ceiling(sqrt(N))
    Nx = ceiling(N/Ny)
    NN = c(1:N,rep(N+1,Nx*Ny-N))

    fig.start(fdir=fdir, ofn=paste(fdir, "/element_summary.pdf", sep=""), type="pdf", width=10, height=10)
    layout(matrix(NN, Nx, Ny, byrow=T))
    par(mai=c(0.05, 0.1, 0.2, 0.1))
    for (element in elements)
        plot.element(element, multi=T)
    fig.end()

    ffdir = paste(fdir, "/elements", sep="")
    for (element in elements) {
        fig.start(fdir=ffdir, ofn=paste(ffdir, "/", element, ".pdf", sep=""), type="pdf", width=6, height=3)
        plot.element(element, multi=F)
        fig.end()
    }
}
