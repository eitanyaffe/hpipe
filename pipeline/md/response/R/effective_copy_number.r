plot.effective.copy.number=function(
    ifn.anchors, ifn.elements, ifn.norm,
    ifn.detection, disturb.ids, base.ids, base.min.correlation, labels, fdir)
{
    df = load.table(ifn.elements)
    anchor.table = load.table(ifn.anchors)
    norm = load.table(ifn.norm)
    norm = norm[,-2]
    min.score = log10(load.table(ifn.detection)[1,1])

    ids = colnames(norm)[-1]
    dindex = which(is.element(ids, disturb.ids))
    drange = c(min(dindex)-1, max(dindex))
    base.index = which(is.element(ids, base.ids))

#    df = df[df$host.count == 1,]
    df$anchor.id = df$hosts
    df$anchor = anchor.table$set[match(df$anchor.id, anchor.table$id)]

    get.response=function(clusters)
    {
        log10(norm[match(clusters, norm$cluster), -1])
    }

    cluster.response = get.response(df$cluster)
    anchor.response = get.response(df$anchor)

    # limit to clusters which were correlated with host before disturbance
    M = dim(df)[1]
    delta = (cluster.response[correlated,] - anchor.response[correlated,])

    mm = apply(delta, 2, mean)
    sd = apply(delta, 2, sd)
    middle = mm
    top = mm + sd
    bottom = mm - sd

    qnt = apply(delta, 2, function(x) { quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95))})
    bottom = qnt[1,]
    q25 = qnt[2,]
    middle = qnt[3,]
    q75 = qnt[4,]
    top = qnt[5,]

    N = length(labels)
    ylim = range(c(top, bottom))
    coord = 1:N-0.5
    width = 0.5
    xleft = coord-width/2
    xright = coord+width/2
    col = "lightblue"

    fig.start(fdir=fdir, ofn=paste(fdir, "/single_host_element_trend.pdf", sep=""), type="pdf", height=4, width=5)

    plot.init(xlim=c(0,N), ylim=ylim, x.axis=F, y.axis=F)
    grid()
    abline(h=0, col=1, lty=3)
    abline(v=drange, lwd=2, lty=2)

    # rect(xleft, bottom, xright, top, col=col, border=col)
    # segments(x0=xleft, x1=xright, y0=middle, y1=middle, col=1, lwd=2)

    segments(x0=coord, x1=coord, y0=bottom, y1=top, col=1, lwd=1)
    rect(xleft, q25, xright, q75, col=col, border=col)
    segments(x0=xleft, x1=xright, y0=middle, y1=middle, col=1, lwd=2)

    axis(2, las=2)
    axis(side=1, labels=labels, at=1:N-0.5, las=2)

    fig.end()
}
