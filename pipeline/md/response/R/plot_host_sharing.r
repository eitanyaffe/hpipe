rotate.rects=function(sm)
{
    rotate=function(x,y) { list(x=(x+y)/2, y=(y-x)/2) }
    p1 = rotate(x=sm$j-1, y=sm$i-1)
    p2 = rotate(x=sm$j, y=sm$i-1)
    p3 = rotate(x=sm$j, y=sm$i)
    p4 = rotate(x=sm$j-1, y=sm$i)
    x = as.vector(t(cbind(p1$x, p2$x, p3$x, p4$x, NA)))
    y = as.vector(t(cbind(p1$y, p2$y, p3$y, p4$y, NA)))
    list(x=x, y=y)
}

plot.host.sharing=function(ifn.anchors, ifn.network, fdir)
{
    network = load.table(ifn.network)
    anchor.table = load.table(ifn.anchors)
    anchor.ids = anchor.table$id
    N = length(anchor.ids)

    network$anchor.id = anchor.table$id[match(network$anchor, anchor.table$set)]
    network$anchor.index = match(network$anchor.id, anchor.ids)
    ss = split(network$anchor.index, network$cluster)

    mm = matrix(0, N, N)
    for (i in 1:length(ss)) {
        gg = expand.grid(ss[[i]], ss[[i]])
        gg = gg[gg[,1] != gg[,2],]
        for (j in 1:dim(gg)[1])
            mm[gg[j,1], gg[j,2]] = mm[gg[j,1], gg[j,2]] + 1
    }
    diag(mm) = 0.5
    aa = matrix2smatrix(mm)
    aa$col = ifelse(aa$value == 0.5, "gray", ifelse(aa$value == 1, "red", "orange"))

    raa.bg = rotate.rects(aa)
    aa = aa[aa$value>0,]
    raa = rotate.rects(aa)

    fig.start(fdir=fdir, ofn=paste(fdir, "/anchor_sharing.pdf", sep=""), type="pdf", width=10, height=10/sqrt(2))
    plot.new()
    par(yaxs="i")
    par(mai=c(0.5, 0.5, 0, 0))
    plot.window(xlim=c(0,N), ylim=c(0,N/2))
    polygon(x=raa.bg$x, y=raa.bg$y, border="lightgray", col=NA)
    polygon(x=raa$x, y=raa$y, border=aa$col, col=aa$col)
    mtext(text=anchor.ids, side=1, at=1:N-0.5, las=2, line=0.5)

    fig.end()
}
