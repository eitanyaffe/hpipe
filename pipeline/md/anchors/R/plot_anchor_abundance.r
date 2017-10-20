plot.anchor.abundance=function(ifn.coverage, ifn.ca, ifn.order, fdir)
{
    anchor.table = load.table(ifn.order)
    table.cov = load.table(ifn.coverage)
    ca = load.table(ifn.ca)

    table = data.frame(contig=table.cov$contig, a=table.cov$abundance.enrichment, l=table.cov$length)
    table$mass = table$l * 10^table$a
    table$mass = 100 * table$mass / sum(table$mass)

    N = 200 + 1
    breaks = seq(min(table$a), max(table$a), length.out=N)
    delta = breaks[2] - breaks[1]

    dense.f=function(table, total) {
        s = sapply(split(table$l, cut(table$a, breaks=breaks, include.lowest=T)), sum)
        result = data.frame(ab=breaks[-N], bp=s, cs=cumsum(s))
        # result$dd = (result$bp / total) / delta
        result$dd = (result$bp / total)
        result
    }

    xlim.all = range(table$a[is.element(table$contig, ca$contig)])

    plot.anchor=function(anchor.id, multi) {

        anchor = anchor.table$set[match(anchor.id, anchor.table$id)]

        ucontigs = ca$contig[ca$anchor == anchor]
        acontigs = ca$contig[ca$anchor == anchor & ca$contig_anchor == anchor]

        atable = table[is.element(table$contig, acontigs),]
        utable = table[is.element(table$contig, ucontigs),]

        total = sum(utable$l)
        df.u = dense.f(utable, total=total)
        df.a = dense.f(atable, total=total)

        xlim = range(df.u$a[df.u$dd != 0])
        ylim = range(c(df.a$dd, df.u$dd))

        plot.init(xlim=xlim.all, ylim=ylim, xlab="abundance", ylab="dd", axis.las=1, main=anchor.id, x.axis=!multi, y.axis=!multi)
        # abline(v=0, col=1)
        lines(x=df.u$a, y=df.u$dd, col="gray", lwd=4)
        lines(x=df.a$a, y=df.a$dd, col="red", lwd=2)
    }

    anchor.ids = anchor.table$id
    N = length(anchor.ids)
    Ny = ceiling(sqrt(N))
    Nx = ceiling(N/Ny)

    fig.start(fdir=fdir, ofn=paste(fdir, "/united.pdf", sep=""), type="pdf", width=8, height=8)
    par(mai=c(0.05, 0.05, 0.15, 0.05))
    layout(matrix(1:(Nx*Ny), Nx, Ny, byrow=T))
    for (anchor.id in anchor.ids) {
        plot.anchor(anchor.id=anchor.id, multi=T)
    }
    fig.end()

    ffdir = paste(fdir, "/anchors", sep="")
    for (anchor.id in anchor.ids) {
        fig.start(fdir=ffdir, ofn=paste(ffdir, "/", anchor.id, ".pdf", sep=""), type="pdf", width=8, height=4)
        plot.anchor(anchor.id=anchor.id, multi=F)
        fig.end()
    }

}
