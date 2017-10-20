plot.snp.density=function(ifn.anchors, ifn.ca, assembly.dir, dataset1, dataset2, fdir)
{
    table = load.table(ifn.anchors)
    ca = load.table(ifn.ca)
    contigs = ca$contig[ca$anchor == ca$contig_anchor]
    anchor.ids = table$id

    get.var=function(dataset) {
        ifn.var = paste(assembly.dir, "/datasets/", dataset, "/map_F_10_40/var_contig_summary_bins", sep="")
        var = load.table(ifn.var)
        var = var[is.element(var$contig, contigs),]
        s = split(var[,-1], ca$anchor[match(var$contig, ca$contig)])

        mm = t(as.matrix(sapply(s, function(x) { colSums(x[,-1]) / sum(x[,1]) }))) * 1000
        mm[match(names(s), table$set),]
    }
    mm1 = get.var(dataset1)
    mm2 = get.var(dataset2)

    N = length(anchor.ids)
    Ny = ceiling(sqrt(N))
    Nx = ceiling(N/Ny)
    NN = c(1:N,rep(N+1,Nx*Ny-N))

    ylim = c(0, 1.1*max(c(mm1,mm2)))

    plot.anchor=function(anchor.id, multi) {
        anchor = table$set[match(anchor.id, table$id)]
        mma1 = mm1[match(anchor.id,anchor.ids),]
        mma2 = mm2[match(anchor.id,anchor.ids),]
        mma = rbind(mma1, mma2)
        ylim = c(0, 1.1*max(mma))
        # barplot(mma, beside=T, ylim=ylim, border=NA, col="darkgreen", axisnames=F, axes=F)
        barplot(mma, beside=T, ylim=ylim, border=NA, axisnames=F, axes=F, col=c("blue", "red"), las=2)
        title(main=anchor.id)
        axis(1, labels=F, las=2)
        axis(2, labels=T, las=2)
    }

    # united anchor plot
    fig.start(fdir=fdir, ofn=paste(fdir, "/compare.pdf", sep=""), type="pdf", height=10, width=10)
    layout(matrix(NN, Nx, Ny, byrow=T))
    par(mai=c(0.05, 0.3, 0.15, 0.05))
    for (anchor.id in anchor.ids)
        plot.anchor(anchor.id=anchor.id, multi=T)
    fig.end()

    fdir = paste(fdir, "/anchors", sep="")
    for (anchor.id in anchor.ids) {
        fig.start(fdir=fdir, ofn=paste(fdir, "/", anchor.id, ".pdf", sep=""), type="pdf", height=4, width=4)
        par(mai=c(0.5, 0.5, 0.3, 0.05))
        plot.anchor(anchor.id=anchor.id, multi=F)
        fig.end()
    }
}
