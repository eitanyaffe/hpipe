plot.assembly=function(ids, titles, idir, min.length, fdir, cols)
{
    N = length(ids)
    cols = cols[1:N]

    fig.start(fdir=fdir, width=600, height=200 + N*40, ofn=paste(fdir, "/id_legend.png", sep=""))
    par(mai=c(2,1.5,1,0.5))
    plot.new()
    legend("center", fill=cols, legend=ids)
    fig.end()

    contigs = NULL
    bps = NULL

    contigs.p = NULL
    bps.p = NULL

    breaks = c(1000, 5000, 10000, 20000, 40000, 100000, 200000, 500000)
    break.titles = paste("b", breaks, sep="_")

    for (id in ids) {
        x = load.table(paste(idir, "/", id, "/fold_0/contig_table", sep=""))

        df = data.frame(total=sum(x$length))
        for (i in 1:length(breaks))
            df[,break.titles[i]] = sum(x$length[x$length>breaks[i]])
        bps = rbind(bps, df)
        bps.p = rbind(bps.p, 100 * df[,-1]/df[1,1])

        df = data.frame(total=dim(x)[1])
        for (i in 1:length(breaks))
            df[,break.titles[i]] = sum(x$length>breaks[i])
        contigs = rbind(contigs, df)
        contigs.p = rbind(contigs.p, 100 * df[,-1]/df[1,1])
    }

    ylim.bpp = c(0, 1.1 * max(c(max(as.matrix(bps.p)), 70)))
    ylim.cp = c(0, 1.1 * max(c(max(as.matrix(contigs.p)), 12)))

    mm = max(bps/10^6)
    fig.start(fdir=fdir, width=200+N*50, ofn=paste(fdir, "/bp_count.png", sep=""))
    barplot(as.matrix(bps/10^6), beside=T, names.arg=c("all", paste(">", breaks/1000, "k", sep="")), col=cols, border=NA, main="nt (M)", las=2, ylim=c(0,mm*1.1))
    fig.end()

    fig.start(fdir=fdir, width=200+N*50, ofn=paste(fdir, "/bp_percentage.png", sep=""))
    barplot(as.matrix(bps.p), beside=T, names.arg=paste(">", breaks/1000, "k", sep=""), col=cols, border=NA, main="nt (%)", las=2, ylim=ylim.bpp)
    fig.end()

    fig.start(fdir=fdir, width=200+N*50, ofn=paste(fdir, "/contig_count.png", sep=""))
    barplot(as.matrix(contigs/1000), beside=T, names.arg=c("all", paste(">", breaks/1000, "k", sep="")), col=cols, border=NA, main="contigs (K)", las=2)
    fig.end()

    fig.start(fdir=fdir, width=200+N*50, ofn=paste(fdir, "/contig_percentage.png", sep=""))
    barplot(as.matrix(contigs.p), beside=T, names.arg=paste(">", breaks/1000, "k", sep=""), col=cols, border=NA, main="contigs (%)", las=2, ylim=ylim.cp)
    fig.end()
}
