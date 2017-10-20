plot.gene.project=function(fwd.source.ifn, bck.source.ifn, order.ifn1, order.ifn2, fdir)
{
    fwd.source = load.table(fwd.source.ifn)
    bck.source = load.table(bck.source.ifn)

    fwd.sets = load.table(order.ifn1)[,1]
    bck.sets = load.table(order.ifn2)[,1]

    types = c("found", "weak", "missing")
    tcols = c("lightgreen", "cyan", "red")

    fig.dir(fdir)

    psource=function(table, title, sets.order) {
        sets = sort(unique(table$set))
        N = length(sets)

        tt = table(table$set, table$type)
        df = data.frame(set=rownames(tt), found=tt[,1], weak=tt[,2], missing=tt[,3])
        df = df[match(sets.order, df$set),]

        fig.start(ofn=paste(fdir, "/breakdown_", title, ".png", sep=""), width=300 + N*30, height=600)
        par(mai=c(3,1,1,1))
        mp = barplot(t(1+tt), beside=T, las=2, col=tcols, log="y", plot=F)
        barplot(t(1+tt), beside=T, las=2, col=tcols, log="y", xlim=c(0, max(mp) * 1.2))
        grid()
        barplot(t(1+tt), beside=T, las=2, col=tcols, log="y", xlim=c(0, max(mp) * 1.2))
        legend("topright", fill=tcols, legend=types)
        fig.end()

        fig.start(ofn=paste(fdir, "/sensetivity_", title, ".png", sep=""), width=300 + N*30, height=600)
        par(mai=c(3,1,1,1))
        barplot(100 * tt[,1] / rowSums(tt), las=2, main="matching genes", ylab="%", ylim=c(0,100))
        fig.end()

    }
    psource(table=fwd.source, title="source_fwd", sets.order=fwd.sets)
    psource(table=bck.source, title="source_bck", sets.order=bck.sets)
}
