plot.uniref.summary=function(ifn.summary, ifn.order, fdir)
{
    anchor.table = load.table(ifn.order)
    summary = load.table(ifn.summary)
    summary = summary[match(anchor.table$set, summary$anchor),]
    N = dim(summary)[1]

    # anchor	no_hit	no_taxa	masked	poor_hit	hit_70	hit_80	hit_90	hit_98	hit_100
    ll = list(
        none=list(fields=c("no_hit", "no_taxa", "masked", "poor_hit"), color="black"),
        hit.70=list(fields=c("hit_70"), color="darkblue"),
        hit.80=list(fields=c("hit_80"), color="blue"),
        hit.90=list(fields=c("hit_90"), color="red"),
        hit.97=list(fields=c("hit_98"), color="orange"),
        hit.100=list(fields=c("hit_100"), color="yellow"))
    df = data.frame(id=anchor.table$id[match(summary$anchor,anchor.table$set)])
    colors = NULL
    for (i in 1:length(ll)) {
        colors = c(colors, ll[[i]]$color)
        name = names(ll)[i]
        rr = rep(0,N)
        for (field in ll[[i]]$fields)
            rr = rr + summary[,field]
        df[,name] = rr
    }

    # genes
    fig.start(fdir=fdir, ofn=paste(fdir, "/gene_summary.pdf", sep=""), type="pdf", width=12, height=4)
    barplot(t(as.matrix(df[,-1])), col=colors, ylab="#genes", xlab="anchors", las=2, border=colors, names.arg=df$id)
    fig.end()
    wlegend(fdir=fdir, title="gene_summary", cols=colors, names=names(df)[-1])

    # fraction of genes
    fig.start(fdir=fdir, ofn=paste(fdir, "/gene_percent_summary.pdf", sep=""), type="pdf", width=12, height=4)
    barplot(t(as.matrix(100*df[,-1]/rowSums(df[,-1]))), col=colors, ylab="%genes", xlab="anchors", las=2, border=colors, names.arg=df$id)
    fig.end()

    cs = colSums(df[,-1])
    # o = order(cs, decreasing=T)
    # cs = cs[o]
    # colors = colors[o]

    css = round(100*cs/sum(cs), 2)
    main = sprintf(">70%% = %.1f%%", sum(css[-(1:2)]))
    fig.start(fdir=fdir, ofn=paste(fdir, "/united_summary.pdf", sep=""), type="pdf", width=4, height=4)
    pie(cs, col=colors, labels=NA, border=NA, main=main)
    fig.end()

    wlegend(fdir=fdir, title="united_counts", cols=colors, names=paste(colnames(cs), " N=", cs, ", P=", css, "%", sep=""))
}
