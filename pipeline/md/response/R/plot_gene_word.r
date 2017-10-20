plot.gene.word=function(ifn, filter.words, min.count, fdir)
{
    df = load.table(ifn)
    df = df[!is.element(df$word, filter.words) & df$gene_count >= min.count,]
    df$enrichment = log2(df$enrichment)
    N = dim(df)[1]

    pos = rep(4, N)
    ss = c("Recombinase", "Recombination", "Tail", "Chromosome", "Hydrolase", "30S")
    pos[is.element(df$word, ss)] = 2
    pos[which.max(df$gene_count)] = 2

    xcolor = 0.75
    xlim = range(df$gene_count)
    xlim[1] = 2.5
    xlim[2] = xlim[2]*1.1
    ylim = range(df$enrichment)
    fig.start(fdir=fdir, ofn=paste(fdir, "/gene_diagram.pdf", sep=""), type="pdf", width=14, height=8)
    par(xaxs="i")
    par(yaxs="i")
    plot.init(xlim=xlim, ylim=ylim, xlab="gene count (log scale)", ylab="enrichment (log2)", log="x", axis.las=1)
    rect(xleft=0.01, xright=xlim[2]*100, ybottom=0, ytop=ylim[2]+1, col=rgb(1,xcolor,xcolor), border=NA)
    rect(xleft=0.01, xright=xlim[2]*100, ybottom=ylim[1]-1, ytop=0, col=rgb(xcolor,xcolor,1), border=NA)
    abline(h=0, col="lightgray")
    grid()
    box()
    points(df$gene_count, df$enrichment, col="lightgray", pch=19, cex=1.1)
    points(df$gene_count, df$enrichment, col="darkgray", cex=1.1)
    text(df$gene_count, df$enrichment, df$word, pos=pos, cex=1.2)
    fig.end()
}
