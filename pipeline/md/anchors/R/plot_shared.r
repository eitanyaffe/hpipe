plot.shared=function(ifn, ifn.genes, ifn.contigs, fdir)
{
    ca = load.table(ifn)
    cc = load.table(ifn.contigs)
    genes = load.table(ifn.genes)

    df = field.count(ca, field="contig")
    df$length = cc$length[match(df$contig,cc$contig)]

    df.genes = field.count(genes, field="contig")
    ix = match(df$contig, df.genes$contig)
    df$gene.count = ifelse(!is.na(ix), df.genes$count[ix],0)

    # length of contigs
    s = round(sapply(split(df$length, df$count), sum)/1000)
    s = s[names(s) != "1"]
    if (length(s) == 0)
        return (NULL)

    ylim = c(0, 1.4*max(s))
    fig.start(fdir=fdir, ofn=paste(fdir, "/shared.pdf", sep=""), type="pdf", width=3, height=2)
    par(mai=c(0.5, 1, 0.5, 0.1))
    m = barplot(s, names.arg=names(s), col="darkblue", border=NA, ylim=ylim, las=1, ylab="kb")
    text(x=m, y=s, pos=3, labels=s)
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/shared_clean.pdf", sep=""), type="pdf", width=1.8, height=2)
    par(mai=c(0.5, 1, 0.5, 0.1))
    m = barplot(s, names.arg=names(s), col="darkblue", border=NA, ylim=ylim, las=1, ylab="kb")
    fig.end()

    # number of contigs
    s = round(sapply(split(df$contig, df$count), length))
    scount = sum(s[names(s) != "1"])
    if (scount == 0)
        return (NULL)

    s = s[names(s) != "1"]
    ylim = c(0, 1.4*(max(s)))
    width = 1 + 0.3*length(s)
    height = 4
    fig.start(fdir=fdir, ofn=paste(fdir, "/shared_count.pdf", sep=""), type="pdf", width=width, height=height)
    par(mai=c(0.5, 1, 0.5, 0.1))
    m = barplot(s, names.arg=names(s), col="darkblue", border=NA, ylim=ylim, las=1, ylab="#contigs")
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/shared_count_labels.pdf", sep=""), type="pdf", width=width, height=height)
    par(mai=c(0.5, 1, 0.5, 0.1))
    m = barplot(s, names.arg=names(s), col="darkblue", border=NA, ylim=ylim, las=1, ylab="#contigs")
    text(x=m, y=s, pos=3, labels=s, cex=0.75)
    fig.end()

    # number of genes
    s = round(sapply(split(df$gene.count, df$count), sum))
    s = s[names(s) != "1"]
    if (length(s) == 0)
        return (NULL)

    ylim = c(0, 1.4*max(s))
    fig.start(fdir=fdir, ofn=paste(fdir, "/genes_shared.pdf", sep=""), type="pdf", width=3, height=2)
    par(mai=c(0.5, 1, 0.5, 0.1))
    m = barplot(s, names.arg=names(s), col="darkblue", border=NA, ylim=ylim, las=1, ylab="kb")
    text(x=m, y=s, pos=3, labels=s)
    fig.end()
}
