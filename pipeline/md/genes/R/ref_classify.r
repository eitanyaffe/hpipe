plot.classify=function(ref.genes.ifn, genes.ifn, contigs.ifn, order.ifn, fdir)
{
    rgenes = load.table(ref.genes.ifn)

    genes = load.table(genes.ifn)
    ngenes = dim(genes)[1]
    contigs = load.table(contigs.ifn)

    refs = read.delim(order.ifn)$set
    nrefs = length(refs)

    fig.dir(fdir)

    p.types = c("match", "match.weak", "spurious", "switch", "no_ref")
    r.types = c("found", "lost", "weak")

    # gene summary
    tg = 100 * table(factor(genes$type, levels=p.types)) / ngenes
    ngene.chimerics = sum(genes$type == "switch")
    fig.start(ofn=paste(fdir, "/pre_genes_pie.png", sep=""))
    pie(tg, labels=paste(names(tg), ":", round(tg,1), "%", sep=""),
        main=paste("genes: n=", ngenes, ", switched=", ngene.chimerics, sep=""))
    fig.end()

    # gene summary, without match
    tg = 100 * table(factor(genes$type, levels=setdiff(p.types, "match"))) / ngenes
    ngene.chimerics = sum(genes$type == "switch")
    fig.start(ofn=paste(fdir, "/pre_genes_pie_no_match.png", sep=""))
    pie(tg, labels=paste(names(tg), ":", round(tg,1), "%", sep=""),
        main=paste("genes: n=", ngenes, ", switched=", ngene.chimerics, sep=""))
    fig.end()

    # type, by ref
    types = c("switch", "spurious")
    for (type in types) {
        xgenes = genes[genes$type == type, ]
        if (dim(xgenes)[1] == 0)
            next
        size = as.vector(table(factor(genes$ref, levels=refs)))
        tg = as.vector(table(factor(xgenes$ref, levels=refs), xgenes$type))
        x = 100 * tg / size
        ylim = c(0, max(x)*1.1)
        fig.start(ofn=paste(fdir, "/pre_genes_", type, ".png", sep=""), width=(nrefs*15 + 200))
        mp = barplot(x, names.arg=refs, ylab="%", col="grey", border=NA, las=2, main=type, ylim=ylim)
        text(x=mp, y=x, labels=tg, pos=3)
        fig.end()
    }

    ncontigs = dim(contigs)[1]

    # genes, by contig
    s = split(contigs,contigs$type)
    raw.tg = sapply(s, function(x) sum(x$total))
    tg = 100 * raw.tg / ngenes
    fig.start(ofn=paste(fdir, "/pre_genes_by_contig_pie.png", sep=""))
    pie(tg, labels=paste(names(tg), ":", round(tg,1), "%", sep=""),
        main=paste("(contig annotation)\ngenes: n=", ngenes, ", on_chimerics=", raw.tg[names(tg) == "chimeric"], sep=""))
    fig.end()


    # ref genes
    nref.genes = dim(rgenes)[1]
    tg = 100 * table(factor(rgenes$type, levels=r.types)) / nref.genes
    fig.start(ofn=paste(fdir, "/ref_genes_pie.png", sep=""))
    pie(tg, labels=paste(names(tg), ":", round(tg,1), "%", sep=""),
        main=paste("genes: n=", nref.genes, sep=""))
    fig.end()

    xgenes = rgenes[rgenes$type != "found", ]
    size = as.vector(table(factor(rgenes$ref, levels=refs)))
    tg = as.matrix(table(factor(xgenes$ref, levels=refs), xgenes$type))
    x = 100 * tg / size
    fig.start(ofn=paste(fdir, "/ref_genes.png", sep=""), width=(nrefs*15 + 200))
    barplot(t(x), names.arg=refs, col=c("darkblue", "darkgrey"), ylab="%", border=NA, las=2, main="ref genes fate", xlim=c(0,nrefs+5))
    legend("topright", fill=c("darkblue", "darkgrey"), legend=c("<70%", "<95%"))
    fig.end()
}
