get.mixmdl.values=function(table)
{
    table$value = log10(table$contacts / table$factor)
    z = table[table$factor > quantile(table$factor, 0.9) & table$contact > 0 & table$contig1 != table$contig2,]
    z$value
}

compute.mixmdl=function(table)
{
    library(mixtools)
    values = get.mixmdl.values(table)
    normalmixEM(values)
}

get.threshold=function(mdl)
{
    weights = 1 - (mdl$sigma / sum(mdl$sigma))
    weights = weights / sum(weights)
    threshold = sum(weights * mdl$mu)

    A = 1
    if (diff(mdl$mu) < A*sum(mdl$sigma)) {
        cat("Note: Two identified guasians are too close, cannot filter out noise at this early stage\n")
        threshold = -Inf
    }

    threshold
}

plot.filter.noise=function(ifn, ifn.mdl, fdir)
{
    library(mixtools)
    table = load.table(ifn)
    table$value = log10(table$contacts / table$factor)
    load(ifn.mdl)

    values = mixmdl$x
    threshold = get.threshold(mixmdl)

    system(paste("mkdir -p", fdir))

    ofn = paste(fdir, "/mixture_density.png", sep="")
    cat(sprintf("plotting figure: %s\n", ofn))
    png(ofn, 600, 400)
    plot(mixmdl, which=2)
    lines(density(values), lty=2, lwd=2)
    if (!is.infinite(threshold))
        abline(v=threshold, lty=2)
    dev.off()

    ofn = paste(fdir, "/mixture_scatter.png", sep="")
    cat(sprintf("plotting figure: %s\n", ofn))
    png(ofn, 600, 400)
    table = table[table$contact > 0 & table$contig1 != table$contig2,]
    table$noise = ifelse(table$value < threshold, T, F)
    table$col = ifelse(table$noise, "blue", "red")
    plot(table$contacts, table$factor, type="p", col=table$col, pch=".", log="xy", cex=2)
    dev.off()
}

filter.noise=function(ifn, ofn, ofn.filter.params, ofn.stats)
{
    cat(sprintf("reading file: %s\n", ifn))
    table = read.delim(ifn)
    table$value = log10(table$contacts / table$factor)
    mixmdl = compute.mixmdl(table)
    save(mixmdl, file=ofn.filter.params)

    threshold = get.threshold(mixmdl)
    noise = ifelse(table$value < threshold, T, F)

    cat(sprintf("filtering out %d/%d contig contacts (%.1f%%)\n", sum(noise), sum(!noise), 100*sum(noise)/dim(table)[1]))
    table = table[!noise,]
    cat(sprintf("saving file: %s\n", ofn))
    write.table(table, ofn, quote=F, col.names=T, row.names=F, sep="\t")

    table$type = ifelse(table$contig1 != table$contig2, "inter", "intra")

    stats = data.frame(
        intra.reads = sum(table$contacts[table$type == "intra"]),
        inter.masked.reads = sum(table$masked_contacts[table$type == "inter"]),
        inter.ok.reads = sum(table$contacts[table$type == "inter"]),
        inter.contig.pairs = sum(table$type == "inter"),
        inter.ok.contig.pairs = sum(table$contacts > 0 & table$type == "inter"))
    write.table(stats, ofn.stats, quote=F, col.names=T, row.names=F, sep="\t")
}

plot.contig.matrix=function(adir, fdir, aid, ids, titles)
{
    tt = NULL
    for (i in 1:length(ids)) {
        id = ids[i]
        title = titles[i]
        tt = rbind(tt, data.frame(id=id, load.table(paste(adir, "/datasets/", id, "/contig_mat/stats", sep=""))))
    }
    contigs = data.frame(ok=tt$inter.ok.contig.pairs, masked=tt$inter.contig.pairs - tt$inter.ok.contig.pairs)
    reads = data.frame(ok=tt$inter.ok.reads, masked=tt$inter.masked.reads)
    intra.reads = data.frame(intra=tt$intra.reads)
    inter.reads.P = data.frame(R=100 * tt$inter.ok.reads/tt$intra.reads)

    ccols = c("darkgreen", "orange")
    rcols = c("darkgreen", "orange")

    wbarplot(m=t(intra.reads/1000000), names=titles, main="intra reads", fdir=fdir, ofn=paste(fdir, "/intra_reads.png", sep=""), beside=F, normalize.columns=F, cols="grey", ylab="M reads")
    wbarplot(m=t(inter.reads.P), names=titles, main="inter/intra ratio (%)", fdir=fdir, ofn=paste(fdir, "/inter_over_intra_percentage.png", sep=""), beside=F, normalize.columns=F, cols="grey", ylab="%")

    wbarplot(m=t(reads/1000000), names=titles, main="inter reads", fdir=fdir, ofn=paste(fdir, "/reads.png", sep=""), beside=F, normalize.columns=F, cols=rcols, ylab="M reads")
    wbarplot(m=t(reads), names=titles, main="inter reads (%)", fdir=fdir, ofn=paste(fdir, "/reads_N.png", sep=""), beside=F, normalize.columns=T, cols=rcols, ylab="%")

    wbarplot(m=t(contigs/1000000), names=titles, main="contig pairs", fdir=fdir, ofn=paste(fdir, "/contigs.png", sep=""), beside=F, normalize.columns=F, cols=ccols, ylab="M pairs")
    wbarplot(m=t(contigs/1000000), names=titles, main="contig pairs", fdir=fdir, ofn=paste(fdir, "/contigs_N.png", sep=""), beside=F, normalize.columns=T, cols=ccols, ylab="%")

    wlegend(fdir=fdir, names=names(reads), cols=rcols, title="reads")
    wlegend(fdir=fdir, names=names(contigs), cols=ccols, title="contigs")

    save.table(tt, paste(fdir, "/contig_stats.txt", sep=""))
    save.table(inter.reads.P, paste(fdir, "/inter_read_percentage.txt", sep=""))


}

plot.contig.metric=function(ifn, fdir)
{
    x = load.table(ifn)
    x$cs = cumsum(x$count)
    x$csn = x$cs / sum(x$count)
    fig.start(fdir=fdir, ofn=paste(fdir, "/shared_cs.png", sep=""))
    plot(x$score, 100*x$csn, xlim=c(0,100), type="l", ylab="%", xlab="#shared neighbours", ylim=c(0,105))
    grid()
    abline(h=100, lty=2)
    abline(h=0, lty=2)
    fig.end()
}
