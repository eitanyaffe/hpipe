detailed.inter.scatter.plots=function(ifn, ifn.ca, ifn.order, ifn.coverage,
    min.contacts, min.enrichment, min.contig.coverage, min.anchor.contigs, fdr, fdir)
{
    table = load.table(ifn)
    anchor.table = load.table(ifn.order)
    ca = load.table(ifn.ca)
    cov = load.table(ifn.coverage)

    table = table[table$contig_expected > 0,]

    # discard low quality contig/anchor pairs
    keep.ii = table$contig_coverage >= min.contig.coverage & table$anchor_contig_count >= min.anchor.contigs
    table = table[keep.ii,]

    df = field.count(ca, "contig")
    table$multi.anchor = ifelse(is.na(match(table$contig, df$contig)), F, df$count[match(table$contig, df$contig)] > 1)

    table$exp = table$contig_expected
    table$obs = table$contig_total_count
    table$type = ifelse(table$contig_anchor == 0,
        ifelse(table$multi.anchor, "extended.multi", "extended.single"),
        ifelse(table$contig_anchor == table$anchor, "intra", "inter"))

    table$enr = log10(table$obs / table$exp)

    anchors = anchor.table$set
    contigs = unique(table$contig)

    T = length(anchors) * length(contigs)
    afdr = fdr / T
    M = sum(table$obs[table$type == "inter"])

    # threshold line
    threshold.exp = 10^seq(log10(min(table$exp)), log10(max(table$exp)), length.out=100)
    threshold.obs = qbinom(p=(1-afdr), size=M, prob=threshold.exp/M, lower.tail=T)
    dup = duplicated(threshold.obs)
    threshold.exp = threshold.exp[!dup]
    threshold.obs = threshold.obs[!dup]
    threshold.obs = ifelse(threshold.obs < min.contacts, min.contacts,
        ifelse(log10(threshold.obs/threshold.exp) < min.enrichment, threshold.exp*10^min.enrichment, threshold.obs))
    threshold.enr = log10(threshold.obs / threshold.exp)

    # add max point
    threshold.obs = c(threshold.obs, max(table$obs))
    threshold.enr = c(threshold.enr, min.enrichment)

    table = table[table$contig_anchor != 0,]

    xlim = c(min.contacts, max(table$obs))
    xlim = c(min.anchor.contigs, max(table$obs))
    ylim = range(table$enr)

    N = length(anchors)
    fig.start(ofn=paste(fdir, "/inter_anchor_scatter.pdf", sep=""), type="pdf", fdir=fdir, width=1 + 1.2*N, height=1 + 1.2*N)
    par(mai=c(0.05, 0.05, 0.15, 0.05))
    layout(matrix(1:(N^2), N, N))
    for (src.anchor in anchors) {
    for (tgt.anchor in anchors) {
        main = sprintf("src=%s tgt=%s", make.anchor.id(src.anchor, anchor.table), make.anchor.id(tgt.anchor, anchor.table))

        ttable = table[table$anchor == tgt.anchor & table$contig_anchor == src.anchor,]
        plot.init(xlim=xlim, ylim=ylim, log="x", main=main, x.axis=F, y.axis=F, grid.lty=1)
        lines(x=threshold.obs, y=threshold.enr, lty=3)
        abline(h=0, lty=2)
        points(ttable$obs, ttable$enr, pch=".", cex=2)
    } }
    fig.end()
}
