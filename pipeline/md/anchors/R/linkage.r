anchor.params=function(ifn.ca, ifn.cov, ofn)
{
    ca = load.table(ifn.ca)
    cov = load.table(ifn.cov)

    ca = ca[ca$anchor == ca$contig_anchor,]
    s = sapply(split(ca$enrichment, ca$anchor), median)
    df = data.frame(anchor=names(s), score=s)
    df$abundance = cov$median[match(df$anchor, cov$anchor)]

    save.table(df, ofn)
}

contig.linkage=function(ifn, ifn.anchors, ifn.cov, min.contig.coverage, min.anchor.contigs, ofn)
{
    table = load.table(ifn)
    atable = load.table(ifn.anchors)
    cov = load.table(ifn.cov)

    keep.ii = table$contig_coverage >= min.contig.coverage & table$anchor_contig_count >= min.anchor.contigs
    table = table[keep.ii,]

    df = data.frame(contig=table$contig, anchor=table$anchor, contig_anchor=table$contig_anchor,
        observed=table$contig_total_count, score=log10(table$contig_total_count / table$contig_expected))
    df$abundance = cov$abundance.enrichment[match(df$contig, cov$contig)]
    df$anchor.abundance = atable$abundance[match(df$anchor, atable$anchor)]
    df$anchor.score = atable$score[match(df$anchor, atable$anchor)]
    df$linkage = (df$score - df$anchor.score) + pmax((df$abundance - df$anchor.abundance),0)

    # table$factor = atable$abundance[match(table$anchor, atable$anchor)] / atable$coef[match(table$anchor, atable$anchor)]
    # table$linkage = (10^table$enrichment - 1) * table$factor
    # table$abundance = cov$abundance.enrichment[match(table$contig, cov$contig)]
    # table$observed = table$contig_total_count

    save.table(df, ofn)
}

plot.linkage=function(ifn, ifn.ca, ifn.anchors, fdir)
{
    table = load.table(ifn)
    atable = load.table(ifn.anchors)
    ca = load.table(ifn.ca)
    all.anchor.contigs = ca$contig[ca$anchor == ca$contig_anchor]

    df = field.count(ca, "contig")
    table$anchor.count = ifelse(is.na(match(table$contig, df$contig)), 0, df$count[match(table$contig, df$contig)])
    table = table[table$anchor.count>0,]

    anchors = atable$set

    # plot linkage vs observed
    ylim = range(table$linkage)
    xlim = c(0, log10(max(table$observed)))
    for (anchor in anchors) {
        aid = atable$id[match(anchor,atable$set)]
        ltable = table[table$anchor == anchor,]
        anchor.contigs = ca$contig[ca$anchor == anchor & ca$contig_anchor == anchor]
        anchored.contigs = ca$contig[ca$anchor == anchor]

        ltable$type =
            ifelse(is.element(ltable$contig, anchor.contigs), "anchor",
                   ifelse(is.element(ltable$contig, anchored.contigs) & ltable$anchor.count>1, "anchored.multi",
                          ifelse(is.element(ltable$contig, anchored.contigs), "anchored",
                                 ifelse(is.element(ltable$contig, all.anchor.contigs), "inter.anchor", "other"))))
        ltable$col =
            ifelse(ltable$type == "anchor", "red",
                   ifelse(ltable$type == "anchored.multi", "orange",
                          ifelse(ltable$type == "anchored", "black",
                                 ifelse(ltable$type == "inter.anchor", "blue", "gray"))))

                   fig.start(fdir=fdir, ofn=paste(fdir, "/", aid, ".pdf", sep=""), width=6, height=6, type="pdf")
        plot.init(xlim=xlim, ylim=ylim, ylab="linkage", xlab="observed", main=aid)
        points(log10(ltable$observed), ltable$linkage, col=ltable$col, pch=".",cex=ifelse(ltable$type == "extended.multi", 4, 2))
        fig.end()
    }
}

# old
anchor.factor=function(ifn.ca, ifn.cov, ofn)
{
    ca = load.table(ifn.ca)
    cov = load.table(ifn.cov)

    ca = ca[ca$anchor == ca$contig_anchor,]
    s = sapply(split(10^ca$enrichment, ca$anchor), median)
    df = data.frame(anchor=names(s), score=s)
    df$abundance = 10^cov$median[match(df$anchor, cov$anchor)]
    df$coef = (df$score - 1) * df$abundance

    save.table(df, ofn)
}

