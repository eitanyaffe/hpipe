select.ca=function(table, min.contacts=10, min.enrichment, min.contig.coverage=0.5, min.anchor.contigs=2, fdr=0.0001)
{
  table$type = ifelse(table$contig_anchor == 0, "extended", ifelse(table$contig_anchor == table$anchor, "intra", "inter"))

  exp = table$contig_expected
  obs = table$contig_total_count
  N = sum(obs[table$type == "inter"])
  prob = exp / N
  obs = ifelse(obs < exp, exp, obs)
  table$pscore = pbinom(q=obs, size=N, prob=pmin(1,prob), lower.tail=F)
  table$enrichment = log10(table$contig_total_count/table$contig_expected)

  table$selected = (table$pscore < fdr &
                    table$contig_total_count >= min.contacts &
                    table$enrichment >= min.enrichment &
                    table$contig_coverage >= min.contig.coverage &
                    table$anchor_contig_count >= min.anchor.contigs)

  # remove anchor if no anchored are selected
  tt = table(table[table$anchor==table$contig_anchor & table$selected,"anchor"])
  ix = match(table$anchor,names(tt))
  table$anchor.count = ifelse(!is.na(ix), tt[ix], 0)
  table$selected = table$selected & table$anchor.count > 0

  table
}

anchor.contigs=function(ifn, min.contacts=10, min.enrichment, min.contig.coverage=0.5, min.anchor.contigs=2, fdr=0.0001, ofn)
{
    table = load.table(ifn)

    # remove contig/anchor pairs that had observed contacts
    table = table[table$any_observed,]

    anchors = unique(table$anchor)
    contigs = unique(table$contig)
    T = length(anchors) * length(contigs)
    fdr = fdr / T

    xtable = select.ca(table=table,
        min.enrichment=min.enrichment,
        min.contacts=min.contacts,
        min.contig.coverage=min.contig.coverage,
        min.anchor.contigs=min.anchor.contigs,
        fdr=fdr)

    table$enrichment = xtable$enrichment
    table = table[xtable$selected,]

    save.table(table, ofn)
}

plot.en=function(ifn, fdir)
{
    library(gplots)
    options(stringsAsFactors=F)
    table = load.table(ifn)

    table$noise = table$marginal - table$count
    table$col = ifelse(table$anchor == table$contig_anchor, "red", ifelse(table$contig_anchor != 0, "black", "gray"))
    anchors = unique(table$anchor)

    table = table[table$count > 0 & table$noise > 0,]
    xlim = range(table$noise)
    ylim = range(table$count)

    for (anchor in anchors) {
        ttable = table[table$anchor == anchor,]

        fig.start(ofn=paste(fdir, "/", anchor, ".png", sep=""), fdir=fdir, width=600, height=600)
        plot(ttable$noise, ttable$count, xlab="trans", ylab="cis", pch=".", col=ttable$col, cex=4, log="xy", xlim=xlim, ylim=ylim)
        legend("topleft", fill=c("red", "black", "grey"), legend=c("anchor", "trans", "other"))
        fig.end()
    }
}

plot.inter.anchor.matrix=function(ifn, ifn.order, min.contacts, min.enrichment, fdir)
{
    table = load.table(ifn)
    anchor.table = load.table(ifn.order)
    ii = table$anchor_src == table$anchor_tgt
    table$type = ifelse(ii, "intra", "inter")
    table$o[ii] = table$o[ii]/2
    table$e[ii] = table$e[ii]/2
    table$en = log10(table$o/table$e)
    table$en[!is.finite(table$en)] = 0
    anchors = sort(intersect(anchor.table$set, unique(table$anchor_src)))
    table$index1 = match(table$anchor_src, anchors)
    table$index2 = match(table$anchor_tgt, anchors)

    # sort by anchor id
    anchor.ids = anchor.table$id
    table$index.id1 = match(table$anchor_src, anchor.table$set)
    table$index.id2 = match(table$anchor_tgt, anchor.table$set)
    table$col = ifelse(ii, "red", "black")

    vtable = table[table$o > min.contacts & table$anchor_src <= table$anchor_tgt,]

    fig.start(ofn=paste(fdir, "/anchor_anchor_distrib.pdf", sep=""), fdir=fdir, width=8, height=8, type="pdf")
    plot.new()
    plot.window(xlim=range(vtable$o), ylim=range(vtable$en), log="x")
    title(xlab="obs", ylab="ca enrichment")
    grid()
    abline(h=0, lty=1)
    abline(h=min.enrichment, lty=2)
    points(vtable$o, vtable$en, col=vtable$col, pch="+")
    box()
    axis(1)
    axis(2)
    fig.end()

    wlegend(fdir=fdir, names=c("inter", "intra"), cols=1:2, title="inter_anchor")

    inter = table[table$anchor_src != table$anchor_tgt,]
    cc = cor(inter$o, inter$e)
    main = paste("cor=", round(cc, 3))

    mai = c(0.5, 0.5, 0.5, 0.5)

    table.n = table
#    table.n$en[table.n$o == 0] = -100
    colors = c("blue", "white", "red")
    breaks = c(-3, 0, 3)
    wheat(sm=table.n, field.x="index1", field.y="index2", field.value="en", mai=mai, main=main,
          fdir=fdir, ofn.prefix=paste(fdir, "/anchor_anchor_enrich", sep=""), width=8, height=8,
          plot.names=F, colors=colors, breaks=breaks, add.box=T,
          labels=anchors)
    wheat(sm=table.n, field.x="index.id1", field.y="index.id2", field.value="en", mai=mai, main=main,
          fdir=fdir, ofn.prefix=paste(fdir, "/id_anchor_anchor_enrich", sep=""), width=12, height=12,
          plot.names=T, colors=colors, breaks=breaks, add.box=T,
          labels=anchor.ids)
    wlegend2(fdir=fdir, panel=make.color.panel(colors=colors), breaks=breaks, title="matrix_enr")

 #   table = table[table$index1 != table$index2,]

#    colors = c("white", "blue", "red", "orange")
#    breaks = c(0, 10, 100, 1000)
    colors = c("white", "red")
    breaks = c(0, 1000)
    wheat(sm=table, field.x="index1", field.y="index2", field.value="o", mai=mai,
          fdir=fdir, ofn.prefix=paste(fdir, "/anchor_anchor_obs", sep=""), width=8, height=8,
          plot.names=F, colors=colors, breaks=breaks, add.box=T, labels=anchors)
    wlegend2(fdir=fdir, panel=make.color.panel(colors=colors), breaks=breaks, title="matrix_obs")

    wheat(sm=table, field.x="index1", field.y="index2", field.value="e", mai=mai,
          fdir=fdir, ofn.prefix=paste(fdir, "/anchor_anchor_exp", sep=""), width=8, height=8,
          plot.names=F, colors=colors, breaks=breaks, add.box=T, labels=anchors)
    wlegend2(fdir=fdir, panel=make.color.panel(colors=colors), breaks=breaks, title="matrix_exp")
}

plot.ca=function(ifn, fdir, min.contacts=10, min.enrichment, min.contig.coverage=0.5, min.anchor.contigs=2, fdr=0.0001)
{
    library(gplots)
    options(stringsAsFactors=F)
    table = load.table(ifn)
    table = table[table$contig_expected > 0 & table$any_observed,]

    # first inter scatter
    inter = table[table$contig_anchor != 0 & table$contig_anchor != table$anchor & table$anchor_contig_count > 0,]

    # discard low quality contig/anchor pairs
    keep.ii = table$contig_coverage >= min.contig.coverage & table$anchor_contig_count >= min.anchor.contigs
    table = table[keep.ii,]

    anchors = unique(table$anchor)
    contigs = unique(table$contig)
    T = length(anchors) * length(contigs)
    fdr = fdr / T

    table = select.ca(
        table=table,
        min.enrichment=min.enrichment,
        min.contacts=min.contacts,
        min.contig.coverage=min.contig.coverage,
        min.anchor.contigs=min.anchor.contigs,
        fdr=fdr)

    # check if contig is shared
    fc = field.count(table[table$selected,], field="contig")
    table$shared = ifelse(is.na(match(table$contig,fc$contig)),F,fc$count[match(table$contig,fc$contig)] > 1)

#    table$col =
#        ifelse(!table$selected, "gray",
#               ifelse(table$type == "intra", "red", ifelse(table$type == "extended", "black", "orange")))

    table$col =
        ifelse(table$selected,
               ifelse(table$type == "intra", "red", ifelse(!table$shared, "black", "orange")), "blue")

    color.order = c("black", "blue", "red", "orange")

    table$pch = 19
    table$cex = 0.4

    table$exp = table$contig_expected
    table$obs = table$contig_total_count

    N = sum(table$obs[table$type == "inter"])
    xlim = range(table$contig_expected)
    ylim = c(1, max(table$contig_total_count))

    # threshold line
    xs = 10^seq(log10(min(table$exp)), log10(max(table$exp)), length.out=100)
    ys = qbinom(p=(1-fdr), size=N, prob=pmin(1,xs/N), lower.tail=T)
    dup = duplicated(ys)
    xs = xs[!dup]
    ys = ys[!dup]
    ys = ifelse(ys < min.contacts, min.contacts, ifelse(log10(ys/xs) < min.enrichment, xs*10^min.enrichment, ys))

    pplot=function(x, xlim, ylim, main) {
        plot.new()
        plot.window(xlim=xlim, ylim=ylim ,log="xy")
        title(xlab="expected", ylab="observed", main=main)
        grid()
        box()
        axis(1)
        axis(2)
        for (i in 1:length(color.order)) {
            color = color.order[i]
            xx = x[x$col == color,]
            points(xx$exp, xx$contig_total_count, pch=xx$pch, col=xx$col, cex=xx$cex)
        }
    }

    ################################################################
    # scatter
    ################################################################

    fig.start(ofn=paste(fdir, "/all_scatter.png", sep=""), fdir=fdir, width=400, height=400)
    pplot(table, xlim=range(table$exp), ylim=range(table$contig_total_count), main="all")
    lines(xs, ys, lty=2)
    abline(a=0,b=1, lty=3)
    legend("topleft", fill=c("red", "black", "orange", "grey"), legend=c("anchor", "extended", "shared", "bg"))
    fig.end()

    fig.start(ofn=paste(fdir, "/inter_scatter.png", sep=""), fdir=fdir, width=600, height=600)
    cc = cor(inter$contig_expected, inter$contig_total_count)
    plot.init(xlab="exp", ylab="obs", xlim=range(inter$contig_expected),
              ylim=range(inter$contig_total_count),
              log="xy", main=paste("n=", dim(inter)[1], ", pearson:", round(cc,2)))
    abline(a=0, b=1, col="gray")
    points(inter$contig_expected, inter$contig_total_count, pch=".", cex=1)
    fig.end()

    ################################################################
    # density
    ################################################################

    pdensity=function(field, type="density", width, height) {

        fin = is.finite(table[,field])
        inter.ind = table$type == "inter" & fin & table$obs>=min.contacts
        intra.ind = table$type == "intra" & fin & table$obs>=min.contacts
        if (sum(inter.ind) < 3 || sum(intra.ind) < 3) {
            return (NULL)
        }


        v.inter = round(median(table[inter.ind, field]),3)
        v.intra = round(median(table[intra.ind, field]),3)

        fig.start(ofn=paste(fdir, "/", field, "_", type, ".png", sep=""), width=width, height=height)
        main = paste(field, " inter=", v.inter, ", intra=", v.intra, sep="")
        if (type == "density") {
            d.inter = density(table[inter.ind, field])
            d.intra = density(table[intra.ind, field])
            xlim = range(c(d.inter$x, d.intra$x))
            ylim = range(c(d.inter$y, d.intra$y))
            plot(d.inter, xlab=field, lwd=2, xlim=xlim, ylim=ylim, col="black", main=main)
            lines(d.intra, col="red", lwd=2)
        } else if (type == "boxplot") {
            ll = list(inter=table[inter.ind, field], intra=table[intra.ind, field])
            boxplot(ll, col=c("blue", "red"), las=2)
        } else stop()

        # legend("topleft", fill=c("black", "red"), legend=c("inter", "intra"))
        fig.end()
    }
    pdensity("enrichment", width=400, height=400)
    pdensity("enrichment", type="boxplot", width=150, height=400)

    ################################################################
    # compare inter to intra with boxplots
    ################################################################

    fin = is.finite(table$enrichment)
    inter.table = table[table$type == "inter" & fin & table$obs>=min.contacts,]
    intra.table = table[table$type == "intra" & fin & table$obs>=min.contacts,]
    ll = split(intra.table$enrichment, intra.table$anchor)
    oo = order(sapply(ll, median))
    lb = c(list(inter=inter.table$enrichment), ll[oo])
    fig.start(ofn=paste(fdir, "/anchors_inter_vs_intra.png", sep=""), width=800, height=300)
    boxplot(lb, col=c("gray", rep("white", length(ll))), outline=F, las=2)
    fig.end()

    scatter.fdir = paste(fdir, "/all_scatter", sep="")
    anchor.fdir = paste(fdir, "/all_anchor_contig_count", sep="")
    coverage.fdir = paste(fdir, "/all_contig_coverage", sep="")
    fig.dir(scatter.fdir)
    fig.dir(anchor.fdir)
    fig.dir(coverage.fdir)

    wlegend(fdir=fdir, names=c("anchor", "extended", "bg"), cols=c("red", "black", "grey"), title="o_e_scatter")

    anchors = sort(unique(table$anchor))
    for (anchor in anchors) {
        ttable = table[table$anchor == anchor,]
        if (dim(ttable)[1] == 0)
            next
        cat(sprintf("plotting anchor scatters: %d\n", anchor))
        fig.start(paste(scatter.fdir, "/", anchor, ".png", sep=""), width=600, height=600)
        pplot(ttable, main=anchor, xlim=xlim, ylim=ylim)
        # plot(ttable$exp, ttable$contig_total_count, col=ttable$col, pch=ttable$pch, cex=ttable$cex, log="xy", main=anchor, xlim=xlim, ylim=ylim)
        abline(a=0,b=1, lty=3)
        lines(xs, ys, lty=2)
        fig.end()

        fig.start(paste(anchor.fdir, "/", anchor, ".png", sep=""), width=600, height=600)
        plot(ttable$contig_total_count, ttable$anchor_contig_count, pch=".", type="p", log="xy", cex=4, xlim=ylim, ylim=ylim,
             xlab="contig_total_count", ylab="anchor_contig_count", col=ttable$col, main=anchor)
        fig.end()

        fig.start(paste(coverage.fdir, "/", anchor, ".png", sep=""), width=600, height=600)
        plot(ttable$contig_total_count, ttable$contig_coverage, pch=".", type="p", log="x", cex=4, xlim=ylim, ylim=c(0,1),
             xlab="contig_total_count", ylab="contig_coverage", col=ttable$col, main=anchor)
        dev.off()
    }
}

plot.cc.coverage=function(ifn, ifn.coverage, fdir)
{
    table = load.table(ifn)
    ctable = load.table(ifn.coverage)
    table$reads.per.bp = ctable$reads.per.bp[match(table$contig, ctable$contig)]
    table$length = ctable$length[match(table$contig, ctable$contig)]

    anchors = sort(unique(table$anchor))
    for (anchor in anchors) {
        cc = table[table$anchor == anchor,]
        ca = table[table$anchor == anchor & table$contig_anchor == anchor,]

        if (dim(ca)[1] > 0) {
            fig.start(fdir=paste(fdir, "/anchors", sep=""), ofn=paste(fdir, "/anchors/", anchor, ".png", sep=""))
            plot(ecdf(log10(ca$reads.per.bp*1000)), main=anchor, xlab="log10(reads/kbp)", ylab="%", do.points=F, xlim=c(2,5))
            grid()
            fig.end()
        }

        if (dim(cc)[1] > 0) {
            fig.start(fdir=paste(fdir, "/cc", sep=""), ofn=paste(fdir, "/cc/", anchor, ".png", sep=""))
            plot(ecdf(log10(cc$reads.per.bp*1000)), main=anchor, xlab="log10(reads/kbp)", ylab="%", do.points=F, xlim=c(2,5))
            fig.end()

            fig.start(fdir=paste(fdir, "/length_vs_coverage", sep=""), ofn=paste(fdir, "/length_vs_coverage/", anchor, ".png", sep=""))
            plot(log10(cc$reads.per.bp*1000), log10(cc$length), pch="+", main=anchor, xlim=c(2,6), ylim=c(3,6))
            grid()
            fig.end()

            fig.start(fdir=paste(fdir, "/score_vs_coverage", sep=""), ofn=paste(fdir, "/score_vs_coverage/", anchor, ".png", sep=""))
            plot(log10(cc$reads.per.bp*1000), cc$enrichment, pch="+", col=ifelse(cc$anchor == cc$contig_anchor, "red", "black"), main=anchor, xlim=c(2,6), ylim=c(0,6))
            grid()
            fig.end()
        }
    }
}
