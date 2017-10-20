
compute.trim.threshold=function(ifn, min.contacts=10, min.enrichment=1, fdr, ofn)
{
  data = load.table(ifn)

  data$obs = data$observed_contact_count
  data$exp = data$expected

  cells = unique(data$cell)
  contigs = unique(data$contig)
  T = length(cells) * length(contigs)

  data = data[data$cell != -1,]
  data = data[data$cell != data$other_cell,]
  data = data[data$exp != 0,]
  M = range(data$exp)
  N = sum(data$obs)

  exp = pmax(0,pmin(N,10^seq(log10(M[1]), log10(M[2])+min.enrichment+1, length.out=10^6)))
  obs = qbinom(p=(1-fdr), size=N, prob=exp/N, lower.tail=T)
  dup = duplicated(obs)
  exp = exp[!dup]
  obs = obs[!dup]

  # enrichment = log10(obs/exp)
  exp = pmin(exp, obs / 10^min.enrichment)

  a = approx(x=obs, y=exp, xout=min(obs):max(obs), rule=2)
  df = data.frame(obs=a$x, exp=a$y)
  save.table(df, ofn)
}

plot.trim.input=function(ifn, ifn.threshold, min.contacts=10, contigs.fn, fdir)
{
    contigs = load.table(contigs.fn)
    data = load.table(ifn)
    data$obs = data$observed_contact_count
    data$exp = data$expected
    data$col = ifelse(data$cell == data$other_cell, "red", "blue")
    data = data[data$obs > 0,]
    data$en = log10(data$obs/data$exp)
    data$type = ifelse(data$cell == data$other_cell, "intra", "inter")
    inter = data[data$type == "inter",]
    intra = data[data$type == "intra",]
    xlim = range(data$exp)
    ylim = range(data$obs)

    ################################################
    # plot anchor stats
    ################################################

    s = sapply(split(data$contig, data$cell), unique)
    df = data.frame(anchor=names(s), length=sapply(s, length), bp=sapply(s, function(x) { sum(contigs$length[is.element(contigs$contig, x)]) } ))
    fig.start(ofn=paste(fdir, "/anchors_stats_size.png", sep=""), fdir=fdir, width=800, height=400)
    barplot(df$length, names.arg=df$anchor, las=2, main="anchors #contigs", ylim=c(0,max(df$length*1.1)))
    fig.end()

    fig.start(ofn=paste(fdir, "/anchors_stats_bp.png", sep=""), fdir=fdir, width=800, height=400)
    barplot(df$bp, names.arg=df$anchor, las=2, main="anchors bp", ylim=c(0,max(df$bp*1.1)))
    fig.end()

    ################################################
    # plot inter density
    ################################################

    pdensity=function(table, field) {
        fin = is.finite(table[,field])
        inter.ind = table$type == "inter"
        intra.ind = table$type == "intra"

        d.inter = density(table[inter.ind, field])
        d.intra = density(table[intra.ind, field])

        v.inter = round(median(table[inter.ind, field]),3)
        v.intra = round(median(table[intra.ind, field]),3)

        xlim = range(c(d.inter$x, d.intra$x))
        ylim = range(c(d.inter$y, d.intra$y))

        fig.start(ofn=paste(fdir, "/", field, "_density.png", sep=""))
        main = paste(field, " inter=", v.inter, ", intra=", v.intra, sep="")
        plot(d.inter, xlab=field, lwd=2, xlim=xlim, ylim=ylim, col="black", main=main)
        lines(d.intra, col="red", lwd=2)
        legend("topleft", fill=c("black", "red"), legend=c("inter", "intra"))
        fig.end()
    }
    pdensity(data, "en")

    ################################################
    # plot anchor scatters
    ################################################

    df = load.table(ifn.threshold)
    df$obs = df$obs + 1

    fig.start(ofn=paste(fdir, "/scatter_threshold.png", sep=""), fdir=fdir, width=600, height=600)
    plot(data$exp, data$obs, xlab="expected", ylab="observed", pch=".", col=data$col, cex=4, log="xy", xlim=xlim, ylim=ylim)
    lines(df$exp, df$obs, lty=2)
    abline(a=0,b=1, lty=3)
    legend("topleft", fill=c("red", "blue"), legend=c("intra", "inter"))
    fig.end()

    ramp = colorRampPalette(c("white", "black"))
    ramp.intra = colorRampPalette(c("white", "red"))
    ramp.inter = colorRampPalette(c("white", "blue"))

    fig.start(ofn=paste(fdir, "/scatter_smooth_all.png", sep=""), fdir=fdir, width=600, height=600)
    smoothScatter(log10(data$exp), log10(data$obs), colramp=ramp, xlab="expected", ylab="observed", xlim=log10(xlim), ylim=log10(ylim))
    lines(log10(df$exp), log10(df$obs), lty=2)
    legend("topleft", fill=c("red", "blue"), legend=c("intra", "inter"))
    fig.end()

    fig.start(ofn=paste(fdir, "/scatter_smooth_inter.png", sep=""), fdir=fdir, width=600, height=600)
    smoothScatter(log10(inter$exp), log10(inter$obs), colramp=ramp.inter, xlab="expected", ylab="observed", xlim=log10(xlim), ylim=log10(ylim))
    lines(log10(df$exp), log10(df$obs), lty=2)
    fig.end()

    fig.start(ofn=paste(fdir, "/scatter_smooth_intra.png", sep=""), fdir=fdir, width=600, height=600)
    smoothScatter(log10(intra$exp), log10(intra$obs), colramp=ramp.intra, xlab="expected", ylab="observed", xlim=log10(xlim), ylim=log10(ylim))
    lines(log10(df$exp), log10(df$obs), lty=2)
    fig.end()

    fdir.anchors = paste(fdir, "/anchors", sep="")
    fig.dir(fdir.anchors)
    anchors = sort(unique(data$other_cell))
    for (anchor in anchors) {
        dd = data[data$other_cell == anchor,]
        fig.start(ofn=paste(fdir.anchors, "/", anchor, ".png", sep=""), width=600, height=600)
        plot(dd$exp, dd$obs, xlab="expected", ylab="observed", pch=".", col=dd$col, cex=4, log="xy", xlim=xlim, ylim=ylim)
        lines(df$exp, df$obs, lty=2)
        abline(a=0,b=1, lty=3)
        legend("topleft", fill=c("red", "blue"), legend=c("intra", "inter"))
        fig.end()
    }
}

plot.trim=function(pre.cc, post.cc, ifn.contigs, fdir, min.contacts, ifn.threshold)
{
    threshold.table = load.table(ifn.threshold)
    fig.dir(fdir)
    post = load.table(post.cc)
    pre = load.table(pre.cc)

    pre = pre[,c("contig", "cell", "other_cell", "observed_contact_count", "expected")]
    names(pre) = c("contig", "original_anchor", "anchor", "observed", "expected")

    contigs = load.table(ifn.contigs)

    result = NULL
    s = split(contigs, contigs$anchor)
    for (i in 1:length(s)) {
        anchor = names(s)[i]
        all = s[[i]]
        dis = sum(all$discarded)
        kept = dim(all)[1] - sum(all$discarded)
        cell = as.numeric(names(s)[i])
        result = rbind(result, data.frame(anchor=anchor, kept=kept, dis=dis))
        # cat(sprintf("anchor=%d kept=%d discarded=%d\n", dim(all)[1], kept, dis))
    }

    N = sum(contigs$discarded)
    M = dim(contigs)[1]

    ofn = paste(fdir, "/filter_stats.png", sep="")
    cat(sprintf("generating figure: %s\n", ofn))
    png(ofn, width=(200 + length(s)*15), height=400)
    barplot(t(as.matrix(result[,-1])), main=sprintf("multi filtered %d out of %d (%.1f%%)", N, M, 100*N/M), las=2, border=NA, col=c("black", "red"), ylab="#contigs")
    dev.off()

    dropped.anchors = result$anchor[result$kept == 0]

    clean=function(data) {
        data = data[data$original_anchor != -1,]
        data = data[data$observed > 0 & data$expected > 0,]
        data
    }
    pre = clean(pre)
    post = clean(post)
    # post$discarded = contigs$discarded[match(post$contig, contigs$contig)]
    # post = post[!post$discarded,]

    xlim = range(c(pre$observed, post$observed))
    ylim = range(c(pre$expected, post$expected))

    pscatter=function(table, title, ofn, main) {
        fig.start(ofn, width=600, height=600)
        plot(table$expected, table$observed, type="p", pch=".", cex=4, col=NA, xlim=ylim, ylim=xlim, main=main, log="xy")
        grid()
        T = dim(threshold.table)[1]
        lines(x=threshold.table$exp, y=threshold.table$obs, col="grey")
#        segments(x0=threshold.table$exp[-T], x1=threshold.table$exp[-1], y0=threshold.table$obs[-T], y1=threshold.table$obs[-T], col="grey")
#        segments(x0=threshold.table$exp[-1], x1=threshold.table$exp[-1], y0=threshold.table$obs[-T], y1=threshold.table$obs[-1], col="grey")
        abline(a=0, b=1, lty=3)
        points(table$expected, table$observed, pch=".", cex=4, col=table$col)
        fig.end()
    }

    mplot=function(data, title, fdir) {
        data$type = contigs$type[match(data$contig, contigs$contig)]
        data$col = ifelse(data$original_anchor == data$anchor,
            ifelse(data$type == "normal", "red", ifelse(data$type == "multi", "orange", "black")),
            ifelse(data$type == "normal", "blue", ifelse(data$type == "multi", "orange", "grey")))

        ofn = paste(fdir, "/", title, "_summary.png", sep="")
        pscatter(table=data, title=title, ofn=ofn, main="all")
        fdir.anchors = paste(fdir, "/", title, "_anchors", sep="")
        fig.dir(fdir.anchors)
        s = split(data, data$anchor)
        for (i in 1:length(s)) {
            tdata = s[[i]]
            anchor = names(s)[i]
            ofn = paste(fdir.anchors, "/", anchor, ".png", sep="")
            pscatter(table=tdata, title=title, ofn=ofn, main=anchor)
        }
    }

    post = post[!is.element(post$anchor, dropped.anchors) & !is.element(post$original_anchor, dropped.anchors),]
    mplot(pre, title="pre", fdir)
    mplot(post, title="post", fdir)
}

get.anchor.range=function(values, wt, factor)
{
    wt = wt / sum(wt)
    xmean = weighted.mean(values, wt)
    xsd = sqrt(sum(wt * (values - xmean)^2))
    from = xmean - xsd * factor
    to = xmean + xsd * factor

    ix = values>from & values<to
    in.range = round(100 * sum(wt[ix]),1)
    ssd = round(sd(values[ix]),3)

    list(ssd=ssd, in.range=in.range, from=from, to=to)
}

plot.anchor.coverage=function(ifn, ifn.coverage, factor, fdir)
{
    table = load.table(ifn)
    table = table[!table$discarded,]
    ctable = load.table(ifn.coverage)
    table$abun = ctable$abundance.enrichment[match(table$contig, ctable$contig)]
    table$length = ctable$length[match(table$contig, ctable$contig)]

    xlim = range(table$abun)

    anchors = sort(unique(table$anchor))
    for (anchor in anchors) {
        fig.start(fdir=fdir, ofn=paste(fdir, "/", anchor, ".png", sep=""))

        # compute range
        values = table$abun[table$anchor == anchor]
        wt = table$length[table$anchor == anchor]
        ar = get.anchor.range(values, wt, factor)
        plot(ecdf(values), main=paste(anchor, ", in.range=", ar$in.range, "%", ", selected.sd=", ar$ssd, sep=""), xlab="abundance", ylab="%", do.points=F, xlim=xlim)
        abline(v=c(ar$from,ar$to), lty=1, col="red")

        fig.end()
    }
}

select.contigs=function(ifn, ifn.coverage, factor, max.sd, ofn.contigs, ofn.lookup)
{
    table = load.table(ifn)
    table = table[!table$discarded, c("contig", "anchor")]

    # sort anchors by abundance
    cov = load.table(ifn.coverage)
    table$abun = cov$abundance.enrichment[match(table$contig, cov$contig)]
    table$length = cov$length[match(table$contig, cov$contig)]

    # trim by abundance
    anchors = sort(unique(table$anchor))
    result = NULL
    n.discarded.anchors = 0
    n.discarded.contigs = 0
    n.discarded.kbp = 0
    for (anchor in anchors) {
        atable = table[table$anchor == anchor,]
        ar = get.anchor.range(values=atable$abun, wt=atable$length, factor=factor)
        ix = atable$abun > ar$from & atable$abun < ar$to
        if (ar$ssd < max.sd) {
            result = rbind(result, atable[ix,])
            n.discarded.contigs = n.discarded.contigs + sum(!ix)
            n.discarded.kbp = n.discarded.kbp + sum(atable$length[!ix])
        } else {
            n.discarded.anchors = n.discarded.anchors + 1
        }
    }
    cat(sprintf("number of discarded anchors, due to potential non-linkage: %d\n", n.discarded.anchors))
    cat(sprintf("number of discarded contigs, abundance out of range: %d\n", n.discarded.contigs))
    cat(sprintf("discarded Mbp (not including discarded anchors): %.2f (%.2f%%)\n", n.discarded.kbp/10^6, n.discarded.kbp/sum(table$length)))

    ss = sapply(split(result$abun, result$anchor), median)
    dx = data.frame(anchor=names(ss), abundance=ss)
    dx = dx[order(dx$abundance, decreasing=T),]
    sorted.anchors = as.numeric(dx$anchor)
    result$new.anchor = match(result$anchor, sorted.anchors)
    df = data.frame(contig=result$contig, anchor=result$new.anchor, seed_anchor=result$anchor)

    save.table(df, ofn.contigs)

    dx$new.anchor = match(dx$anchor, sorted.anchors)
    save.table(dx[,c("new.anchor", "anchor")], ofn.lookup)
}
