
# scatter of CA contact enrichment vs CA temporal correlation
plot.contact.vs.temporal=function(ifn.contigs, ifn.ca, ifn.norm, ifn.anchors, ifn.matrix, min.contacts, fdir)
{
    contig.table = load.table(ifn.contigs)
    ca = load.table(ifn.ca)
    norm = load.table(ifn.norm)
    mean.response.table = load.table(ifn.anchors)
    mm = load.table(ifn.matrix)
    mm = mm[mm$contig_total_count >= min.contacts,]
    mm$enrichment = log10(mm$contig_total_count / mm$contig_expected)
    anchors = sort(unique(ca$anchor))

    ca$type = ifelse(ca$anchor == ca$contig_anchor, "anchor", "extended")
    xlim = c(-1,1)
    ylim = range(mm$enrichment)

    N = length(anchors)
    Ny = ceiling(sqrt(N))
    Nx = ceiling(N/Ny)
    NN = c(1:N,rep(N+1,Nx*Ny-N))

    # all together
    fig.start(ofn=paste(fdir, "/summary.pdf", sep=""), type="pdf", fdir=fdir, width=8, height=7)
    layout(matrix(NN, Nx, Ny))
    par(mai=c(0.05, 0.05, 0.15, 0.05))

    for (anchor in anchors) {
        mma = mm[mm$anchor == anchor,]
        contigs = mma$contig
        enrichment = mma$enrichment
        col = ifelse(mma$contig_anchor == anchor, "red", ifelse(mma$contig_anchor != 0, "blue", "black"))

        norm.anchor = mean.response.table[mean.response.table$anchor == anchor,-1]
        norm.contigs = norm[match(contigs, norm$contig), -1]

        cc = cor(t(norm.contigs), t(norm.anchor))

        # plot black on bottom on top
        ix = (col == "black")

        main = anchor
        plot.init(xlim=xlim, ylim=ylim, xlab="temporal correlation", ylab="contact enrichment", main=main, x.axis=F, y.axis=F)
        abline(h=0, v=0, lty=3)
        points(cc[ix], enrichment[ix], pch=".", col=col[ix], cex=2)
        points(cc[!ix], enrichment[!ix], pch=".", col=col[!ix], cex=2)
    }
    fig.end()
}

plot.anchor.majors=function(ifn, fdir)
{
    x = load.table(ifn)
    fig.start(fdir=fdir, ofn=paste(fdir, "/anchor_majors.png", sep=""), width=800)
    barplot(100*x$anchor.f, names.arg=x$anchor, las=2)
    title(main=paste("median=", median(100*x$anchor.f), sep=""))
    fig.end()
}

get.response.matrix=function(obs, ifn)
{
    library(inline)
    sig = signature(obs="numeric", n_contigs="integer", n_states="integer", result="numeric")
    body = readChar(ifn, file.info(ifn)$size)
    cfun = cfunction(sig, body=body, language="C", convention=".C")

    mat = as.matrix(obs - apply(obs, 1, mean))

    n.contigs = dim(mat)[1]
    n.states = dim(mat)[2]
    result = matrix(0, n.contigs, n.contigs)

    cat(sprintf("computing response correlation matrix...\n"))
    result[] = cfun(as.numeric(mat), n.contigs, n.states, as.numeric(result))$result
    result
}

plot.response.anchor.distrib=function(ifn.obs, ifn.ca, fdir)
{
    ca = load.table(ifn.ca)
    ca = ca[,c("contig", "anchor")]
    anchors = unique(ca$anchor)

    obs = load.table(ifn.obs)
    rownames(obs) = obs$contig

    ll = list()
    for (anchor in anchors) {
        acontigs = ca$contig[ca$anchor == anchor]
        obs.anchor = obs[is.element(obs$contig, acontigs),]
        cc = cor(t(obs.anchor[,-1]))
        cc[is.na(cc) | !is.finite(cc)] = -1
        hh = hclust(as.dist(1-cc), method="average")
        ll[[as.character(anchor)]] = 1-hh$height
        # sm = matrix2smatrix(cc)
        # sm = sm[sm$i != sm$j,]
        # ll[[as.character(anchor)]] = sm$value
    }
    ll = ll[order(sapply(ll, median))]
    fig.start(fdir=fdir, ofn=paste(fdir, "/anchor_response_distrib.png", sep=""), height=800, width=800)
    par(mai=c(2,1,1,0.1))
    boxplot(ll, las=2, outline=F)
    fig.end()
}

plot.response.matrix.anchor=function(ifn.obs, ifn.ca, ifn.cluster, ifn.order, ifn.majors, fdir)
{
    df = load.table(ifn.majors)
    cluster.table = load.table(ifn.cluster)

    breaks = c(-1,0,0.8,1)
    colors = c("blue", "white", "red", "orange")
    icolors = make.image.colors(colors=colors, breaks=breaks)

    contigs = load.table(ifn.order)$contig
    contigs = contigs[is.element(contigs, cluster.table$contig)]
    ca = load.table(ifn.ca)
    ca = ca[,c("contig", "anchor")]
    ca = ca[is.element(ca$contig, cluster.table$contig),]

    anchors = sort(unique(ca$anchor))

    obs = load.table(ifn.obs)
    rownames(obs) = obs$contig
    obs = obs[match(contigs, obs$contig),]

    for (anchor in anchors) {
        acontigs = ca$contig[ca$anchor == anchor]
        obs.anchor = obs[is.element(obs$contig, acontigs),]
        aclusters = cluster.table$cluster[match(obs.anchor$contig, cluster.table$contig)]
        cc = cor(t(obs.anchor[,-1]))
        cc[is.na(cc) | !is.finite(cc)] = -1
        N = length(acontigs)
        coords = 1:N

        mcluster = df$cluster[df$anchor == anchor]
        mportion = df$anchor.f[df$anchor == anchor]
        major.range = range(coords[aclusters == mcluster])

        fig.start(fdir=fdir, ofn=paste(fdir, "/", anchor, ".png", sep=""), height=800, width=800)
        plot.init(xlim=c(0,N), ylim=c(0,N), add.grid=F, add.box=F, x.axis=F, y.axis=F, main=paste(anchor, " fraction=", round(100*mportion, 1), "%", sep=""))
        image(coords, coords, cc, useRaster=T, col=icolors$col, breaks=icolors$breaks, add=T)
        rect(major.range[1], major.range[1], major.range[2]+1, major.range[2]+1)
        fig.end()
    }
}

plot.response.matrix=function(ifn.obs, ifn.ca, ifn.cluster, ifn.order, ifn.majors, fdir)
{
    debug = F

    breaks = c(-1,0,0.8,1)
    colors = c("blue", "white", "red", "orange")
    icolors = make.image.colors(colors=colors, breaks=breaks)
    panel = make.color.panel(colors=colors)
    wlegend2(fdir=fdir, panel=panel, breaks=breaks, title="response_matrix")

    ca = load.table(ifn.ca)
    ca = ca[,c("contig", "anchor")]

    cluster.table = load.table(ifn.cluster)
    contigs = load.table(ifn.order)$contig
    if (debug) contigs = contigs[1:400]
    N = length(contigs)

    cmap = field.count(ca, "contig")
    cmap$anchor = ifelse(cmap$count == 1, ca$anchor[match(cmap$contig, ca$contig)], 0)

    obs = load.table(ifn.obs)
    rownames(obs) = obs$contig
    obs = obs[match(contigs, obs$contig),]
    # mm = get.response.matrix(obs=obs, ifn=ifn.pearson.func)

    anchors = c(0, sort(unique(ca$anchor)))
    M = length(anchors)

    cat(sprintf("computing pearson over contigs: %d\n", length(contigs)))
    cc = cor(t(obs[,-1]))
    cc[is.na(cc) | !is.finite(cc)] = -1

    df = data.frame(
        contig=contigs,
        anchor=cmap$anchor[match(contigs, cmap$contig)],
        cluster=cluster.table$cluster[match(contigs, cluster.table$contig)])
    df$index = 1:N

    anchor.colors = c("white", rainbow(length(anchors)))
    df$anchor.color = anchor.colors[match(df$anchor, anchors)]

    axis.df = NULL
    for (anchor in anchors) {
        if (anchor == 0)
            next
        aclusters = df$cluster[df$anchor == anchor]
        aclusters = aclusters[aclusters != 0]
        if (length(aclusters) == 0)
            next
        mcluster = as.numeric(names(sort(table(aclusters), decreasing=T))[1])
        median.index = median(df$index[df$cluster == mcluster])
        min.index = min(df$index[df$cluster == mcluster])
        max.index = max(df$index[df$cluster == mcluster])
        axis.df = rbind(axis.df, data.frame(anchor=anchor, cluster=mcluster, start=min.index, end=max.index, center=median.index))
    }

    if (!debug) fig.start(fdir=fdir, ofn=paste(fdir, "/response_matrix.png", sep=""), height=1500, width=1500)

    bottom1 = -N/100
    bottom2 = -2*N/100
    bottom3 = -2.5*N/100
    bottom5 = -5*N/100
    plot.init(xlim=c(bottom5,N), ylim=c(bottom5,N), add.grid=F, add.box=F, x.axis=F, y.axis=F)
    image(1:N, 1:N, cc, useRaster=T, col=icolors$col, breaks=icolors$breaks, add=T)

    # majors
    rect(xleft=axis.df$start, xright=axis.df$end, ybottom=axis.df$start, ytop=axis.df$end, border=1, lwd=2)

    # color anchors below
    rect(ybottom=df$index, ytop=df$index+1, xleft=bottom1, xright=bottom2, col=df$anchor.color, border=df$anchor.color)

    # mark majors below
    rect(ybottom=axis.df$start, ytop=axis.df$end, xleft=bottom2, xright=bottom3, col=1, border=1)
    # segments(y0=c(axis.df$start, axis.df$end), y1=c(axis.df$start, axis.df$end), x0=bottom2, x1=bottom3, col=1)
    #abline(h=c(axis.df$start, axis.df$end), col="gray")
    text(y=axis.df$center, x=bottom3, labels=axis.df$anchor, pos=2, offset=0.7, cex=0.75)
    #text(y=axis.df$center, x=axis.df$center, labels=axis.df$anchor, cex=0.8, pos=2, col="white")

    fig.end()
}

plot.anchor.patterns=function(
    ifn.order, ifn.majors,
    ifn.median, ifn.top95, ifn.top75, ifn.bottom25, ifn.bottom05,
    ifn.taxa, ifn.taxa.legend, ifn.detection, labels, base.ids, disturb.ids, fdir)
{
    anchor.table = load.table(ifn.order)
    all.patterns = load.table(ifn.median)
    all.top95 = load.table(ifn.top95)
    all.top75 = load.table(ifn.top75)
    all.bottom25 = load.table(ifn.bottom25)
    all.bottom05 = load.table(ifn.bottom05)
    df = load.table(ifn.majors)

    min.score = log10(load.table(ifn.detection)[1,1])
    min.drop = -2

    xx = as.matrix(all.patterns[match(df$cluster, all.patterns$cluster),-(1:2)])
    detected = as.matrix(log10(xx) > min.score)
    base = rowSums(xx[,base.ids]) / length(base.ids)
    patterns = ifelse(detected, log10(xx / base), min.drop)
    M = dim(patterns)[2]

    dindex = which(is.element(colnames(patterns), disturb.ids))
    drange = c(min(dindex)-1, max(dindex))

    cc = cor(t(patterns))
    cc[is.na(cc)] = -1
    colnames(cc) = make.anchor.id(df$anchor, anchor.table)
    rownames(cc) = make.anchor.id(df$anchor, anchor.table)

    hh = hclust(as.dist(1-cc), method="average")
    df = df[hh$order,]
    cc = cc[hh$order,hh$order]
    patterns = patterns[hh$order,]
    base = base[hh$order]

    anchors = df$anchor
    sorted.anchors = anchor.table$set

    N = length(anchors)

    taxa.legend = load.table(ifn.taxa.legend)
    taxa.legend = taxa.legend[match(anchors, taxa.legend$anchor),]
    taxa.legend$text = taxa.legend$letter

    taxa = load.table(ifn.taxa)
    taxa = taxa[match(anchors, taxa$anchor),]

    ################################################################
    # plot clustering dendrogram
    ################################################################

    fig.start(fdir=fdir, ofn=paste(fdir, "/dendrogram.pdf", sep=""), type="pdf", height=4, width=1+M*0.75)
    plot(hh, hang=-1, ylab="1-r")
    fig.end()

    ################################################################
    # plot correlation matrix
    ################################################################

    colors = c("blue", "white", "red")
    breaks = c(-1, 0, 1)
    panel = make.color.panel(colors=colors)
    wlegend2(fdir=fdir, panel=panel, breaks=breaks, title="anchor_matrix")
    sm = matrix2smatrix(cc)
    sm$col = panel[vals.to.cols(sm$value, breaks)]

    fig.start(fdir=fdir, ofn=paste(fdir, "/anchor_matrix.pdf", sep=""), type="pdf", height=8, width=8)
    plot.new()
    plot.window(xlim=c(0,N), ylim=c(0,N))
    rect(sm$i-1, sm$j-1, sm$i, sm$j, col=sm$col, border=sm$col)
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/anchor_matrix_text.pdf", sep=""), type="pdf", height=16, width=16)
    plot.new()
    plot.window(xlim=c(0,N), ylim=c(0,N))
    rect(sm$i-1, sm$j-1, sm$i, sm$j, col=sm$col, border=sm$col)
    text(sm$i-0.5, sm$j-0.5, round(sm$value,2), cex=0.5)
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/pairwise_host_pearson_analysis.pdf", sep=""), type="pdf", height=6, width=6)
    smx = sm[sm$i != sm$j,]
    s = sapply(split(smx$value, smx$i), max)
    plot(ecdf(s), xlab="pearson", main="closest host (max pearson)")
    grid()
    fig.end()

    ################################################################
    # plot colored summary
    ################################################################

    major.colors = c("white", "darkgray")
    major.breaks = c(0, 1)
    major.panel = make.color.panel(colors=major.colors)
    wlegend2(fdir=fdir, panel=major.panel, breaks=major.breaks, title="major_anchor")

    base.colors = c("blue", "red", "orange")
    base.breaks = c(-1, 0, 1)
    base.panel = make.color.panel(colors=base.colors)
    wlegend2(fdir=fdir, panel=base.panel, breaks=base.breaks, title="base_anchor")

    colors = c("darkblue", "blue", "white", "red")
    breaks = c(-2, -1, 0, 1)
    panel = make.color.panel(colors=colors)
    wlegend2(fdir=fdir, panel=panel, breaks=breaks, title="anchor_trend")

    sm = matrix2smatrix(patterns)
    sm$col = ifelse(sm$value != min.drop, panel[vals.to.cols(sm$value, breaks)], "black")

    fig.start(fdir=fdir, type="pdf", ofn=paste(fdir, "/anchor_trend.pdf", sep=""), height=1+N*0.2, width=1+0.4*M)
    par(mai=c(2, 1, 1, 4))
    plot.new()
    plot.window(xlim=c(0,M+7), ylim=c(1,N))

    rect(sm$j-1+3, sm$i-1, sm$j+3, sm$i, col=sm$col, border=NA)
    rect(xleft=1.5, xright=2.5, ybottom=1:N-1, ytop=1:N, col=base.panel[vals.to.cols(log10(base), base.breaks)], border=NA)
    axis(side=1, labels="base", at=2, las=2)

    rect(xleft=M+3.5, xright=M+5.5, ybottom=1:N-1, ytop=1:N, border=NA, col=taxa.legend$color)
    text(x=M+4.5, y=1:N-0.5, labels=taxa.legend$text)
    axis(side=1, labels="taxa", at=M+4.5, las=2)

    text(x=M+5.7, y=1:N-0.5, adj=0, labels=make.anchor.id(anchors, anchor.table))
    axis(side=1, labels="anchor", at=M+6.5, las=2)

    abline(v=drange+3, lwd=2)
    axis(side=1, labels=labels, at=1:M-0.5+3, las=2)
    axis(side=4, labels=taxa$name, at=1:N-0.5, las=2)
    fig.end()

    ################################################################
    # plot detailed line plots
    ################################################################

    ylim = c(min.score, 1.5)
    plot.anchor=function(i, multi=F, sorted=F) {
        if (sorted)
            anchor = sorted.anchors[i]
        else
            anchor = anchors[i]
        cluster = df$cluster[match(anchor, df$anchor)]
        x.median = all.patterns[match(cluster, all.patterns$cluster),-(1:2)]
        x.top95 = all.top95[match(cluster, all.top95$cluster),-(1:2)]
        x.top75 = all.top75[match(cluster, all.top75$cluster),-(1:2)]
        x.bottom25 = all.bottom25[match(cluster, all.bottom25$cluster),-(1:2)]
        x.bottom05 = all.bottom05[match(cluster, all.bottom05$cluster),-(1:2)]

        high.5 = log10(x.top95)
        low.5 = log10(x.bottom05)
        high.25 = log10(x.top75)
        low.25 = log10(x.bottom25)

        main = make.anchor.id(anchor, anchor.table)
        x = 1:M
        plot.init(xlim=c(1,M), ylim=ylim,
                  main=main,
                  x.axis=F, y.axis=F, xaxs="i", yaxs="i")

        if (!multi)
            axis(side=2, las=2)

        abline(h=0, lty=3)

        color.5 = colors()[121]
        color.25 = colors()[563]
        polygon(x=c(x,rev(x)), y=c(high.5, rev(low.5)), col=color.5, border=NA)
        polygon(x=c(x,rev(x)), y=c(high.25, rev(low.25)), col=color.25, border=NA)
        lines(x=x, y=log10(x.median), lwd=2)
        at = drange+0.5
        abline(v=at, lwd=2, lty=2)
        if (!multi) {
            segments(x0=x, y0=t(low.5), x1=x, y1=t(high.5))
            axis(side=1, labels=labels, at=1:M, las=2)
        }
    }

    Ny = ceiling(sqrt(N))
    Nx = ceiling(N/Ny)
    NN = c(1:N,rep(N+1,Nx*Ny-N))

    # all together
    fig.start(fdir=fdir, ofn=paste(fdir, "/anchor_summary.pdf", sep=""), type="pdf", width=10, height=10)
    layout(matrix(NN, Nx, Ny, byrow=T))
    par(mai=c(0.05, 0.1, 0.2, 0.1))
    for (i in N:1)
        plot.anchor(i, multi=T)
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/anchor_summary_sorted.pdf", sep=""), type="pdf", width=10, height=10)
    layout(matrix(NN, Nx, Ny, byrow=T))
    par(mai=c(0.05, 0.1, 0.2, 0.1))
    for (i in 1:N)
        plot.anchor(i, multi=T, sorted=T)
    fig.end()

    ffdir = paste(fdir, "/anchors", sep="")
    for (i in 1:N) {
        anchor = df$anchor[i]
        fig.start(fdir=ffdir, ofn=paste(ffdir, "/", make.anchor.id(anchor, anchor.table), ".pdf", sep=""), type="pdf", width=6, height=3)
        plot.anchor(i, multi=F)
        fig.end()
    }
}

plot.response.vs.identity=function(ifn.majors, ifn.mean, ifn.identity, base.ids, fdir)
{
    all.patterns = load.table(ifn.mean)
    df = load.table(ifn.majors)
    anchors = df$anchor

    identity = load.table(ifn.identity)
    identity = identity[is.element(identity$set1, anchors) & is.element(identity$set2, anchors),]

    xx = as.matrix(all.patterns[match(df$cluster, all.patterns$cluster),-(1:2)])
    base = rowSums(xx[,base.ids]) / length(base.ids)
    patterns = log10(xx / base)
    M = dim(patterns)[2]
    cc = cor(t(patterns))
    cc[is.na(cc)] = -1

    hh = hclust(as.dist(1-cc), method="average")

    mcc = matrix2smatrix(cc)
    mcc$set1 = anchors[mcc$i]
    mcc$set2 = anchors[mcc$j]
    mcc = mcc[mcc$set1 != mcc$set2,]

    m = merge(mcc[,c("set1", "set2", "value")], identity[,c("set1", "set2", "identity")])

    # boxplots
#    s = split(m$value, cut(m$identity, breaks=c(0,40,60,80,100)))
    breaks = unique(quantile(m$identity, 0:10/10))
    s = split(m$value, cut(m$identity, breaks=breaks))
    fig.start(fdir=fdir, ofn=paste(fdir, "/response_vs_identity_boxplot.png", sep=""), height=300, width=100+length(breaks)*12)
    boxplot(s, col="gray", las=2, main="response rho \nvs % identity")
    fig.end()

    # scatter
    m = m[m$identity > 50,]
    fig.start(fdir=fdir, ofn=paste(fdir, "/response_vs_identity.png", sep=""), height=800, width=800)
    plot.init(xlim=range(m$identity), ylim=range(m$value), ylab="response pearson", xlab="% identity")
    plot(m$identity, m$value, pch="+")
    fig.end()
}

plot.anchor.control=function(ifn.order, ifn.contigs, ifn.mean, ifn.ca, ifn.norm, fdir)
{
    anchor.table = load.table(ifn.order)
    contig.table = load.table(ifn.contigs)
    ca = load.table(ifn.ca)
    norm = load.table(ifn.norm)
    patterns = load.table(ifn.mean)
    anchor.ids = anchor.table$id

    fc = field.count(ca, "contig")
    fc$shared = fc$count > 1

    breaks = c(-3, -2, 0.5, 0.8, 0.9, 0.99, 2)
    labels = c("shared", "<0.5", "0.5-0.8", "0.8-0.9", "0.9-0.99", ">0.99")
    colors = c("gray", "orange", "red", "black", "darkblue", "blue")
    wlegend(fdir=fdir, names=rev(labels), cols=rev(colors), title="response_matrix")

    compute.control=function(only.anchor) {
        result = NULL
        for (anchor.id in anchor.ids) {
            anchor = anchor.table$set[match(anchor.id, anchor.table$id)]
            ix = ca$anchor == anchor
            if (only.anchor)
                ix = ix & ca$contig_anchor != 0
            contigs = ca$contig[ix]
            norm.contigs = norm[match(contigs, norm$contig), -1]
            norm.anchor = patterns[patterns$cluster == anchor,-(1:2)]

            cc = cor(t(norm.contigs), t(norm.anchor))
            cc[is.na(cc)] = -1
            lengths = contig.table$length[match(contigs, contig.table$contig)]

            df = data.frame(cc=as.vector(cc), length=lengths, shared=fc$shared[match(contigs,fc$contig)])
            df$cc = ifelse(df$shared, -3, df$cc)
            x = sapply(split(df$length, cut(df$cc, breaks=breaks, include.lowest=T)), sum)
            x = 100*x/sum(x)
            result.anchor = data.frame(anchor.id=anchor.id, t(as.matrix(x)))
            colnames(result.anchor) = c("anchor.id", labels)
            result = rbind(result, result.anchor)
        }
        result
    }
    width = 10
    height = 5

    over.anchor = compute.control(only.anchor=T)
    fig.start(ofn=paste(fdir, "/control_over_anchor.pdf", sep=""), type="pdf", fdir=fdir, width=width, height=height)
    barplot(t(over.anchor[,-1]), col=colors, names.arg=anchor.ids, las=2, border=NA)
    fig.end()

    over.cloud = compute.control(only.anchor=F)
    fig.start(ofn=paste(fdir, "/control_over_cloud.pdf", sep=""), type="pdf", fdir=fdir, width=width, height=height)
    barplot(t(over.cloud[,-1]), col=colors, names.arg=anchor.ids, las=2, border=NA)
    fig.end()
}

plot.classes=function(ifn.median, ifn.class, fdir)
{
    system(paste("rm -rf", fdir))
    patterns = load.table(ifn.median)
    table = load.table(ifn.class)
    s = split(table$anchor, table$class)
    M = dim(patterns)[2]-2
    x = 1:M
    ylim = c(-3.5, 1.5)
    for (i in 1:length(s)) {
        anchors = s[[i]]

        fig.start(ofn=paste(fdir, "/", i, ".pdf", sep=""), fdir=fdir, type="pdf", width=4, height=4)

        plot.init(xlim=c(1,M), ylim=ylim, main=i, x.axis=T, y.axis=T)
        for (anchor in anchors) {
            x.median = patterns[match(anchor, patterns$cluster),-(1:2)]
            lines(x=x, y=log10(x.median), lwd=2)
        }
        fig.end()
    }
}

plot.anchors=function(rep.table, anchors, medians, min.score, main, drange, labels)
{
    M = dim(medians)[2] - 2
    ylim = c(min.score, 1.5)
    coords = 1:M

    names = paste(rep.table$anchor.id[match(anchors, rep.table$anchor)], ": ", rep.table$name[match(anchors, rep.table$anchor)], sep="")

    layout(matrix(1:2, 1, 2), widths=c(1,1.3))
    par(mai=c(0.6,0.6,0.4,0.05))
    plot.init(xlim=c(1,M), ylim=ylim,
              main=main,
              x.axis=F, y.axis=F, xaxs="i", yaxs="i")
    abline(h=0, lty=3)
    at = drange+0.5
    abline(v=at, lwd=2, lty=2)
    axis(side=1, labels=labels, at=1:M, las=2)
    axis(side=2, las=2)

    N = length(anchors)
    colors = if (N <= 3) 1:N+1 else rainbow(N)
    for (i in 1:N) {
        anchor = anchors[i]
        x.median = medians[match(anchor, medians$cluster),-(1:2)]
        lines(x=coords, y=log10(x.median), lwd=2, col=colors[i])
    }
    if (length(anchors) > 3) {
        group.median = apply(medians[match(anchors, medians$cluster),-(1:2)], 2, median)
        lines(x=coords, y=log10(group.median), lwd=4)
        colors = c(colors, 1)
        names = c(names, "median")
    }

    par(mai=c(0,0,0,0))
    plot.new()
    plot.window(0:1, 0:1)
    legend(x=0.05, y=0.95, fill=colors, legend=names, border=NA, bty="n")
}

plot.family.patterns=function(
    ifn.rep, ifn.majors, ifn.median,
    ifn.taxa, ifn.taxa.legend, ifn.detection, disturb.ids, labels, fdir)
{
    system(paste("rm -rf", fdir))
    rep.table = load.table(ifn.rep)
    medians = load.table(ifn.median)
    df = load.table(ifn.majors)
    taxa = load.table(ifn.taxa)
    min.score = log10(load.table(ifn.detection)[1,1])
    dindex = which(is.element(colnames(medians[,-(1:2)]), disturb.ids))
    drange = c(min(dindex)-1, max(dindex))

    taxa.legend = load.table(ifn.taxa.legend)
    taxa.legend$text = taxa.legend$letter
    s = split(taxa.legend$anchor, taxa.legend$group.id)

    for (i in 1:length(s)) {
        id = names(s)[i]
        name = taxa$name[match(id, taxa$tax_id)]
        level = taxa$level[match(id, taxa$tax_id)]
        anchors = s[[i]]
        fig.start(fdir=fdir, ofn=paste(fdir, "/", name, "_", level, ".pdf", sep=""), type="pdf", width=9, height=3)
        main = paste(name, " ", level, ", n=", length(anchors), sep="")
        plot.anchors(rep.table=rep.table, anchors=anchors, medians=medians, min.score=min.score, main=main, labels=labels, drange=drange)
        fig.end()
    }
}

plot.genus.patterns=function(
    ifn.rep, ifn.majors, ifn.median,
    ifn.taxa, ifn.taxa.legend, ifn.detection, disturb.ids, labels, fdir)
{
    system(paste("rm -rf", fdir))
    rep.table = load.table(ifn.rep)
    medians = load.table(ifn.median)
    df = load.table(ifn.majors)
    taxa = load.table(ifn.taxa)
    min.score = log10(load.table(ifn.detection)[1,1])
    dindex = which(is.element(colnames(medians[,-(1:2)]), disturb.ids))
    drange = c(min(dindex)-1, max(dindex))

    taxa.legend = load.table(ifn.taxa.legend)
    s = split(taxa.legend$anchor, taxa.legend$sub.group.id)

    for (i in 1:length(s)) {
        id = names(s)[i]
        anchors = s[[i]]
        name = taxa$name[match(id, taxa$tax_id)]
        level = taxa$level[match(id, taxa$tax_id)]
        fig.start(fdir=fdir, ofn=paste(fdir, "/", name, "_", level, ".pdf", sep=""), type="pdf", width=14, height=3)
        main = paste(name, " ", level, ", n=", length(anchors), sep="")
        plot.anchors(rep.table=rep.table, anchors=anchors, medians=medians, min.score=min.score, main=main, labels=labels, drange=drange)
        fig.end()
    }
}

plot.group.patterns=function(
    ifn.rep, ifn.majors, ifn.median,
    ifn.detection, ifn.groups, disturb.ids, labels, fdir)
{
    rep.table = load.table(ifn.rep)
    medians = load.table(ifn.median)
    df = load.table(ifn.majors)
    min.score = log10(load.table(ifn.detection)[1,1])
    group.table = load.table(ifn.groups)

    dindex = which(is.element(colnames(medians[,-(1:2)]), disturb.ids))
    drange = c(min(dindex)-1, max(dindex))

    s = split(group.table$anchor, group.table$group)
    for (i in 1:length(s)) {
        group = names(s)[i]
        anchors = s[[i]]
        fig.start(fdir=fdir, ofn=paste(fdir, "/g", group, ".pdf", sep=""), type="pdf", width=9, height=3)

        main = paste("g", group, ", n=", length(anchors), sep="")
        plot.anchors(rep.table=rep.table, anchors=anchors, medians=medians, min.score=min.score, main=main, labels=labels, drange=drange)
        fig.end()
    }
}
