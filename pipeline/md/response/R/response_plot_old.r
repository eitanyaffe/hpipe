plot.anchor.majors=function(ifn, fdir)
{
    x = load.table(ifn)
    fig.start(fdir=fdir, ofn=paste(fdir, "/anchor_majors.png", sep=""), width=800)
    barplot(100*x$anchor.f, names.arg=x$anchor, las=2)
    title(main=paste("median=", median(100*x$anchor.f), sep=""))
    fig.end()
}

plot.hclust=function(ifn, ifn.ca, threshold, fdir)
{
    library("fastcluster")
    debug = F

    ca = load.table(ifn.ca)
    ca = ca[,c("contig", "anchor")]
    cmap = field.count(ca, "contig")
    cmap$anchor = ifelse(cmap$count == 1, ca$anchor[match(cmap$contig, ca$contig)], 0)

    table = load.table(ifn)
    rownames(table) = table$contig
    table = table[is.element(table$contig, ca$contig),]
    if (debug) table = table[1:400,]

    contigs = table$contig
    N = length(contigs)

    anchors = c(0, sort(unique(ca$anchor)))
    M = length(anchors)

    cat(sprintf("computing pearson over contigs: %d\n", length(contigs)))
    cc = cor(t(table[,-1]))
    cc[is.na(cc) | !is.finite(cc)] = -1

    cat(sprintf("computing hclust...\n"))
    dd = as.dist(1 - cc)
    hh = fastcluster::hclust(dd, method="average")

    contig.sets = list()
    anchor.count = NULL
    get.contigs=function(index) {
        if (index < 0)
            return(contigs[-index])
        else
            return(contig.sets[[as.character(index)]])
    }
    for (i in 1:(N-1)) {
        icontigs = unique(c(get.contigs(hh$merge[i,1]), get.contigs(hh$merge[i,2])))
        contig.sets[[as.character(i)]] = icontigs
        anchors = unique(cmap$anchor[match(icontigs, cmap$contig)])
        anchor.count = rbind(anchor.count, data.frame(index=i, value=length(anchors)))
    }

    flagged.nodes = NULL
    cluster.table = NULL
    cindex = 1
    for (i in (N-1):1) {
        i1 = hh$merge[i,1]
        i2 = hh$merge[i,2]
        flagged = is.element(i, flagged.nodes)
        cluster.root = hh$height < (1-threshold) && anchor.count$value[i] == 1 && !flagged
        if (cluster.root || flagged) {
            if (i1 > 0) flagged.nodes = c(flagged.nodes, i1)
            if (i2 > 0) flagged.nodes = c(flagged.nodes, i2)
        }
        if (cluster.root) {
            icontigs = get.contigs(i)
            cluster.table = rbind(cluster.table, data.frame(contig=icontigs, cluster=cindex))
            cindex = cindex + 1
        }
    }


    df = data.frame(contig=contigs, anchor=cmap$anchor[match(contigs, cmap$contig)], cut.cluster=cutree(hh, h=1-threshold))
    mx = match(contigs, cluster.table$contig)
    df$cluster = ifelse(!is.na(mx), cluster.table$cluster[mx], -1)

    df = df[hh$order,]
    df$index = 1:N
    # df = add.field.count(df, "cluster")
    # df$cluster[df$cluster_count == 1] = -1

    cluster.borders = c(which(diff(c(-10, df$cluster))!=0),N)
    anchor.borders = c(which(diff(c(-10, df$anchor))!=0),N)

    clusters = sort(unique(df$cluster))
    anchor.colors = c("white", sample(rainbow(length(anchors))))
    df$anchor.color = anchor.colors[match(df$anchor, anchors)]

    # cols = colors()[c(153, 201, 232, 250)]
    # cols = c("red", "purple", "blue", "gray")

    ## get.colors=function(ix) {
    ##     col1 = cols[bitwAnd(ix,3) + 1]
    ##     ix = bitwShiftR(ix, 2)
    ##     col2 = cols[bitwAnd(ix,3) + 1]
    ##     ix = bitwShiftR(ix, 2)
    ##     col3 = cols[bitwAnd(ix,3) + 1]
    ##     list(col1=col1, col2=col2, col3=col3)
    ## }
    ## ac = get.colors(match(df$anchor, anchors))
    ## df$col1 = ac$col1
    ## df$col2 = ac$col2
    ## df$col3 = ac$col3

    ## df$anchor.color1 = anchor.colors1[match(df$anchor, anchors)]
    ## anchor.colors2 = c("white", sample(rainbow(length(anchors))))
    ## df$anchor.color2 = anchor.colors2[match(df$anchor, anchors)]

    ## cluster.colors = sample(rainbow(length(clusters)))
    ## cluster.colors[clusters == -1] = "white"
    ## df$cluster.color = cluster.colors[match(df$cluster, clusters)]



    dg = as.dendrogram(hh)


    lwd = 1
    if (!debug) fig.start(fdir=fdir, ofn=paste(fdir, "/cluster_threshold_", threshold, ".png", sep=""), height=800, width=1800)
    plot(dg, yaxt="n", center=F, ylim=c(-0.2,attr(dg, "height")), leaflab="none")
    at = round(seq(0, attr(dg, "height"), by=0.1),2)
    axis(side=2, at=at, labels=1-at, las=2)
    segments(x0=df$index, x1=df$index, y0=-0.10, y1=-0.05, col=df$anchor.color, lwd=lwd)

    ## segments(x0=df$index, x1=df$index, y0=-0.10, y1=-0.05, col=df$col1, lwd=lwd)
    ## segments(x0=df$index, x1=df$index, y0=-0.15, y1=-0.1, col=df$col2, lwd=lwd)
    ## segments(x0=df$index, x1=df$index, y0=-0.2, y1=-0.15, col=df$col3, lwd=lwd)

    # segments(x0=df$index, x1=df$index, y0=-0.15, y1=-0.05, col=df$cluster.color, lwd=lwd)
    segments(x0=cluster.borders, x1=cluster.borders, y0=-0.05, y1=0, col="gray")

    vanchors = anchors[is.element(anchors, df$anchor)]
    # text(x=match(vanchors, df$anchor), y=-0.3, labels=vanchors)

    # grid()
    # abline(h=1-threshold, lty=2)
    # title(main=paste("threshold rho=", threshold, sep=""))
    if (!debug) fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/cluster_threshold_", threshold, "legend.png", sep=""), width=M*15+200, height=200)
    plot.init(ylim=c(-1,4), xlim=c(0,M), add.grid=F, add.box=F, x.axis=F, y.axis=F)
    rect(ybottom=0, ytop=1, xleft=1:M-1+0.1, xright=1:M-0.1, col=anchor.colors, border=1)
    text(y=0, pos=1, x=1:M-0.5, labels=anchors)
    fig.end()

    ## prects=function(y, col) {
    ##     rect(ybottom=y, ytop=y+1, xleft=1:M-1+0.1, xright=1:M-0.1, col=col, border=NA)
    ## }
    ## lac = get.colors(1:M)
    ## prects(y=0, col=lac$col1)
    ## prects(y=1, col=lac$col2)
    ## prects(y=2, col=lac$col3)
    ## rect(ybottom=0, ytop=3, xleft=1:M-1+0.1, xright=1:M-0.1, col=NA, border=1)
    ## text(y=0, pos=1, x=1:M-0.5, labels=anchors)
    ## fig.end()
}
