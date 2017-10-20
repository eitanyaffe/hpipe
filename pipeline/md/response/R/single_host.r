plot.single.host.elements=function(lid, ifn.anchors, ifn.elements, ifn.norm,
    ifn.detection, ifn.classify, ifn.taxa.legend, labels, fdir)
{
    get.response.cluster=function(clusters)
    {
        log10(norm[match(clusters, norm$cluster), -1])
    }
    df = load.table(ifn.elements)
    anchor.table = load.table(ifn.anchors)
    norm = load.table(ifn.norm)
    norm = norm[,-2]
    min.score = log10(load.table(ifn.detection)[1,1])
    classify = load.table(ifn.classify)
    taxa.legend = load.table(ifn.taxa.legend)

    df = df[df$host.count == 1 & df$gene.count>2,]
    df$anchor.cluster = anchor.table$set[match(df$hosts, anchor.table$id)]
    m.element = get.response.cluster(df$cluster)
    m.host = get.response.cluster(df$anchor.cluster)

    eclass = data.frame(cluster=df$cluster,
        is.phage=classify$phage[match(df$cluster, classify$cluster)] == 1,
        is.mobile=rowSums(classify[match(df$cluster, classify$cluster),-1]))

    family = data.frame(cluster=df$cluster, color=taxa.legend$color[match(df$anchor, taxa.legend$anchor)])

    mm = m.element - m.host
    mm[m.element == min.score & m.host == min.score] = -30

    cc = cor(t(mm))
    hh = hclust(dist(1-cc))
    mm.order = mm[hh$order,]
    clusters = df$cluster[hh$order]

    eclass$y.coord = match(eclass$cluster, clusters)
    family$y.coord = match(family$cluster, clusters)

    colors = c("blue", "white", "red")
    breaks = c(-2, 0, 2)
    panel = make.color.panel(colors=colors)
    wlegend2(fdir=fdir, panel=panel, breaks=breaks, title="element_delta")

    colors.host = c("blue", "white", "red")
    breaks.host = c(-2, 0, 2)
    panel.host = make.color.panel(colors=colors.host)
    wlegend2(fdir=fdir, panel=panel, breaks=breaks.host, title="host")

    M = dim(norm)[2] - 1
    N = dim(df)[1]

    sm = matrix2smatrix(as.matrix(mm.order))
    sm$col = panel[vals.to.cols(sm$value, breaks)]
    sm$col[sm$val == -20] = "black"
    sm$col[sm$val == -30] = "gray"

    sm.host = matrix2smatrix(as.matrix(m.host[hh$order,]))
    sm.host$col = panel.host[vals.to.cols(sm.host$value, breaks.host)]

    fig.start(fdir=fdir, type="pdf", ofn=paste(fdir, "/element_over_host_trend.pdf", sep=""), height=1+N*0.02, width=12)
    layout(matrix(1:2,1,2))

    par(mai=c(1, 1, 1, 1))

    plot.new()
    plot.window(xlim=c(0,M), ylim=c(0,N))
    rect(sm$j-1, sm$i-1, sm$j, sm$i, col=sm$col, border=sm$col)

    plot.new()
    plot.window(xlim=c(0,M), ylim=c(0,N))
    rect(sm.host$j-1, sm.host$i-1, sm.host$j, sm.host$i, col=sm.host$col, border=sm.host$col)

    fig.end()
}

plot.single.host.element.detailed=function(
    lid, ifn.anchors, ifn.elements, ifn.norm,
    ifn.detection, disturb.ids, base.ids, base.min.correlation, labels, fdir)
{
    get.response.cluster=function(clusters)
    {
        log10(norm[match(clusters, norm$cluster), -1])
    }

    df = load.table(ifn.elements)
    anchor.table = load.table(ifn.anchors)
    norm = load.table(ifn.norm)
    norm = norm[,-2]
    min.score = log10(load.table(ifn.detection)[1,1])

    ids = colnames(norm)[-1]
    dindex = which(is.element(ids, disturb.ids))
    drange = c(min(dindex)-1, max(dindex))

    base.index = which(is.element(ids, base.ids))

    M = dim(norm)[2] - 1
    coords = 1:M
    ylim = c(-2, 1.5)
    ylim.diff = c(-2, 2)

    get.responses=function(anchor.id, filter.base) {
        anchor = anchor.table$set[match(anchor.id, anchor.table$id)]
        clusters = s[[anchor.id]]
        cluster.response = get.response.cluster(clusters)
        host.response = get.response.cluster(anchor)
        if (filter.base) {
            cc = cor(t(cluster.response[,base.index]), t(host.response[,base.index]))
            if (all(cc<base.min.correlation))
                return (NULL)
            cluster.response = cluster.response[cc>=base.min.correlation,]
        }
        list(host.response=host.response, cluster.response=cluster.response)
    }
    plot.response=function(anchor.id, multi=F, filter.base=T) {
        ll = get.responses(anchor.id=anchor.id, filter.base=filter.base)
        if (is.null(ll)) {
            plot.empty(title="no elements")
            return (NULL)
        }
        K = dim(ll$cluster.response)[1]
        plot.init(xlim=c(1,M), ylim=ylim, main=anchor.id, x.axis=F, y.axis=F)
        abline(h=0)
        at = drange+0.5
        abline(v=at, lwd=2, lty=2)
        if (K > 0)
            for (i in 1:K)
                lines(x=coords, y=ll$cluster.response[i,], lwd=1, col="gray")
        lines(x=coords, y=ll$host.response, lwd=2, col=1)

        if (!multi) {
            axis(side=1, labels=labels, at=1:M, las=2)
            axis(side=2, las=2)
        }
    }

    plot.response.diff=function(anchor.id, multi=F, filter.base=T) {
        ll = get.responses(anchor.id=anchor.id, filter.base=filter.base)
        if (is.null(ll)) {
            plot.empty(title="no elements")
            return (NULL)
        }
        K = dim(ll$cluster.response)[1]
        clusters = rownames(ll$cluster.response)
        delta.cluster.response = as.matrix(ll$cluster.response) - matrix(rep(unlist(ll$host.response), K), nrow=K, ncol=M, byrow=T)
        plot.init(xlim=c(1,M), ylim=ylim.diff, main=anchor.id, x.axis=F, y.axis=F)
        abline(h=0, lty=3)
        at = drange+0.5
        abline(v=at, lwd=2, lty=2)
        # abline(h=0, lty=3)
        colors = colorpanel(K, "orange", "purple")
        if (K > 0)
            for (i in 1:K)
                lines(x=coords, y=delta.cluster.response[i,], lwd=1, col=colors[i])

        if (!multi) {
            axis(side=1, labels=labels, at=1:M, las=2)
            axis(side=2, las=2)
            legend("topright", legend=clusters, fill=colors)
        }
    }

    df = df[df$host.count == 1,]
    s = split(df$cluster, df$hosts)
    anchor.ids = anchor.table$id[is.element(anchor.table$id,names(s))]

    # united plot
    N = length(anchor.ids)
    Ny = ceiling(sqrt(N))
    Nx = ceiling(N/Ny)
    NN = c(1:N,rep(N+1,Nx*Ny-N))

    ####################################
    # elements and host
    ####################################

    fig.start(fdir=fdir, ofn=paste(fdir, "/all.pdf", sep=""), type="pdf", height=10, width=10)
    layout(matrix(NN, Nx, Ny, byrow=T))
    par(mai=c(0.05, 0.05, 0.15, 0.05))
    for (anchor.id in anchor.ids)
        plot.response(anchor.id=anchor.id, multi=T)
    fig.end()

    # anchor plots
    ffdir = paste(fdir, "/anchors", sep="")
    for (anchor.id in anchor.ids) {
        fig.start(fdir=ffdir, ofn=paste(ffdir, "/", anchor.id, ".pdf", sep=""), type="pdf", height=3, width=6)
        plot.response(anchor.id=anchor.id, multi=F)
        fig.end()
    }

    ####################################
    # elements over host response
    ####################################

    fig.start(fdir=fdir, ofn=paste(fdir, "/all_diff.pdf", sep=""), type="pdf", height=10, width=10)
    layout(matrix(NN, Nx, Ny, byrow=T))
    par(mai=c(0.05, 0.05, 0.15, 0.05))
    for (anchor.id in anchor.ids)
        plot.response.diff(anchor.id=anchor.id, multi=T)
    fig.end()

    # anchor plots
    ffdir = paste(fdir, "/anchors_diff", sep="")
    for (anchor.id in anchor.ids) {
        fig.start(fdir=ffdir, ofn=paste(ffdir, "/", anchor.id, ".pdf", sep=""), type="pdf", height=3, width=6)
        plot.response.diff(anchor.id=anchor.id, multi=F)
        fig.end()
    }
}
