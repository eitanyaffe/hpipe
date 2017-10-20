plot.sample.matrix=function(ifn.order, ifn.response, ifn.contigs, ifn.median, labels, disturb.ids, fdir)
{
    anchor.table = load.table(ifn.order)
    contig.table = load.table(ifn.contigs)
    response = load.table(ifn.response)
    min.abundance = min(response[,-1])
    method = "pearson"

    # disturbance
    ids = colnames(response)[-1]
    dindex = which(is.element(ids, disturb.ids))
    drange = c(min(dindex)-1, max(dindex))

    patterns = load.table(ifn.median)
    patterns = patterns[,-2]
    patterns = patterns[is.element(patterns$cluster, anchor.table$set),]

    colors=c("blue", "white", "orange", "darkgreen")
    breaks=c(-1, 0, 0.5, 1)
    panel = make.color.panel(colors=colors)
    wlegend2(fdir=fdir, panel=panel, breaks=breaks, title="pearson")

    colors.limited=c("white", "orange", "darkgreen")
    breaks.limited=c(0, 0.5, 1)
    panel.limited = make.color.panel(colors=colors.limited)
    wlegend2(fdir=fdir, panel=panel.limited, breaks=breaks.limited, title="limited_pearson")

    M = length(ids)
    lim = c(0, M)

    par(xaxs="i")
    par(yaxs="i")
    plot.matrix=function(mat, fdir, title, multi=F, add.text=F, cex, check.detected=F) {
        tcc = cor(mat, method=method)
        tcc[is.na(tcc)] = -2

        if (check.detected) {
            detected = apply(mat, 2, median) > min.abundance
            tcc[!detected,] = -2
            tcc[,!detected] = -2
        }

        sm = matrix2smatrix(tcc)
        sm$x = sm[,1]
        sm$y = sm[,2]
        sm = sm[sm$i != sm$j,]
        sm$cc = ifelse(sm$value != -2, panel[vals.to.cols(vals=sm$value, breaks=breaks)], "gray")

        plot.new()
        plot.window(xlim=lim, ylim=lim)
        title(main=title)
        rect(sm$x-1,sm$y-1,sm$x,sm$y,border=NA,col=sm$cc)
        if (!multi) {
            abline(v=drange, lwd=2)
            abline(h=drange, lwd=2)
            box()
            axis(1, at=1:M - 0.5, labels=labels, las=2)
            axis(2, at=1:M - 0.5, labels=labels, las=2)
        }
        if (add.text) {
            text(sm$x-0.5,sm$y-0.5, round(sm$value,2), cex=cex)
        }
    }

    height = 6
    width = 5

    fig.start(ofn=paste(fdir, "/all_clusters.pdf", sep=""), type="pdf", fdir=fdir, width=width, height=height)
    plot.matrix(mat=patterns[,-1], fdir=fdir, title="all", check.detected=F)
    fig.end()

    fig.start(ofn=paste(fdir, "/all_clusters_text.pdf", sep=""), type="pdf", fdir=fdir, width=width*2, height=height*2)
    plot.matrix(mat=patterns[,-1], fdir=fdir, title="all", add.text=T, cex=0.75, check.detected=F)
    fig.end()

    anchor.ids = anchor.table$id
    ffdir = paste(fdir, "/anchors", sep="")
    for (anchor.id in anchor.ids) {
        anchor = anchor.table$set[match(anchor.id, anchor.table$id)]
        contigs = contig.table$contig[contig.table$anchor == anchor]
        mat = response[is.element(response$contig,contigs),-1]

        title = paste(anchor.id, round(median(cor(t(mat), method=method)),2))
        fig.start(ofn=paste(ffdir, "/", anchor.id, ".pdf", sep=""), type="pdf", fdir=ffdir, width=width, height=height)
        plot.matrix(mat=mat, fdir=ffdir, title=title, check.detected=T)
        fig.end()
    }

    N = length(anchor.ids)
    Ny = ceiling(sqrt(N))
    Nx = ceiling(N/Ny)
    NN = c(1:N,rep(N+1,Nx*Ny-N))
    fig.start(ofn=paste(fdir, "/all_anchors.pdf", sep=""), type="pdf", fdir=fdir, width=12, height=12)
    par(mai=c(0.1, 0.2, 0.2, 0.1))
    layout(matrix(NN, Nx, Ny, byrow=T))
    for (anchor.id in anchor.ids) {
        anchor = anchor.table$set[match(anchor.id, anchor.table$id)]
        contigs = contig.table$contig[contig.table$anchor == anchor]
        mat = response[is.element(response$contig,contigs),-1]
        title = paste(anchor.id, round(median(cor(t(mat), method=method)),2))
        plot.matrix(mat=mat, fdir=ffdir, title=title, multi=T, check.detected=T)
    }
    fig.end()
}
