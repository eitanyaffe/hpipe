plot.response.matrix.anchor=function(ifn.order, ifn.norm, ifn.ca, fdir)
{
    system(paste("rm -rf", fdir))
    anchor.table = load.table(ifn.order)
    anchor.ids = anchor.table$id

    breaks = c(-1,0,0.8,1)
    colors = c("blue", "white", "red", "orange")
    icolors = make.image.colors(colors=colors, breaks=breaks)

    ca = load.table(ifn.ca)
    ca$is.anchor = ca$anchor == ca$contig_anchor
    ca = ca[,c("contig", "anchor", "is.anchor")]

    norm = load.table(ifn.norm)
    rownames(norm) = norm$contig

    plot.anchor=function(anchor.id, only.anchor)
    {
        anchor = anchor.table$set[match(anchor.id, anchor.table$id)]

        caa = ca[ca$anchor == anchor,]
        if(only.anchor)
            caa = caa[caa$is.anchor,]
        acontigs = caa$contig
        is.anchor = caa$is.anchor

        norm.anchor = norm[match(acontigs,norm$contig),]
        cc = cor(t(norm.anchor[,-1]))
        cc[is.na(cc) | !is.finite(cc)] = -1
        N = length(acontigs)
        coords = 1:N
        marg = max(4,N/20)

        hh = hclust(as.dist(1-cc), method="average")
        cc = cc[hh$order,hh$order]
        is.anchor = is.anchor[hh$order]

        if (!only.anchor)
            lim = c(-marg,N)
        else
            lim = c(0,N)
        plot.init(xlim=lim, ylim=lim, add.grid=F, add.box=F, x.axis=F, y.axis=F, main=anchor.id)
        image(coords, coords, cc, useRaster=T, col=icolors$col, breaks=icolors$breaks, add=T)
        if (!only.anchor) {
            rect(xleft=coords-0.5,xright=coords+0.5,ybottom=-marg,ytop=0.5, border=NA, col=ifelse(is.anchor,"black","white"))
            rect(ybottom=coords-0.5,ytop=coords+0.5,xleft=-marg,xright=0.5, border=NA, col=ifelse(is.anchor,"black","white"))
        }
    }

    N = length(anchor.ids)
    Ny = ceiling(sqrt(N))
    Nx = ceiling(N/Ny)
    NN = c(1:N,rep(N+1,Nx*Ny-N))

    # united anchor plot
    fig.start(fdir=fdir, ofn=paste(fdir, "/all_anchors.pdf", sep=""), type="pdf", height=10, width=10)
    layout(matrix(NN, Nx, Ny, byrow=T))
    par(mai=c(0.05, 0.05, 0.15, 0.05))
    for (anchor.id in anchor.ids)
        plot.anchor(anchor.id=anchor.id, only.anchor=T)
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/all_clouds.pdf", sep=""), type="pdf", height=10, width=10)
    layout(matrix(NN, Nx, Ny, byrow=T))
    par(mai=c(0.05, 0.05, 0.15, 0.05))
    for (anchor.id in anchor.ids)
        plot.anchor(anchor.id=anchor.id, only.anchor=F)
    fig.end()

    # separate plot per anchor
    for (anchor.id in anchor.ids) {
        fig.start(fdir=paste(fdir,"/anchor", sep=""), ofn=paste(fdir, "/anchor/", anchor.id, ".pdf", sep=""), type="pdf", height=4, width=4)
        par(mai=c(0.2, 0.2, 0.2, 0.2))
        plot.anchor(anchor.id=anchor.id, only.anchor=T)
        fig.end()

        fig.start(fdir=paste(fdir,"/cloud", sep=""), ofn=paste(fdir, "/cloud/", anchor.id, ".pdf", sep=""), type="pdf", height=4, width=4)
        par(mai=c(0.2, 0.2, 0.2, 0.2))
        plot.anchor(anchor.id=anchor.id, only.anchor=F)
        fig.end()
    }

}
