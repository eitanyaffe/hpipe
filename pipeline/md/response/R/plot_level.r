plot.level.summary=function(ifn.elements, ifn.matrix, ifn.summary, fdir)
{
    elements = load.table(ifn.elements)
    mat = load.table(ifn.matrix)
    summary = load.table(ifn.summary)
    names = colnames(elements)[-1]

    # summary
    fig.start(fdir=fdir, ofn=paste(fdir, "/level_summary.pdf", sep=""), type="pdf", width=1+dim(summary)[1]*0.14, height=5)
    par(mai=c(2, 0.5, 0.5, 0.5))
    barplot(summary$count, names.arg=summary$level, col="darkblue", border=NA, las=2)
    fig.end()

    # cluster elements
    ## hh = hclust(dist(as.matrix(elements[,-1])))
    ## elements = elements[hh$order,]

    colors = c("white", "blue", "red", "orange")
    panel = make.color.panel(colors=colors)
    breaks = c(0, 1, 2, 3)
    wlegend(fdir=fdir, names=breaks, cols=colors, title="matrix_count")

    Nx = dim(elements)[2]-1
    Ny = dim(elements)[1]

    sm = matrix2smatrix(as.matrix(elements[,-1]))
    sm$col = panel[vals.to.cols(sm$value, breaks)]

    fig.start(fdir=fdir, ofn=paste(fdir, "/element_matrix.pdf", sep=""), type="pdf", width=1+Nx*0.1, height=1+Ny*0.14)
    par(mai=c(1.8, 0.25, 0.25, 0.25))
    plot.new()
    plot.window(xlim=c(0,Nx), ylim=c(0,Ny))
    rect(xleft=sm[,2]-1, xright=sm[,2], ybottom=Ny-sm[,1], ytop=Ny-sm[,1]+1, col=sm$col, border=NA)
    abline(v=0:Nx, h=0:Ny)
    mtext(side=1, at=1:Nx-0.5, names, las=2, cex=0.75)
    mtext(side=2, at=Ny:1-0.5, elements$index, las=2, cex=0.75)
    fig.end()

    # plot matrix
    colors = c("white", "blue", "red")
    panel = make.color.panel(colors=colors)
    breaks = c(0, 1, max(mat))
    N = dim(mat)[1]
    sm = matrix2smatrix(as.matrix(mat))
    sm$col = panel[vals.to.cols(sm$value, breaks)]
    sm = sm[sm[,1] < sm[,2],]
    ix = sm$value > 0
    fig.start(fdir=fdir, ofn=paste(fdir, "/level_matrix.pdf", sep=""), type="pdf", height=2+N*.25, width=2+N*.25)
    par(mai=c(2, 2, 0.5, 0.5))
    plot.new()
    plot.window(xlim=c(0,N), ylim=c(0,N))
    rect(xleft=sm[,1]-1, xright=sm[,1], ybottom=sm[,2]-1, ytop=sm[,2], col=sm$col, border=NA)
    text(x=sm[ix,1]-0.5, y=sm[ix,2]-0.5, sm$value[ix])
    abline(v=0:N, h=0:N)
    mtext(side=1, at=1:N-0.5, names, las=2)
    mtext(side=2, at=1:N-0.5, names, las=2)
    fig.end()
}
