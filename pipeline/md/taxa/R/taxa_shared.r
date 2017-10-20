matrix2smatrix=function(matrix)
{
  dim1 = dim(matrix)[1]
  dim2 = dim(matrix)[2]
  v = as.vector(matrix)
  indices = 1:(dim1*dim2)
  i = (indices-1) %% dim1 + 1
  j = floor((indices-1) / dim1) + 1
  data.frame(i=i, j=j, value=v, stringsAsFactors=F)
}

make.color.panel=function(colors, ncols=256)
{
  panel = NULL
  for (i in 2:length(colors))
    panel = c(panel, colorpanel(ncols, colors[i-1], colors[i]))
  panel
}

# for example vals.to.cols(1:10, c(1, 3, 10), ncols=10) returns:
# [1] 1  6 11 12 14 15 16 17 19 20
vals.to.cols=function(vals, breaks, ncols=256)
{
  min = breaks[1]
  max = breaks[length(breaks)]
  vals = ifelse(vals < min, min, ifelse(vals>max, max, vals))
  n = length(breaks)-1
  cols = rep(-1, length(vals))
  for (i in 1:n)
  {
    ind = (breaks[i] <= vals) & (vals <= breaks[i+1])
    if (!any(ind))
      next
    # normalize to [0,1]
    cols[ind] = (vals[ind] - breaks[i]) / (breaks[i+1] - breaks[i])
    # normalize to [i*ncols,i*(ncols+1)]
    cols[ind] = (i-1)*ncols + cols[ind]*(ncols-1) + 1
    # round
    cols[ind] = round(cols[ind])
  }
  return (cols)
}

shared.matrix=function(ifn, ifn.rep, order.dir, fdir)
{
    p=function(names, ofn, size) {
        cat(sprintf("plotting file: %s\n", ofn))
        png(ofn, width=1000, height=1000)
        par(mai=c(size,size,0.5,0.2))
        plot.new()
        lim = c(-1, N+1)
        plot.window(xlim=lim, ylim=lim)
        rect(sm$i-1,sm$j-1,sm$i,sm$j,col=sm$col, border="lightgray")
        mtext(text=names, side=1, at=1:N-0.5, las=2, line=-2, cex=0.8)
        mtext(text=names, side=2, at=1:N-0.5, las=2, line=-2, cex=0.8)
        title(main="number of shared contigs between clusters")
        dev.off()
    }

    library(gplots)

    cat(sprintf("reading ga table: %s\n", ifn))
    genes = read.delim(ifn)

    s = split(genes$gene, genes$anchor)
    N = length(s)
    m = matrix(1:N^2,N,N)
    for (i in 1:N) {
        for (j in 1:N) {
            m[i,j] = ifelse(i==j,0,length(intersect(s[[i]], s[[j]])))
        }
    }

    cat(sprintf("reading anchor rep table: %s\n", ifn.rep))
    rep = read.delim(ifn.rep)

    otypes = list.files(order.dir)
    for (i in 1:length(otypes)) {
        title = otypes[i]
        ifn.order = paste(order.dir, "/", otypes[i], sep="")
        otable = read.delim(ifn.order)
        anchors = otable$set
        rep.o = match(rep$anchor, anchors)

        o = match(anchors, as.numeric(names(s)))
        mx = m[o,o]

        sm = matrix2smatrix(mx)

        ncols = 256
        colors = c("white", "blue", "red")
        panel = make.color.panel(colors, ncols)
        breaks = c(0, 10, 100)
        sm$col = panel[vals.to.cols(sm$value, breaks, ncols=ncols)]

        system(paste("mkdir -p", fdir))

        p(ofn=paste(fdir, "/", title, "_shared_matrix_species.png", sep=""), names=paste(rep$name[rep.o], anchors), size=3)
        p(ofn=paste(fdir, "/", title, "_shared_matrix_lineage.png", sep=""), names=paste(rep$lineage[rep.o], anchors), size=5)
    }
}
