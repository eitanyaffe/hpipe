smatrix2matrix=function(smatrix, dim, i.field="i", j.field="j", value.field="value", default.value=0)
{
  indices = smatrix[,i.field] + (smatrix[,j.field]-1) * dim[1]
  v = rep(default.value, dim[1]*dim[2])
  v[indices] = smatrix[,value.field]
  matrix(v, dim[1], dim[2])
}

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

compute.matrix=function(ifn, ofn)
{
    tree = read.delim(ifn)
    anchors = sort(unique(tree$anchor))
    N = length(anchors)

    fields = c("tax_id", "weight")
    cat(sprintf("going over pairwise anchors to count number of hits which meet criteria...\n"))
    result = NULL
    for (i1 in 1:N) {
    for (i2 in 1:N) {
        anchor1 = anchors[i1]
        anchor2 = anchors[i2]
        if (i1 == i2) {
            dist = 0
        } else {
            tree1 = tree[tree$anchor == anchor1,]
            tree2 = tree[tree$anchor == anchor2,]
            m = merge(tree1[,fields], tree2[,fields], by="tax_id", all=T)
            m[is.na(m)] = 0
            x1 = m$weight.x
            x2 = m$weight.y
            x1 = x1/sum(x1)
            x2 = x2/sum(x2)
            dist = sum(abs(x1-x2))
        }
        result = rbind(result, data.frame(anchor1=anchor1, anchor2=anchor2, dist=dist))
    } }

    cat(sprintf("saving distance matrix: %s\n", ofn))
    write.table(result, ofn, sep="\t", quote=F, col.names=T, row.names=F)
}

cluster.anchors=function(ifn, ofn)
{
    table = read.delim(ifn)
    anchors = sort(unique(table$anchor1))
    N = length(anchors)
    dim = c(N,N)
    table$anchor1 = match(table$anchor1, anchors)
    table$anchor2 = match(table$anchor2, anchors)

    m = smatrix2matrix(table, dim, i.field="anchor1", j.field="anchor2", value.field="dist", default.value=0)

    lim = c(0,N)
    method = "average"

    d = as.dist(m)
    hc = hclust(d, method=method)
    o = hc$order

    cat(sprintf("saving sorted anchors into table: %s\n", ofn))
    otable = data.frame(set=anchors[o])
    write.table(otable, ofn, sep="\t", quote=F, col.names=T, row.names=F)
}

plot.matrix=function(ifn, order.dir, fdir)
{
    library(gplots)
    library(smacof)
    table = read.delim(ifn)
    anchors = sort(unique(table$anchor1))
    N = length(anchors)
    dim = c(N,N)
    table$anchor1 = match(table$anchor1, anchors)
    table$anchor2 = match(table$anchor2, anchors)

    m = smatrix2matrix(table, dim, i.field="anchor1", j.field="anchor2", value.field="dist", default.value=0)

    # plot gene count
    fig.fn = paste(fdir, "/dendrogram.png", sep="")
    png(fig.fn, 1000, 400)
    method = "average"
    d = as.dist(m)
    hc = hclust(d, method=method)
    plot(hc, labels=anchors, hang=-1)
    dev.off()

    lim = c(0,N)

    ncols = 256
    colors = c("orange", "red", "blue", "white")
    panel = make.color.panel(colors, ncols)
    breaks = c(0, 0.2, 0.5, 1.1)

    otypes = list.files(order.dir)
    for (i in 1:length(otypes)) {
        title = otypes[i]
        ifn.order = paste(order.dir, "/", otypes[i], sep="")
        cat(sprintf("reading order table: %s\n", ifn.order))
        otable = read.delim(ifn.order)
        o = match(otable$set, anchors)

        mo = m[o,o]
        sm = matrix2smatrix(mo)
        sm$col = panel[vals.to.cols(sm$value, breaks, ncols)]

        system(paste("mkdir -p", fdir))
        fig.fn = paste(fdir, "/", title, "_matrix.png", sep="")
        cat(sprintf("plotting: %s\n", fig.fn))
        png(fig.fn, 1000, 1000)
        plot.new()
        plot.window(xlim=lim, ylim=lim)
        rect(sm$i-1, sm$j-1, sm$i, sm$j, col=sm$col, border=NA)
        mtext(text=anchors[o], side=1, at=1:N-0.5, las=2, line=0)
        mtext(text=anchors[o], side=2, at=1:N-0.5, las=2, line=0)
        dev.off()
    }
}
