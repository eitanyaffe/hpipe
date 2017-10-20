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

compute.anchor.matrix=function(mat.dir, ga.ifn, ofn)
{
    library(gplots)

    cat(sprintf("reading ga: %s\n", ga.ifn))
    ga = read.delim(ga.ifn)
    anchors = sort(unique(ga$anchor))
    s = split(ga$cgene, ga$anchor)
    N = length(s)
    result = NULL
    for (i1 in 1:N) {
        for (i2 in 1:N) {
            if (i1 == i2)
                next
            anchor1 = names(s)[i1]
            anchor2 = names(s)[i2]

            ifn = paste(mat.dir, "/", anchor1, "_", anchor2, sep="")
            if (!file.exists(ifn))
                next
            map = read.delim(ifn)

            cgenes1 = s[[i1]]
            cgenes2 = s[[i2]]

            n1 = length(cgenes1)
            n2 = length(cgenes2)

            # shared cgenes
            shared = unique(intersect(cgenes1, cgenes2))
            n.shared = length(shared)

            cgenes1 = setdiff(cgenes1, shared)
            cgenes2 = setdiff(cgenes2, shared)

            # unmapped cgenes
            n.unmapped1 = length(cgenes1) - length(unique(map$cgene1))
            n.unmapped2 = length(cgenes2) - length(unique(map$cgene2))

            n.strong = sum(map$identity >= 70 & map$coverage >= 0.7)
            n.weak = sum(map$identity < 70 & map$coverage < 0.7)

            # compute average identity score
            map = map[map$identity >= 60 & map$coverage >= 0.7,]

            identity = mean(c(map$identity, rep(100, n.shared)))
            identity.count = length(map$identity)

            df = data.frame(anchor1=anchor1, anchor2=anchor2, identity=identity,
                identity.count=identity.count, shared=n.shared, strong=n.strong, weak=n.weak,
                unmapped1=n.unmapped1, unmapped2=n.unmapped2)
            result = rbind(result, df)
        }
    }

    cat(sprintf("saving result table: %s\n", ofn))
    write.table(result, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}

cluster.anchors=function(mat.ifn, ofn)
{
    smat = read.delim(mat.ifn)
    anchors = sort(unique(smat$anchor1))
    N = length(anchors)
    smat$anchor1.i = match(smat$anchor1, anchors)
    smat$anchor2.i = match(smat$anchor2, anchors)
    m = smatrix2matrix(smat, dim=c(N,N), i.field="anchor1.i", j.field="anchor2.i", value.field="identity", default.value=0)
    d = max(m) - as.dist(m)
    hc = hclust(d)

    df = data.frame(anchor=anchors[hc$order])
    cat(sprintf("saving result: %s\n", ofn))
    write.table(df, ofn, sep="\t", quote=F, col.names=T, row.names=F)
}

plot.anchor.clusters=function(mat.ifn, order.ifn, fdir)
{
    library(gplots)

    smat = read.delim(mat.ifn)
    anchors = read.delim(order.ifn)[,1]
    N = length(anchors)

    smat$anchor1.i = match(smat$anchor1, anchors)
    smat$anchor2.i = match(smat$anchor2, anchors)
    m.identity = smatrix2matrix(smat, dim=c(N,N), i.field="anchor1.i", j.field="anchor2.i", value.field="identity", default.value=0)
    m.shared = smatrix2matrix(smat, dim=c(N,N), i.field="anchor1.i", j.field="anchor2.i", value.field="shared", default.value=0)

    cat(sprintf("plotting results in: %s\n", fdir))
    system(paste("mkdir -p", fdir))

    ncols = 256
    lim = c(0,N)

    plot.mat=function(m,breaks,colors,title) {
        if (length(breaks) != length(colors))
            stop("#breaks != #colors")
        sm = matrix2smatrix(m)
        panel = make.color.panel(colors, ncols)
        sm$col = panel[vals.to.cols(sm$value, breaks, ncols)]
        ofn = paste(fdir, "/matrix_", title, ".png", sep="")
        cat(sprintf("plotting figure: %s\n", ofn))
        png(ofn, 800, 800)
        plot.new()
        plot.window(xlim=lim, ylim=lim)
        title(main=title)
        rect(sm$i-1, sm$j-1, sm$i, sm$j, col=sm$col, border=NA)
        mtext(text=anchors, side=1, at=1:N-0.5, las=2, line=0)
        mtext(text=anchors, side=2, at=1:N-0.5, las=2, line=0)
        dev.off()
    }

    plot.mat(m=m.identity, breaks=c(60,80,95,100), colors=c("white","blue","red","orange"), title="identity")
    plot.mat(m=m.shared, breaks=c(0,10,100), colors=c("white", "red", "orange"), title="shared")
}
