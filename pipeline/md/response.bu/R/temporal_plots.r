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

get.temporal.table=function(tables, ids, contigs, field="sum", baseline.idx=1)
{
    N = length(tables)

    cat(sprintf("uniting %d profiles\n", length(tables)))
    table = data.frame(contig=contigs)

    for (i in 1:N) {
        p = load.table(tables[i])
        p$value = p[,field]
        df = data.frame(contig=p$contig, value=p$value)
        id = ids[i]
        table[[id]] = p$value[match(contigs, p$contig)]
    }
    table
}

get.anchor.response=function(table, ids, func=mean)
{
    # anchor response
    s = split(table[,ids], table$anchor)
    response = NULL
    for (i in 1:length(s)) {
        anchor = names(s)[i]
        x = s[[i]]
        response = rbind(response, cbind(data.frame(anchor=anchor), t(as.matrix(sapply(x, func)))))
    }
    response
}

compute.contig.response=function(ifn, assembly.dir, ids, ofn)
{
    options(stringsAsFactors=F)
    contigs = load.table(ifn)$contig

    ttables = paste(assembly.dir, "/datasets/", ids, "/coverage/table", sep="")
    xt = get.temporal.table(ttables, ids=ids, contigs=contigs, field="relative.abundance")

    save.table(xt, ofn)
}

compute.response=function(cc.table, assembly.dir, ids, ofn.response, ofn.order)
{
    options(stringsAsFactors=F)

    cat(sprintf("reading cc table: %s\n",cc.table))
    ctable = read.delim(cc.table)

    contigs = unique(ctable$contig)

    ttables = paste(assembly.dir, "/datasets/", ids, "/coverage/table", sep="")
    xt = get.temporal.table(ttables, ids=ids, contigs=contigs, field="sum")
    xt$anchor = ctable$anchor[match(xt$contig, ctable$contig)]
    cr = get.anchor.response(table=xt, ids=ids, func=sum)
    anchors = as.numeric(cr$anchor)

    cat(sprintf("saving anchor response: %s\n", ofn.response))
    write.table(cr, ofn.response, quote=F, col.names=T, row.names=F, sep="\t")

    # cluster anchors
    marginal = colSums(cr[,-1])
    mratio = marginal / marginal[1]
    v = t(t(cr[,-1]) / mratio)
    z = (v[,1]+v[,2]+v[,3])/3
    dd = as.matrix(log10(v / z))
    dd[v<10] = -5
    hc = hclust(dist(dd))
    df = data.frame(set=anchors[hc$order])

    cat(sprintf("saving anchor response order: %s\n", ofn.order))
    write.table(df, ofn.order, quote=F, col.names=T, row.names=F, sep="\t")
}

plot.response=function(ifn.response, order.ifn, taxa.ifn, ids, fdir)
{
    options(stringsAsFactors=F)
    N = length(ids)
    library(gplots)

    system(paste("mkdir -p", fdir))

    cat(sprintf("reading anchor order: %s\n", order.ifn))
    anchors = rev(read.delim(order.ifn)$set)
    M = length(anchors)

    taxa = load.table(taxa.ifn)
    taxa = taxa[match(anchors, taxa$anchor),]

    cat(sprintf("reading anchor response: %s\n", ifn.response))
    cr = read.delim(ifn.response)

    marginal = colSums(cr[,-1])
    mratio = marginal / marginal[1]
    cc = cor(t(t(cr[,-1])/mratio), method="spearman")

    #v = t(t(cr[,-1]) / mratio)
    #z = (v[,1]+v[,2]+v[,3])/3
    #dd = as.matrix(log10(v / z))
    #dd[v<10] = -5

    m = as.matrix(cr[,-1])
    dd = log10((m[,-1]) / (1+m[,1]))


    ncols = 256
    #colors = c("black", "black", "darkgreen", "blue", "white", "red", "orange")
    #panel = make.color.panel(colors, ncols)
    # breaks = c(-5, -3, -3, -1, 0, 1, 2)

    colors = c("black", "blue", "white", "red")
    panel = make.color.panel(colors=colors)
    breaks = c(-3, -1, 0, 1)
    wlegend(fdir=fdir, names=breaks, cols=colors, title="fold_change")

    or = match(anchors, cr$anchor)
    ofn = paste(fdir, "/anchor_time_mat_clust.png", sep="")
    zz = dd[or,]
    labs = anchors

    sm = matrix2smatrix(zz)
    sm$col = ifelse(is.finite(sm$value), panel[vals.to.cols(sm$value, breaks, ncols=ncols)], "black")

    cat(sprintf("generating figure: %s\n", ofn))
    png(ofn, 400, 800)

    par(mai=c(1,2,1,1))
    plot.new()
    xlim = c(0,N+1)
    ylim = c(0,M)
    plot.window(xlim=xlim, ylim=ylim)
    rect(sm$j-1, sm$i-1, sm$j, sm$i, col=sm$col, border=NA)

    # taxa
    rect(N, 1:M-1, N+1, 1:M, border=NA, col=taxa$group.color)
    text(x=N+0.5, y=1:M-0.5, taxa$group.sign)

    abline(v=1:N, lty=3)
    abline(v=c(3,7), lwd=2)
    axis(side=1, labels=ids, at=1:N-0.5, las=2)
    axis(side=2, labels=labs, at=1:M-0.5, las=2)

    dev.off()

    # plot general correlation
    ofn = paste(fdir, "/lib_mat.png", sep="")
    cat(sprintf("generating figure: %s\n", ofn))
    png(ofn, 400, 400)
    par(mai=c(2,2,1,1))
    n = 300
    col = colorpanel(n, "white", "red")
    breaks = seq(min(cc),1,length.out=n+1)
    image(z=cc, zlim=c(-1,1), col=col, breaks=breaks, axes=F, main=round(min(cc),2))
    axis(side=1, labels=ids, at=seq(0,1,length.out=N), las=2)
    axis(side=2, labels=ids, at=seq(0,1,length.out=N), las=2)
    dev.off()
}

plot.temporal.genes=function(cc.table, gene.table, assembly.dir, tdatasets, odir)
{
    options(stringsAsFactors=F)

    cat(sprintf("reading cc table: %s\n",cc.table))
    ctable = read.delim(cc.table)
    anchors = unique(ctable$anchor)

    field.count=function(x, field="contig") {
        tt = table(x[,field])
        result = data.frame(x=names(tt), count=as.vector(tt))
        names(result)[1] = field
        result[order(result$count, decreasing=T),]
    }

    df = field.count(ctable)

    cat(sprintf("reading gene table: %s\n",gene.table))
    gene = read.delim(gene.table)
    gene$count = df$count[match(gene$contig, df$contig)]
    gene = gene[gene$count >= 1,]
    s = split(gene, gene$class)
    classes = names(s)
    classes = setdiff(classes, c("hypothetical protein", "no_hit", "uncharacterized protein", "hit"))
#    classes = "phage"

    ttables = paste(assembly.dir, "/datasets/", tdatasets, "/coverage/table", sep="")
    xt = get.temporal.table(ttables, ids=tdatasets, contigs=unique(ctable$contig), field="median")
    xt$anchor = ctable$secondary_anchor[match(xt$contig, ctable$contig)]
    cr = get.anchor.response(table=xt, ids=tdatasets)

    N = length(tdatasets)

    for (class in classes) {
        x = s[[class]]
        contigs = unique(x$contig)

        # fix spaces in class name
        class.f = gsub(" ", "_", class)

        odir.class = paste(odir, "/", class.f, sep="")
        system(paste("mkdir -p", odir.class))
        cat(sprintf("plotting class %s into dir: %s\n", class, odir.class))

        odir.compare = paste(odir.class, "/compare", sep="")
        system(paste("mkdir -p", odir.compare))

        for(i in 1:length(contigs)) {
            contig = contigs[i]
            anchors = ctable$secondary_anchor[ctable$contig == contig]
            M = length(anchors)

            cvalues = xt[xt$contig == contig,-c(1,14)]

            asum = NULL
            avalues = NULL
            for (j in 1:M) {
                anchor = anchors[j]
                if (!is.element(anchor, cr$anchor))
                    next
                vv = cr[cr$anchor == anchor,-1]
                avalues = rbind(avalues, vv)
                if (j == 1)
                    asum = avalues
                else
                    asum = asum + vv
            }

            xlim = c(1,N+2)
            ylim = range(c(0, asum, avalues, cvalues))
            ofn = paste(odir.class, "/", i, "_", contig, ".png", sep="")
            png(ofn, 600, 600)
            plot(1:N, t(as.vector(cvalues)), type="l", lwd=2, xlim=xlim, ylim=ylim)
            colors = 1:M + 1
            for (j in 1:M)
                lines(1:N, avalues[j,], col=j+1)
            lines(1:N, asum, lty=2, lwd=2)
            legend("topright", fill=colors, legend=anchors)
            dev.off()
        }
    }
}
