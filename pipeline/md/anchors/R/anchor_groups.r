group.cluster=function(smat, anchor.table, min.identity, group.threshold)
{
    smat = smat[smat$set1 != "NONE" & smat$set2 != "NONE",]
    anchors = sort(anchor.table$set)
    N = length(anchors)

    smat$set1.i = match(smat$set1, anchors)
    smat$set2.i = match(smat$set2, anchors)
    m = smatrix2matrix(smat, dim=c(N,N), i.field="set1.i", j.field="set2.i", value.field="identity", default.value=min.identity)
    m[m<min.identity] = min.identity
    rownames(m) = anchor.table$id[match(anchors,anchor.table$set)]
    colnames(m) = anchor.table$id[match(anchors,anchor.table$set)]
    d = 100 - as.dist(m)
    hc = hclust(d, method="average")
    ct = cutree(hc, h=100-group.threshold)
    result = data.frame(anchor.id=names(ct), anchor=anchor.table$set[match(names(ct),anchor.table$id)], group.index=ct)
    result = result[hc$order,]
    indices = unique(result$group.index)
    result$group = match(result$group.index, indices)
    list(hc=hc, ct=ct, table=result)
}

plot.group.dendrogram=function(mat.ifn, group.threshold, order.ifn, min.identity, fdir)
{
    smat = load.table(mat.ifn)
    anchor.table = load.table(order.ifn)

    ll = group.cluster(smat=smat, anchor.table=anchor.table, min.identity=min.identity, group.threshold=group.threshold)
    hc = ll$hc
    ct = ll$ct
    result.anchors = ll$table

    rg = table(result.anchors$group)
    result.groups = data.frame(group=names(rg), count=as.vector(rg))

    group.colors = rainbow(dim(result.groups)[1])
    colLab = function(n) {
        if (is.leaf(n)) {
            a = attributes(n)
            group.index = ct[which(names(ct) == a$label)]
            group = result.anchors$group[match(group.index, result.anchors$group.index)]
            col = group.colors[group]
            attr(n, "nodePar") = c(a$nodePar, list(pch=19, col=col))
            attr(n, "label") = paste(a$label, group)
        }
        n
    }
    dd = dendrapply(as.dendrogram(hc), colLab)

    fig.start(fdir=fdir, paste(fdir, "/dendrogram_", group.threshold, ".pdf", sep=""), type="pdf", width=10, height=6)
    plot(dd, main=paste("threshold=",group.threshold, sep=""))
#    plot(hc, hang=-1, ylab="divergence", las=2)
    abline(h=100-group.threshold, lty=2)
    grid()
    fig.end()
}

group.anchors=function(mat.ifn, group.threshold, order.ifn, min.identity, group.ofn, anchor.ofn)
{
    smat = load.table(mat.ifn)
    anchor.table = load.table(order.ifn)

    ll = group.cluster(smat=smat, anchor.table=anchor.table, min.identity=min.identity, group.threshold=group.threshold)
    hc = ll$hc
    ct = ll$ct

    result.anchors = ll$table
    rg = table(result.anchors$group)
    result.groups = data.frame(group=names(rg), count=as.vector(rg))

    save.table(result.anchors, anchor.ofn)
    save.table(result.groups, group.ofn)
}
