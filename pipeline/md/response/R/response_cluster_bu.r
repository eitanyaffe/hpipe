local.cluster.contigs=function(hh, anchors, contigs, cmap, threshold)
{
    N = length(contigs)
    contig.sets = list()
    get.contigs=function(index)  {
        if (index < 0)
            return(contigs[-index])
        else
            return(contig.sets[[as.character(index)]])
    }
    get.anchors=function(contig.set) {
        anchors[colSums(cmap[match(contig.set, contigs),]) > 0]
    }
    all.in.anchor=function(contig.set, anchor) {
        all(cmap[match(contig.set, contigs),match(anchor, anchors)] > 0)
    }
    on.anchor=function(anchor) {
        contigs[which(cmap[,match(anchor, anchors)] > 0)]
    }
    # compute contigs per node
    for (i in 1:(N-1))
        contig.sets[[as.character(i)]] = unique(c(get.contigs(hh$merge[i,1]), get.contigs(hh$merge[i,2])))

    # compute number of unique anchors per node
    cat(sprintf("computing single anchor support...\n"))
    single.anchor = NULL
    for (i in 1:(N-1)) {
        icontigs = contig.sets[[as.character(i)]]
        panchors = get.anchors(icontigs)
        value = T
        for (anchor in panchors) {
            if (!all.in.anchor(contig.set=icontigs, anchor=anchor)) {
                value = F
                break
            }
        }
        single.anchor = rbind(single.anchor, data.frame(index=i, value=value))
    }

    cat(sprintf("assigning clusters...\n"))
    flagged.nodes = NULL
    cluster.table = NULL
    cindex = 1
    for (i in (N-1):1) {
        flagged = is.element(i, flagged.nodes)
        cluster.root = hh$height[i] < (1-threshold) && single.anchor$value[i] && !flagged
        if (cluster.root || flagged) {
            i1 = hh$merge[i,1]
            i2 = hh$merge[i,2]
            if (i1 > 0) flagged.nodes = c(flagged.nodes, i1)
            if (i2 > 0) flagged.nodes = c(flagged.nodes, i2)
        }
        if (cluster.root) {
            cluster.table = rbind(cluster.table, data.frame(contig=get.contigs(i), cluster=cindex))
            cindex = cindex + 1
        }
    }
    singles = setdiff(contigs, cluster.table$contig)
    if (length(singles) > 0) {
        nx = seq_along(singles) + dim(cluster.table)[1]
        cluster.table = rbind(cluster.table, data.frame(contig=singles, cluster=nx))
        cluster.table = add.field.count(cluster.table, "cluster")
    }

    cluster.table = cluster.table[match(contigs[hh$order], cluster.table$contig),]
    cluster.table$cluster = match(cluster.table$cluster, unique(cluster.table$cluster))
    cluster.table
}

cluster.contigs=function(ifn, ifn.ca, threshold, ofn.table, ofn.order)
{
    library("fastcluster")

    ca = load.table(ifn.ca)
    ca = ca[,c("contig", "anchor")]

    table = load.table(ifn)
    rownames(table) = table$contig
    table = table[is.element(table$contig, ca$contig),]

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

    ordered.contigs = contigs[hh$order]
    save.table(data.frame(contig=ordered.contigs), ofn.order)

    ca$value = 1
    ca$contig.i = match(ca$contig, contigs)
    ca$anchor.i = match(ca$anchor, anchors)
    cmap = smatrix2matrix(smatrix=ca, dim=c(N,M), i.field="contig.i", j.field="anchor.i", value.field="value", default.value=0)
    cluster.table = local.cluster.contigs(hh=hh, anchors=anchors, contigs=contigs, cmap=cmap, threshold=threshold)
    save.table(cluster.table, ofn.table)
}
