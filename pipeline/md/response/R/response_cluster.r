# we compute two cluster types on the same hclust tree:
#  single:     contigs of each cluster have a single anchor (required for the anchor major elements)
#  consistent: contigs of each cluster have a common anchor (for all other elements)
cluster.contigs.style=function(hh, anchors, contigs, cmap, threshold, single.anchor)
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
    cat(sprintf("computing node validity...\n"))
    inner.node = NULL
    for (i in 1:(N-1)) {
        icontigs = contig.sets[[as.character(i)]]
        panchors = get.anchors(icontigs)

        if (single.anchor) {
            # single anchor for all contigs
            value = (length(panchors) == 1)
        } else {
            # there exists a single anchor common to all contigs
            value = F
            for (anchor in panchors) {
                if (all.in.anchor(contig.set=icontigs, anchor=anchor)) {
                    value = T
                    break
                }
            }
        }
        inner.node = rbind(inner.node, data.frame(index=i, value=value))
    }

    cat(sprintf("assigning clusters...\n"))
    flagged.nodes = NULL
    cluster.table = NULL
    cindex = 1
    for (i in (N-1):1) {
        flagged = is.element(i, flagged.nodes)
        cluster.root = hh$height[i] < (1-threshold) && inner.node$value[i] && !flagged
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

cluster.contigs=function(ifn, ifn.ca, threshold, ofn.table, ofn.majors)
{
    library("fastcluster")

    table = load.table(ifn)
    rownames(table) = table$contig

    ca = load.table(ifn.ca)
    anchors = sort(unique(ca$anchor))
    M = length(anchors)

    ca.intra = ca[ca$contig_anchor == ca$anchor,c("contig", "anchor")]
    ca.inter = ca[ca$contig_anchor != ca$anchor,c("contig", "anchor")]

    table = table[is.element(table$contig, ca.inter$contig),]

    contigs = table$contig
    N = length(contigs)

    cat(sprintf("computing pearson over contigs: %d\n", length(contigs)))
    cc = cor(t(table[,-1]))
    cc[is.na(cc) | !is.finite(cc)] = -1

    cat(sprintf("computing hclust...\n"))
    dd = as.dist(1 - cc)
    hh = fastcluster::hclust(dd, method="average")

    ca.inter$value = 1
    ca.inter$contig.i = match(ca.inter$contig, contigs)
    ca.inter$anchor.i = match(ca.inter$anchor, anchors)
    cmap = smatrix2matrix(smatrix=ca.inter, dim=c(N,M), i.field="contig.i", j.field="anchor.i", value.field="value", default.value=0)

    # cluster non-anchor contigs
    result = cluster.contigs.style(hh=hh, anchors=anchors, contigs=contigs, cmap=cmap, threshold=threshold, single.anchor=F)
    result$anchor = F
    result$cluster = result$cluster + M
    result = result[order(result$cluster),]

    # add clusters for anchor
    df = field.count(ca.intra, "anchor")
    aresult = data.frame(contig=ca.intra$contig, cluster=ca.intra$anchor, cluster_count=df$count[match(ca.intra$anchor, df$anchor)], anchor=T)
    aresult = aresult[order(aresult$cluster),]

    save.table(rbind(aresult, result), ofn.table)

    majors = data.frame(anchor=anchors, cluster=anchors)
    save.table(majors, ofn.majors)
}
