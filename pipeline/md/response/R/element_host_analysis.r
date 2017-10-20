get.name.of.level=function(ids, taxa, level) {
    result = rep("other", length(ids))
    for (i in 1:length(ids)) {
        id = ids[i]
        tid = id
        while (tid != -1) {
            ix = match(tid, taxa$tax_id)
            if (taxa$level[ix] == level) {
                result[i] = taxa$name[ix]
                break
            }
            tid = taxa$parent_id[ix]
        }
    }
    result
}

get.lca=function(ids, taxa) {
    all.ids = NULL
    for (id in ids) {
        tid = id
        while (tid != -1) {
            all.ids = c(all.ids, tid)
            tid = taxa$parent_id[match(tid, taxa$tax_id)]
        }
        all.ids = c(all.ids, tid)
    }
    df = data.frame(id=unique(all.ids), count=0)
    for (id in ids) {
        tid = id
        while (tid != -1) {
            df$count[df$id == tid] = df$count[df$id == tid] + 1
            tid = taxa$parent_id[match(tid, taxa$tax_id)]
        }
        df$count[df$id == tid] = df$count[df$id == tid] + 1
    }
    tid = ids[1]
    while (df$count[df$id == tid] < length(ids))
        tid = taxa$parent_id[match(tid, taxa$tax_id)]
    return (tid)
}

element.to.level=function(ifn.elements, ifn.network, ifn.taxa, ifn.reps, level, ofn.elements, ofn.matrix, ofn.summary)
{
    elements = load.table(ifn.elements)
    taxa = load.table(ifn.taxa)
    taxa$level[taxa$level == "no rank"] = "strain"
    network = load.table(ifn.network)
    reps = load.table(ifn.reps)
    reps$level.name = get.name.of.level(ids=reps$tax_id, taxa=taxa, level=level)
    network$level.name = reps$level.name[match(network$anchor,reps$anchor)]
    network$index = elements$index[match(network$cluster, elements$cluster)]
    s = split(network$level.name, network$index)

    # summary table
    tt = table(reps$level.name)
    tt = tt[order(tt, decreasing=T)]
    tt = c(tt[names(tt) != "other"], tt[names(tt) == "other"])
    save.table(data.frame(level=names(tt), count=as.vector(tt)), ofn.summary)

    # all.names = c(setdiff(sort(unique(network$level.name)), "other"), "other")
    all.names = names(tt)
    N = length(all.names)

    # per element
    result.elements = NULL
    for (i in 1:length(s)) {
        tt = table(s[[i]])
        vv = rep(0, N)
        vv[match(names(tt), all.names)] = tt
        df = data.frame(index=names(s)[i], t(vv))
        names(df)[-1] = all.names
        result.elements = rbind(result.elements, df)
    }
    save.table(result.elements, ofn.elements)

    # summary matrix
    result.mat = matrix(rep(0, N^2), N, N)
    rownames(result.mat) = all.names
    colnames(result.mat) = all.names
    for (i in 1:length(s)) {
        ids = match(unique(s[[i]]), all.names)
        if (length(ids) < 2)
            next
        gg = expand.grid(ids, ids)
        gg = gg[gg[,1] != gg[,2],]
        gg = unique(gg)
        if (dim(gg)[1] == 0)
            next
        for (j in 1:dim(gg)[1])
            result.mat[gg[j,1], gg[j,2]] = result.mat[gg[j,1], gg[j,2]] + 1
    }
    save.table(result.mat, ofn.matrix)
}

plot.element.host.lca=function(ifn.network, ifn.taxa, ifn.reps, identity.threshold=99, frac.threshold=60, fdir)
{
    taxa = load.table(ifn.taxa)

    taxa$level[taxa$level == "no rank"] = "strain"
    network = load.table(ifn.network)
    reps = load.table(ifn.reps)
    network$tax_id = reps$tax_id[match(network$anchor,reps$anchor)]

    reps = reps[reps$identity > identity.threshold & reps$frac > frac.threshold,]
    reps = reps[is.element(reps$level, c("genus", "species", "strain")),]
    network = network[is.element(network$anchor, reps$anchor),]

    s = split(network$tax_id, network$cluster)
    result = NULL
    for (i in 1:length(s)) {
        ids = s[[i]]
        lca.id = get.lca(ids=ids, taxa=taxa)
        ix = match(lca.id, taxa$tax_id)
        result = rbind(result, data.frame(cluster=names(s)[i], name=taxa$name[ix], level=taxa$level[ix]))
    }

    tt = table(result$level)
    levels = c("phylum", "order", "family", "genus", "species", "strain")
    ix = match(levels, names(tt))
    tt = tt[ix]

    fig.start(fdir=fdir, ofn=paste(fdir, "/levels.png", sep=""), width=200, height=300)
    par(mai=c(1,0.5,0.5,0.5))
    barplot(tt, names.arg=names(tt), las=2, main=sum(tt), col="darkblue", border=NA)
    fig.end()

    tt = table(result$name)
    tt = tt[order(tt)]
    fig.start(fdir=fdir, ofn=paste(fdir, "/names.png", sep=""), width=600)
    par(mai=c(0.5,4,0.5,0.5))
    barplot(tt, names.arg=names(tt), las=2, main=sum(tt), horiz=T, col="darkblue", border=NA)
    fig.end()

}
