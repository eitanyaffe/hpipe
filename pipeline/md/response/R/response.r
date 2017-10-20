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

compute.contig.response=function(ifn, min.detected, assembly.dir, map.tag, ids, ofn.observed, ofn.expected, ofn.norm, ofn.min.score)
{
    contigs = load.table(ifn)
    ttables = paste(assembly.dir, "/datasets/", ids, "/coverage_", map.tag, "/table", sep="")

    observed.table = get.temporal.table(ttables, ids=ids, contigs=contigs$contig, field="sum")
    save.table(observed.table, ofn.observed)

    min.score = 1

    # uniform expected values
    expected.table = data.frame(contig=contigs$contig)
    total.length = sum(contigs$length)
    for (i in 2:dim(observed.table)[2]) {
        total.reads = sum(observed.table[,i])
        expected = total.reads * (contigs$length / total.length)
        expected.table = cbind(expected.table, expected)
        min.score = min(min.score, min.detected/min(expected))
    }
    names(expected.table) = names(observed.table)
    save.table(expected.table, ofn.expected)

    save.table(data.frame(score=min.score), ofn.min.score)

    # uniform expected values
    scores = (observed.table[,-1]) / (expected.table[,-1])
    scores[scores<min.score] = min.score

    norm.table = data.frame(contig=contigs$contig, scores)
    save.table(norm.table, ofn.norm)
}

compute.anchor.response.global=function(ifn.contigs, ifn.ca, ifn.norm, ofn)
{
    contig.table = load.table(ifn.contigs)
    ca = load.table(ifn.ca)
    norm = load.table(ifn.norm)
    anchors = sort(unique(ca$anchor))

    result = NULL
    for (anchor in anchors) {
        contigs = ca$contig[ca$contig_anchor == anchor]
        anchor.norm = as.matrix(norm[match(contigs, norm$contig),-1])
        anchor.length = contig.table$length[match(contigs, contig.table$contig)]
        anchor.length = anchor.length / sum(anchor.length)
        df = data.frame(anchor=anchor, as.data.frame(round(t(colSums(anchor.norm * anchor.length)),5)))
#        df = data.frame(anchor=anchor, as.data.frame(round(t(colSums(anchor.norm)),5)))
        result = rbind(result, df)
    }
    save.table(result, ofn)
}

response.anchor.order=function(ifn.ca, ifn.median, ifn.detection, class.count, max.height, type, ofn, base.ids)
{
    all.patterns = load.table(ifn.median)
    df = load.table(ifn.ca)
    anchors = sort(unique(df$anchor))

    min.score = log10(load.table(ifn.detection)[1,1])
    min.drop = -2

    xx = as.matrix(all.patterns[match(anchors, all.patterns$cluster),-(1:2)])
    detected = as.matrix(log10(xx) > min.score)
    base = rowSums(xx[,base.ids]) / length(base.ids)
    patterns = ifelse(detected, log10(xx / base), min.drop)
    M = dim(patterns)[2]

    cc = cor(t(patterns))
    cc[is.na(cc)] = -1
    hh = hclust(as.dist(1-cc), method="average")

    ct = switch(type,
        count=cutree(hh, k=class.count),
        height=cutree(hh, h=max.height),
        stop(paste("unknown type")))
    anchors = anchors[hh$order]
    result = data.frame(anchor=anchors, class=ct[match(anchors, names(ct))])
    cat(sprintf("number of classes: %d\n", length(unique(result$class))))
    save.table(result, ofn)
}

mean.cluster=function(ifn, ifn.ca, thresholds, ofn.order, ofn.prefix)
{
    library("fastcluster")

    ca = load.table(ifn.ca)
    table = load.table(ifn)
    table = table[is.element(table$contig, ca$contig),]

    contigs = table$contig

    cat(sprintf("computing pearson over contigs: %d\n", length(contigs)))
    cc = cor(t(table[,-1]))
    cc[is.na(cc) | !is.finite(cc)] = -1

    cat(sprintf("computing hclust...\n"))
    hh = fastcluster::hclust(as.dist(1 - cc), method="average")

    ordered.contigs = contigs[hh$order]
    save.table(data.frame(contig=ordered.contigs), ofn.order)

    for (threshold in thresholds) {
        ofn = paste(ofn.prefix, "_", threshold, sep="")
        cat(sprintf("cutting tree, threshold: %f\n", threshold))
        result = data.frame(contig=contigs, cluster=cutree(hh, h=1-threshold))
        result = add.field.count(result, "cluster")
        result$cluster[result$cluster_count == 1] = -1
        save.table(result[,1:2], ofn)
    }
}

element.genes=function(ifn.clusters, ifn.genes, ifn.uniref, ofn)
{
    table = load.table(ifn.clusters)
    uniref = load.table(ifn.uniref)
    gene = load.table(ifn.genes)
    gene$index = 1:dim(gene)[1]

    gene = lookup.append(table=gene, lookup.table=uniref, lookup.field="gene", value.field="uniref", omit.na=F)
    gene = lookup.append(table=gene, lookup.table=uniref, lookup.field="gene", value.field="identity", omit.na=F)
    gene = lookup.append(table=gene, lookup.table=uniref, lookup.field="gene", value.field="prot_desc", omit.na=F)
    gene = lookup.append(table=gene, lookup.table=uniref, lookup.field="gene", value.field="tax", omit.na=F)

    df = lookup.append(table=gene, lookup.table=table, lookup.field="contig", value.field="cluster")
    df = df[order(df$cluster, df$index),]
    df = df[,c("cluster", "contig", "gene", "uniref", "identity", "prot_desc", "tax")]
    save.table(df, ofn)
}

pattern.response=function(ifn.norm, ifn.obs, ifn.exp, ifn.clusters,
    ofn.obs, ofn.exp, ofn.mean, ofn.sd, ofn.top95, ofn.top75, ofn.median, ofn.bottom25, ofn.bottom05, ofn.bottom0, ofn.top100)
{
    norm = load.table(ifn.norm)
    obs = load.table(ifn.obs)
    exp = load.table(ifn.exp)
    cluster.table = load.table(ifn.clusters)
    cluster.table = cluster.table[cluster.table$cluster > 0,]
    clusters = sort(unique(cluster.table$cluster))
    N = dim(norm)[2] - 1
    names = names(norm)[-1]

    result.mean = NULL
    result.top100 = NULL
    result.top95 = NULL
    result.top75 = NULL
    result.median = NULL
    result.bottom0 = NULL
    result.bottom05 = NULL
    result.bottom25 = NULL
    result.sd = NULL

    result.obs = NULL
    result.exp = NULL

    for (cluster in clusters) {
        contigs = cluster.table$contig[cluster.table$cluster == cluster]
        size = length(contigs)

        # normalized matrix
        rx = as.matrix(norm[is.element(norm$contig,contigs),-1])
        result.mean = rbind(result.mean, data.frame(cluster=cluster, size=size, t(apply(rx, 2, mean))))
        result.top100 = rbind(result.top100, data.frame(cluster=cluster, size=size, t(apply(rx, 2, function(x) quantile(x, 1)))))
        result.top95 = rbind(result.top95, data.frame(cluster=cluster, size=size, t(apply(rx, 2, function(x) quantile(x, 0.95)))))
        result.top75 = rbind(result.top75, data.frame(cluster=cluster, size=size, t(apply(rx, 2, function(x) quantile(x, 0.75)))))
        result.median = rbind(result.median, data.frame(cluster=cluster, size=size, t(apply(rx, 2, function(x) quantile(x, 0.5)))))
        result.bottom0 = rbind(result.bottom0, data.frame(cluster=cluster, size=size, t(apply(rx, 2, function(x) quantile(x, 0)))))
        result.bottom05 = rbind(result.bottom05, data.frame(cluster=cluster, size=size, t(apply(rx, 2, function(x) quantile(x, 0.05)))))
        result.bottom25 = rbind(result.bottom25, data.frame(cluster=cluster, size=size, t(apply(rx, 2, function(x) quantile(x, 0.25)))))
        ssd = if (size > 1) t(apply(rx, 2, sd)) else t(rep(0, N))
        colnames(ssd) = names
        result.sd = rbind(result.sd, data.frame(cluster=cluster, size=size, ssd))

        # obs
        ox = as.matrix(obs[is.element(obs$contig,contigs),-1])
        result.obs = rbind(result.obs, data.frame(cluster=cluster, t(colSums(ox))))

        # exp
        ex = as.matrix(exp[is.element(exp$contig,contigs),-1])
        result.exp = rbind(result.exp, data.frame(cluster=cluster, t(colSums(ex))))
    }

    save.table(result.obs, ofn.obs)
    save.table(result.exp, ofn.exp)
    save.table(result.mean, ofn.mean)
    save.table(result.top100, ofn.top100)
    save.table(result.top95, ofn.top95)
    save.table(result.top75, ofn.top75)
    save.table(result.median, ofn.median)
    save.table(result.bottom05, ofn.bottom05)
    save.table(result.bottom25, ofn.bottom25)
    save.table(result.bottom0, ofn.bottom0)
    save.table(result.sd, ofn.sd)
}

anchor.major=function(ifn.contigs, ifn.clusters, ifn.cluster.map, ifn.ca, ofn.major, ofn.clusters.final)
{
    contigs = load.table(ifn.contigs)
    cluster.table = load.table(ifn.clusters)

    ca =  load.table(ifn.ca)
    anchors = sort(unique(ca$anchor))

    cluster.table$cluster = cluster.table$cluster.consistent
    cluster.table$length = contigs$length[match(cluster.table$contig, contigs$contig)]
    ca$length = contigs$length[match(ca$contig, contigs$contig)]

    major.df = NULL
    major.result = NULL
    for (anchor in anchors) {
        ca.anchor = ca[ca$anchor == anchor,]
        anchor.size = sum(ca.anchor$length)

        # add cluster
        ca.anchor$cluster = cluster.table$cluster[match(ca.anchor$contig, cluster.table$contig)]

        # limit to valid clusters
        ca.anchor = ca.anchor[ca.anchor$cluster > 0,]

        # limit to anchor contigs
        ca.anchor = ca.anchor[ca.anchor$contig_anchor != 0,]

        ss = sapply(split(ca.anchor$length, ca.anchor$cluster), sum)
        cluster = as.numeric(names(which.max(ss)))
        match.size = max(ss)

        full.size = sum(ca.anchor$length[ca.anchor$contig_anchor == anchor])
        df = data.frame(anchor=anchor, cluster=cluster, anchor.f=match.size/anchor.size, anchor.full.f=full.size/anchor.size)
        major.df = rbind(major.df, df)

        anchor.contigs = ca.anchor$contig[ca.anchor$cluster == cluster]
        major.result = rbind(major.result, data.frame(contig=anchor.contigs, cluster=cluster))
    }
    save.table(major.df, ofn.major)

    result = cluster.table[!is.element(cluster.table$cluster, major.df$cluster),c("contig", "cluster")]
    result = rbind(result, major.result)
    save.table(result, ofn.clusters.final)
}

anchor.elements=function(ifn.matrix, ifn.means, ofn)
{
    ea = load.table(ifn.matrix)
    means = load.table(ifn.means)
    anchors = sort(unique(ea$anchor))

    result = NULL
    for (anchor in anchors) {
        e.major = means[means$cluster == anchor,-(1:2)]
        if (!any(ea$anchor == anchor)) {
            next
        }
        clusters = ea$cluster[ea$anchor == anchor]
        e.minors = means[is.element(means$cluster, clusters),-(1:2)]
        cc = cor(t(e.major), t(e.minors))
        cc[is.na(cc) | !is.finite(cc)] = -1
        df = data.frame(anchor=anchor, cluster=clusters, round(t(cc),4))
        names(df)[3] = "pearson"
        df$major = df$cluster == anchor
        df = df[order(df$pearson),]
        result = rbind(result, df)
    }
    save.table(result, ofn)
}
