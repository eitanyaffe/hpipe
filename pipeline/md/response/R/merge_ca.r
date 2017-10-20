merge.ca=function(ifn.anchors, ifn.contigs, ifn.clusters, ifn.ca,
    min.contacts, min.enrichment, separate.min.contacts, separate.max.enrichment, ofn)
{
    contig.cluster.table = load.table(ifn.clusters)
    anchor.table = load.table(ifn.anchors)
    contig.table = load.table(ifn.contigs)
    ca = load.table(ifn.ca)

    # anchors depend only on ca
    anchors = sort(unique(ca$anchor))

    contig.cluster.table = contig.cluster.table[!contig.cluster.table$anchor,]
    contig.cluster.table$contig.length = contig.table$length[match(contig.cluster.table$contig, contig.table$contig)]

    # contigs are intersection of contigs in clusters and contigs in ca
    contigs = intersect(unique(ca$contig), contig.cluster.table$contig)

    # we trim cluster table according to relevant contigs
    contig.cluster.table = contig.cluster.table[is.element(contig.cluster.table$contig, contigs),]

    # add cluster size
    xs = sapply(split(contig.cluster.table$contig.length, contig.cluster.table$cluster), sum)
    cluster.length = data.frame(cluster=names(xs), length=xs)
    clusters = sort(unique(contig.cluster.table$cluster))

    cat(sprintf("number of anchors: %d\n", length(anchors)))
    cat(sprintf("number of clusters: %d\n", length(clusters)))
    cat(sprintf("number of contigs: %d\n", length(contigs)))

    # add cluster to ca and restrict ca
    ca$cluster = contig.cluster.table$cluster[match(ca$contig, contig.cluster.table$contig)]
    ca = ca[!is.na(ca$cluster) & is.element(ca$cluster,clusters),]

    append.key = function(x) { paste(x$anchor, x$cluster, sep="_") }
    get.score = function(obs,exp) { ifelse(obs>0,log10(obs/exp),0) }
    get.high.score = function(obs,exp) { get.score(pmax(1,obs+sqrt(obs)),exp) }
    get.low.score = function(obs,exp) { get.score(pmax(0,obs-sqrt(obs)),exp) }

    # add info to ca table
    ca$obs = ca$contig_total_count
    ca$exp = ca$contig_expected
    ca$score = get.score(ca$obs, ca$exp)
    ca$high.score = get.high.score(ca$obs, ca$exp)
    ca$low.score = get.low.score(ca$obs, ca$exp)

    ca$contig.length = contig.table$length[match(ca$contig, contig.table$contig)]
    ca$cluster.length = cluster.length$length[match(ca$cluster, cluster.length$cluster)]
    ca$weight = ca$contig.length / ca$cluster.length

    dummy = sapply(split(ca, ca$anchor), function(x) {
        if (!setequal(x$contig, contigs))
            stop(sprintf("anchor %s missing some contigs in ca matrix", x$anchor[1]))
    } )

    ca$key = append.key(ca)
    ss = split(ca, ca$key)
    append.f=function(df, field, f) {
        cat(sprintf("computing field: %s\n", field))
        rs = sapply(ss, f)
        if (!setequal(names(rs),df$key))
            stop("missing keys")
        df[,field] = rs[match(df$key,names(rs))]
        df
    }

    gg = expand.grid(anchors, clusters)
    df = data.frame(anchor=gg[,1], cluster=gg[,2])
    df$key = append.key(df)

    df = append.f(df=df, field="sum.score", f=function(x) { get.score(sum(x$obs),sum(x$exp)) })
    df = append.f(df=df, field="n.connected.contigs", f=function(x) { sum(x$obs>0) } )
    df = append.f(df=df, field="n.contigs", f=function(x) { length(x$obs) } )
    df = append.f(df=df, field="sum.observed", f=function(x) { sum(x$obs) } )
    df = append.f(df=df, field="sum.expected", f=function(x) { sum(x$exp) } )

    # mean score (weighted by contig length)
    df = append.f(df=df, field="mean.score", f=function(x) { sum(x$weight * x$score) })
    df = append.f(df=df, field="high.mean.score", f=function(x) { sum(x$weight * x$high.score) })
    df = append.f(df=df, field="low.mean.score", f=function(x) { sum(x$weight * x$low.score) })

    # the weighted average is used since it is the most conservative
    df$score = df$mean.score
    df$high.score = df$high.mean.score
    df$low.score = df$low.mean.score
    df$observed = df$sum.observed
    df$expected = df$sum.expected

    df$type = ifelse(df$score >= min.enrichment & df$observed >= min.contacts, "connected",
        ifelse(df$score <= separate.max.enrichment & df$expected*10^separate.max.enrichment >= separate.min.contacts, "separated", "unknown"))

    save.table(df, ofn)
}
