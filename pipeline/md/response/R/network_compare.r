network.compare=function(
    ifn.ca, ifn.clusters, ifn.anchors,
    ifn.map1, ifn.map2,
    ifn.network1, ifn.network2,
    ifn.coverage1, ifn.coverage2,
    min.anchor.abundance, min.element.abundance,
    min.contacts, min.enrichment, separate.min.contacts, separate.max.enrichment,
    ofn.anchors, ofn.clusters, ofn.map)
{
    ca = load.table(ifn.ca)
    map1 = load.table(ifn.map1)
    map2 = load.table(ifn.map2)
    anchor.table = load.table(ifn.anchors)
    cluster.table = load.table(ifn.clusters)
    clusters = sort(unique(cluster.table$cluster))
    cov1 = load.table(ifn.coverage1)
    cov2 = load.table(ifn.coverage2)

    # anchors
    result.anchors = NULL
    for (i in 1:dim(anchor.table)[1]) {
        anchor.id = anchor.table$id[i]
        anchor = anchor.table$set[i]
        acontigs = ca$contig[ca$anchor == anchor & ca$anchor == ca$contig_anchor]

        abundance1 = median(cov1$abundance.enrichment[match(acontigs, cov1$contig)])
        abundance2 = median(cov2$abundance.enrichment[match(acontigs, cov2$contig)])

        lost = abundance2 < min.anchor.abundance
        df = data.frame(anchor.id=anchor.id, anchor=anchor, abundance1=abundance1, abundance2=abundance2, lost=lost)
        result.anchors = rbind(result.anchors, df)
    }
    anchors = result.anchors$anchor[!result.anchors$lost]

    # clusters
    get.cluster.abundance = function(clusters, cov) {
        cov$cluster = cluster.table$cluster[match(cov$contig, cluster.table$contig)]
        ss = sapply(split(cov$abundance.enrichment, cov$cluster), median)
        df = data.frame(cluster=names(ss), abundance=ss)
        ix = match(clusters, df$cluster)
        if (any(is.na(ix)))
            stop("not all clusters found")
        df$abundance[ix]
    }
    result.clusters = data.frame(cluster=clusters)
    result.clusters$abundance1 = get.cluster.abundance(clusters, cov1)
    result.clusters$abundance2 = get.cluster.abundance(clusters, cov2)
    result.clusters$delta = result.clusters$abundance2 - result.clusters$abundance1
    result.clusters$lost = result.clusters$abundance2 < min.element.abundance
    clusters = result.clusters$cluster[!result.clusters$lost]

    # map
    mm = merge(map1, map2, by=c("anchor", "cluster"))
    map = data.frame(
        anchor=mm$anchor, cluster=mm$cluster,
        observed1=mm$observed.x, expected1=mm$expected.x, score1=mm$score.x, sd.score1=abs(mm$high.score.x-mm$low.score.x)/2, type1=mm$type.x,
        observed2=mm$observed.y, expected2=mm$expected.y, raw.score2=mm$score.y, sd.score2=abs(mm$high.score.y-mm$low.score.y)/2, type2=mm$type.y)
    map = map[is.element(map$anchor, anchors) & is.element(map$cluster, clusters),]
    map$mid = paste(map$anchor, map$cluster, sep="_")
    map$cluster.abundance.factor = 10^result.clusters$delta[match(map$cluster, result.clusters$cluster)]

    ## add.score=function(map, o.field, e.field, a.field, score.field, score.sd.field) {
    ##     exp = map[,e.field]
    ##     obs = ifelse(map[,o.field]>0, map[,o.field], exp)
    ##     map[,score.field] = log10(obs/exp)
    ##     high.obs = pmax(1,obs+sqrt(obs))
    ##     map[,score.sd.field] = log10(high.obs/exp) - map[,score.field]
    ##     map
    ## }
    ## map = add.score(map=map, o.field="observed1", e.field="expected1", score.field="score1", score.sd.field="sd.score1")
    ## map = add.score(map=map, o.field="observed2", e.field="expected2", score.field="raw.score2", score.sd.field="sd.score2")

    s = split(map, map$anchor)
    result.anchors$score.factor = 0
    for (i in 1:length(s)) {
        anchor = names(s)[i]
        anchor.id = anchor.table$id[match(anchor, anchor.table$set)]
        imap = s[[i]]
        imap = imap[imap$score1 > 0 & imap$raw.score2 > 0,]
        fit = lm(imap$raw.score2 ~ 0 + imap$score1)
        result.anchors$score.factor[match(anchor, result.anchors$anchor)] = fit$coefficient
    }

    # delta values
    map$score2 = map$raw.score2 / result.anchors$score.factor[match(map$anchor, result.anchors$anchor)]
    map$delta.score = map$score2 - map$score1
    # map$delta.z = map$delta.score / pmax(map$sd.score1, map$sd.score2)

    # recompute type2 since scores were corrected
    map$type2 = ifelse(map$score2 >= min.enrichment & map$observed2 >= min.contacts, "connected",
        ifelse(map$score2 <= separate.max.enrichment & map$expected2*10^separate.max.enrichment >= separate.min.contacts, "separated", "unknown"))

    map$type =
        ifelse(map$type1 == "connected" & map$type2 == "connected", "stable",
               ifelse(map$type1 == "connected" & map$type2 == "separated", "lost",
                      ifelse(map$type2 == "connected" & map$type1 == "separated", "gained", "unknown")))
#    map$type =
#        ifelse(map$delta.score>min.delta.score, ifelse(map$connected1, "increase", "gained"),
#               ifelse(-map$delta.score>min.delta.score, ifelse(map$connected2, "decrease", "lost"),
#                      ifelse(map$connected1 & map$connected2, "stable", "unknown")))


    save.table(result.anchors, ofn.anchors)
    save.table(result.clusters, ofn.clusters)
    save.table(map, ofn.map)
}

network.compare.select=function(ifn.anchors, ifn.elements, ifn.map, ofn.map)
{
    anchor.table = load.table(ifn.anchors)
    element.table = load.table(ifn.elements)
    map = load.table(ifn.map)
    map = map[is.element(map$cluster, element.table$cluster),]

    map$element.id = element.table$id[match(map$cluster, element.table$cluster)]
    map$anchor.id = anchor.table$id[match(map$anchor, anchor.table$set)]

    map$connected1 = map$type1 == "connected"
    map$connected2 = map$type2 == "connected"
    # sig = (abs(map$delta.score)>min.delta.score) & (map$connected1 | map$connected2)
    map = map[map$connected1|map$connected2,]

    save.table(map, ofn.map)
}
