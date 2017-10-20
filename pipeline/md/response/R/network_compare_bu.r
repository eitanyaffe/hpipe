network.compare=function(ifn.ca, ifn.clusters,
    ifn.map1, ifn.map2, ifn.network1, ifn.network2,
    ifn.coverage1, ifn.coverage2,
    min.contacts, ifn.elements, ifn.anchors, ofn)
{
    load.map=function(ifn) {
        map = load.table(ifn)
        map$en = log10(map$observed/map$expected)
        map$en[map$observed < min.contacts] = 0
        map
    }
    cluster.table = load.table(ifn.clusters)
    ca = load.table(ifn.ca)
    map1 = load.map(ifn.map1)
    map2 = load.map(ifn.map2)

    element.table = load.table(ifn.elements)
    anchor.table = load.table(ifn.anchors)

    cov1 = load.table(ifn.coverage1)
    cov2 = load.table(ifn.coverage2)

    anchors = anchor.table$set
    anchor.ids = anchor.table$id
    major.cluster = anchors

    network1 = load.table(ifn.network1)
    network2 = load.table(ifn.network2)

    element.ids = element.table$id
    elements = element.table$cluster

    get.cluster.median.coverage=function(cluster, cov) {
        contigs = cluster.table$contig[match(cluster, cluster.table$cluster)]
        median(cov$abundance.enrichment[match(contigs, cov$contig)])
    }
    result = NULL
    for (i in 1:dim(anchor.table)[1]) {
        anchor.id = anchor.table$id[i]
        anchor = anchor.table$set[i]
        acontigs = ca$contig[ca$anchor == anchor & ca$anchor == ca$contig_anchor]

        anchor.abundance1 = median(cov1$abundance.enrichment[match(acontigs, cov1$contig)])
        anchor.abundance2 = median(cov2$abundance.enrichment[match(acontigs, cov2$contig)])

        c1 = network1$cluster[network1$anchor == anchor]
        c2 = network2$cluster[network2$anchor == anchor]
        cc = sort(unique(c(c1,c2)))
        cc = cc[is.element(cc, element.table$cluster)]
        if (length(cc) == 0)
            next

        cluster.cov1 = NULL
        for (cluster in cc) {
        }

        get.obs=function(map) {
            map = map[map$cluster1 == anchor,]
            ix = match(cc, map$cluster2)
            ifelse(!is.na(ix), map$observed[ix], 0)
        }

        get.score=function(map) {
            map = map[map$cluster1 == anchor,]
            ix = match(cc, map$cluster2)
            ifelse(!is.na(ix), map$en[ix], 0)
        }

        df = data.frame(anchor.id=anchor.id, anchor=anchor,
            element.id=element.table$id[match(cc, element.table$cluster)], cluster=cc,
            observed1=get.obs(map1), observed2=get.obs(map2),
            score1=get.score(map1), score2=get.score(map2),
            anchor.abundance1=anchor.abundance1, anchor.abundance2=anchor.abundance2, element.abundance1=0, element.abundance2=0)
        for (j in 1:length(cc)) {
            cluster = cc[j]
            df$element.abundance1[j] = get.cluster.median.coverage(cluster=cluster, cov=cov1)
            df$element.abundance2[j] = get.cluster.median.coverage(cluster=cluster, cov=cov2)
        }
        df$anchor.abudnance.delta = df$anchor.abundance2 - df$anchor.abundance1
        df$element.abudnance.delta = df$element.abundance2 - df$element.abundance1

        result = rbind(result, df)
    }
    save.table(result, ofn)
}
