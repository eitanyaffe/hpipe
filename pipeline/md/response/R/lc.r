collect.reads=function(ifn.elements, ifn.network, ifn.clusters, ifn.anchor.contigs, ifn.anchor.ids, paired.dir, odir)
{
    element.table = load.table(ifn.elements)
    network = load.table(ifn.network)
    cluster.table = load.table(ifn.clusters)
    anchor.table.contigs = load.table(ifn.anchor.contigs)
    anchor.table.ids = load.table(ifn.anchor.ids)

    result = NULL
    clusters = element.table$cluster
    for (cluster in clusters) {
        anchors = network$anchor[network$cluster == cluster]
        for (anchor in anchors) {
            element.id = element.table$id[match(cluster, element.table$cluster)]
            anchor.id = anchor.table.ids$id[match(anchor, anchor.table.ids$set)]
            anchor.contigs =
        }
    }
}

