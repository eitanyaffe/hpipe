element.support=function(ifn.anchors, ifn.clusters, ifn.elements, ifn.network, ifn.ca, ifn.ca.matrix, ofn)
{
    cluster.table = load.table(ifn.clusters)
    ca = load.table(ifn.ca)
    ca.matrix = load.table(ifn.ca.matrix)
    network = load.table(ifn.network)

    clusters = sort(unique(network$cluster))

    result = NULL
    for (cluster in clusters) {
        anchors = network$anchor[network$cluster == cluster]
        contigs = cluster.table$contig[cluster.table$cluster == cluster]
        for (anchor in anchors) {
            anchor.contigs = ca$contig[ca$anchor == anchor & ca$contig_anchor != 0]
            cmat = ca.matrix[ca.matrix$anchor == anchor & is.element(cmat$contig, contigs),]

            # not mentioned at all
            zeros = sum(!is.element(contigs, cmat$contig)) + sum(cmat$contig_total_count == 0)

            # associated
            selected.contigs = cmat$contig
            n.associated = sum(is.element(cmat$contig, selected.contigs))

            df = data.frame(cluster=cluster, anchor=anchor, n.associated=n.associated, n.no.contacts=zeros, n.total=length(contigs))
            result = rbind(result, df)
        }
    }
    save.table(result, ofn)
}
