compute.network=function(
    ifn.anchors, ifn.map, ifn.elements, ifn.contigs, ifn.element.info, ifn.identity,
    ifn.genes, ofn.elements, ofn.genes, ofn.network)
{
    anchor.table = load.table(ifn.anchors)
    info = load.table(ifn.element.info)
    network = load.table(ifn.elements)

    anchors = anchor.table$set
    network = network[!network$major,]

    map = load.table(ifn.map)
    map$connected = map$type == "connected"
    map$id = paste(map$cluster, map$anchor, sep="_")
    network$id = paste(network$cluster, network$anchor, sep="_")
    network$connected = map$connected[match(network$id, map$id)]
    network = network[network$connected,]

    # compute min correlation of each cluster
    smin = sapply(split(network$pearson, network$cluster), min)
    smax = sapply(split(network$pearson, network$cluster), max)
    network = network[!is.element(network$cluster, anchors),]
    network = network[,c("anchor", "cluster", "pearson")]

    elements = sort(unique(network$cluster))

    table.iden = load.table(ifn.identity)
    table.iden$anchor1 = table.iden$set1
    table.iden$anchor2 = table.iden$set2
    table.iden = table.iden[,c("anchor1" ,"anchor2", "identity")]
    table.iden = table.iden[table.iden$anchor1 != "NONE" & table.iden$anchor2 != "NONE",]

    genes = load.table(ifn.genes)
    genes = genes[is.element(genes$cluster, elements),]

    # compute identity diameter
    cat(sprintf("computing identity diameter...\n"))
    s = split(network$anchor, network$cluster)
    element.table = NULL
    for (i in 1:length(s)) {
        element = names(s)[i]
        hosts = s[[i]]
        gg = expand.grid(hosts, hosts)
        gg = gg[gg[,1] != gg[,2],]
        min.identity = 100
        if (dim(gg)[1] > 0) {
            for (j in 1:dim(gg)[1]) {
                identity = table.iden$identity[table.iden$anchor1 == gg[j,1] & table.iden$anchor2 == gg[j,2]]
                min.identity = min(min.identity, identity)
            }
        }
        df = data.frame(
            cluster=element,
            host.count=length(hosts),
            hosts=paste(anchor.table$id[match(hosts, anchor.table$set)], sep="", collapse=","),
            identity.diameter=min.identity)
        element.table = rbind(element.table, df)
    }
    element.table$min.pearson = smin[match(element.table$cluster, names(smin))]
    element.table$max.pearson = smax[match(element.table$cluster, names(smax))]
    element.table$gene.count = info$gene.count[match(element.table$cluster, info$cluster)]
    element.table$gene.count.uniref = info$gene.count.uniref[match(element.table$cluster, info$cluster)]
    element.table$novel.score = info$novel.score[match(element.table$cluster, info$cluster)]
    element.table$median.identity = info$median.identity[match(element.table$cluster, info$cluster)]

    element.table = element.table[,]

    genes = lookup.append2(genes, element.table, lookup.field1="cluster", lookup.field2="cluster", value.field="identity.diameter")
    genes = lookup.append2(genes, element.table, lookup.field1="cluster", lookup.field2="cluster", value.field="gene.count")
    genes = lookup.append2(genes, element.table, lookup.field1="cluster", lookup.field2="cluster", value.field="host.count")
    genes = lookup.append2(genes, element.table, lookup.field1="cluster", lookup.field2="cluster", value.field="hosts")

    save.table(element.table, ofn.elements)
    save.table(genes, ofn.genes)
    save.table(network, ofn.network)
}

get.element.index=function(anchor.table, network, clusters)
{
    anchor.ids = anchor.table$id
    network$anchor.id = anchor.table$id[match(network$anchor, anchor.table$set)]
    network$anchor.index = match(network$anchor.id, anchor.ids)

    # sort by first host
    smin = sapply(split(network$anchor.index, network$cluster), min)
    smax = sapply(split(network$anchor.index, network$cluster), max)
    if (!all(names(smin) == names(smax)))
        stop("SAD")
    df = data.frame(cluster=names(smin), min=smin, max=smax)
    df = df[order(df$min, df$max, decreasing=F),]
    match(clusters, df$cluster)
}

select.network=function(
    ifn.anchors, ifn.network, ifn.genes, ifn.elements,
    min.hosts, min.genes, max.diameter, min.uniref.genes, max.pearson,
    ofn.network, ofn.genes, ofn.elements)
{
    anchor.table = load.table(ifn.anchors)
    network = load.table(ifn.network)
    elements = load.table(ifn.elements)
    genes = load.table(ifn.genes)

    elements = elements[
        elements$max.pearson <= max.pearson &
        elements$identity.diameter <= max.diameter &
        elements$host.count >= min.hosts &
        elements$gene.count >= min.genes &
        elements$gene.count.uniref >= min.uniref.genes,]

    genes = genes[is.element(genes$cluster, elements$cluster),]
    network = network[is.element(network$cluster, elements$cluster),]
    elements = elements[is.element(elements$cluster, network$cluster),]

    elements$id = get.element.index(anchor.table=anchor.table, network=network, clusters=elements$cluster)
    elements = elements[order(elements$id),]
    elements = elements[,c(dim(elements)[2],1:(dim(elements)[2]-1))]

    # add element index to genes
    genes = cbind(data.frame(element.id=elements$id[match(genes$cluster, elements$cluster)]), genes[,-1], data.frame(cluster=genes$cluster))
    genes = genes[order(genes$element.id, genes$contig),]

    save.table(network, ofn.network)
    save.table(elements, ofn.elements)
    save.table(genes, ofn.genes)
}

network.coord=function(ifn, type, anchor.order.ifn, ofn)
{
    library(igraph)

    mat = load.table(ifn)
    e.ids = unique(mat$cluster)
    a.ids = load.table(anchor.order.ifn)$set
    a.ids = a.ids[is.element(a.ids, unique(mat$anchor))]

    vs.anchors = data.frame(id=paste("g", a.ids, sep=""), type="anchor", sid=a.ids)
    vs.elements = data.frame(id=e.ids, type="element", sid=e.ids)
    vs = rbind(vs.anchors, vs.elements)

    N.anchors = dim(vs.anchors)[1]
    N.elements = dim(vs.elements)[1]

    es = data.frame(from=paste("g", mat$anchor, sep=""), to=mat$cluster)
    es = es[is.element(es$from, vs$id) & is.element(es$to, vs$id),]

    g = graph_from_data_frame(es, directed=F, vertices=vs)

    ll = switch(type,
        drl = layout_with_fr(graph=g),
        graphopt = layout_with_graphopt(g),
        stop(paste("unknown network layout type:", type)))


    vs$x = ll[,1]
    vs$y = ll[,2]
    vs = vs[,c("id", "x", "y")]
    save.table(vs, ofn)
}

network.compare.coord=function(ifn1, ifn2, ofn)
{
    library(igraph)
    network1 = load.table(ifn1)
    network2 = load.table(ifn2)
    a.ids = sort(unique(c(network1$anchor, network2$anchor)))
    e.ids = sort(unique(c(network1$cluster, network2$cluster)))
    network = rbind(network1, network2)
    network = network[,c("anchor", "cluster")]
    network = unique(network)

    vs.anchors = data.frame(id=paste("g", a.ids, sep=""), type="anchor", sid=a.ids)
    vs.elements = data.frame(id=e.ids, type="element", sid=e.ids)
    vs = rbind(vs.anchors, vs.elements)

    es = data.frame(from=paste("g", network$anchor, sep=""), to=network$cluster)
    es = es[is.element(es$from, vs$id) & is.element(es$to, vs$id),]

    g = graph_from_data_frame(es, directed=F, vertices=vs)
    ll = layout_with_fr(g)
    vs$x = ll[,1]
    vs$y = ll[,2]
    vs = vs[,c("id", "x", "y")]
    save.table(vs, ofn)
}
