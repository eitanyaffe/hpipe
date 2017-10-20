short.rank=function(ranks) {
    ranks = sub(" ", ".", ranks)
    rank.list = list(strain="sbS", class="C", family="F", forma="FO", genus="G", infraclass="iC", infraorder="iO", kingdom="K", no.rank="nr", order="O", parvorder="pO", phylum="P", species="S", species.group="Sg" , species.subgroup="Ssg" , subclass="sbC", subfamily="sbF", subgenus="sbG", subkingdom="sbK", suborder="sbO", subphylum="sbP", subspecies="sbS", subtribe="sbT", superclass="spC", superfamily="spF", superkingdom="spK", superorder="spO", tribe="T", varietas="V")
    rank.map = data.frame(rank=names(rank.list), short=unlist(rank.list))
    ix = match(ranks, rank.map$rank)
    ifelse(!is.na(ix), rank.map$short[ix], "X")
}

taxa.resolve=function(ifn.order, ifn.tree, ifn.path, ifn.taxa, ifn.summary, min.coverage, min.ratio, ofn)
{
    anchor.table = load.table(ifn.order)
    table = load.table(ifn.tree)
    path = load.table(ifn.path)
    summary = load.table(ifn.summary)

    root.id = 1
    cat(sprintf("reading tax table: %s\n", ifn.taxa))
    taxa = read.delim(ifn.taxa, quote="")
    taxa$level = ifelse(taxa$level == "no rank" & taxa$level[match(taxa$parent_id, taxa$tax_id)] == "species", "subspecies", taxa$level)
    taxa$level[is.na(taxa$level)] = "root"
    nodes = tree.node(id=taxa$tax_id, root.id=root.id, parent.id=taxa$parent_id, name=taxa$name, level=taxa$level)

    ids = anchor.table$id
    result = NULL
    for (anchor.id in ids) {
        anchor = anchor.table$set[match(anchor.id, anchor.table$id)]
        atable = table[table$anchor == anchor,]
        anodes = nodes[is.element(nodes$id, atable$tax_id),]
        gene.total = sum(summary[match(anchor, summary$anchor),-1])
        atable$coverage = 100 * atable$gene_count / gene.total
        anodes$coverage = atable$coverage[match(anodes$id, atable$tax_id)]
        atree = make.tree(anodes)
        apath = path[path$anchor.id == anchor.id,]

        mid.path = F
        type = "ok"
        depth = 1
        violating.id = -1
        cat(".")
        for (path.i in dim(apath)[1]:2) {
            current.id = apath$tax_id[path.i]
            next.id = apath$tax_id[path.i-1]
            cids = get.children(atree, current.id)
            cnodes = anodes[is.element(anodes$id, cids),]
            cnodes = cnodes[order(cnodes$coverage, decreasing=T),]
            children.count = dim(cnodes)[1]
            if (children.count == 1)
                next
            ix = match(next.id, cnodes$id)
            if (is.na(ix)) {
                stop("cannot find path id in children")
            }
            if (ix != 1) {
                type = "ambiguous.children"
                mid.path = T
                depth = 1
                violating.id = next.id
                break
            }
            if (cnodes$coverage[1] < min.coverage) {
                type = "low.coverage"
                mid.path = T
                depth = 0
                violating.id = next.id
                break
            }
            aratio = 100 * cnodes$coverage[1] / cnodes$coverage[2]
            if (aratio < min.ratio) {
                type = "ambiguous.children"
                mid.path = T
                depth = 1
                violating.id = next.id
                break
            }
        }
        if (mid.path) {
            rid = current.id
        } else {
            rid = next.id
        }
        ix = match(rid, apath$tax_id)

        # go up if no rank
        if (apath$level[ix] == "no rank") {
            rid = apath$parent.tax.id[ix]
            depth = depth + 1
            ix = match(rid, apath$tax_id)
        }

        df = data.frame(anchor.id=anchor.id, anchor=anchor, tax.id=rid, violating.id=violating.id,
            mid.path=mid.path, type=type, identity=apath$identity[ix], coverage=atable$coverage[match(rid, atable$tax_id)],
            level=apath$level[ix], short.level=short.rank(apath$level[ix]), name=apath$name[ix], depth=depth)
        result = rbind(result, df)
    }
    cat("\n")
    save.table(result, ofn)
}
