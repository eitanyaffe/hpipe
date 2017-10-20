select.path=function(ifn.order, ifn.tree, ifn.taxa, ifn.gene.table, ifn.summary, ofn)
{
    anchor.table = load.table(ifn.order)

    cat(sprintf("reading tax table: %s\n", ifn.taxa))
    taxa = read.delim(ifn.taxa, quote="")
    taxa$level = ifelse(taxa$level == "no rank" & taxa$level[match(taxa$parent_id, taxa$tax_id)] == "species", "subspecies", taxa$level)
    taxa$level[is.na(taxa$level)] = "root"

    tree = load.table(ifn.tree)

    cat(sprintf("reading gene table: %s\n", ifn.gene.table))
    gene.table = read.delim(ifn.gene.table, quote="")
    names(gene.table)[1] = "accession"

    summary = load.table(ifn.summary)

    tree$has.genome = !is.na(match(tree$tax_id, gene.table$taxid))
    anchors = sort(unique(tree$anchor))

    result = NULL
    for (anchor in anchors) {
        ttree = tree[tree$anchor == anchor,]
        tleaf = ttree[ttree$is_leaf & ttree$has.genome,]
        tleaf = tleaf[order(tleaf$gene_count, decreasing=T),]

        if(dim(tleaf)[1] < 1)
            stop("no rep found")
        id = tleaf$tax_id[1]

        total_genes = rowSums(summary[match(anchor, summary$anchor),-1])

        df = NULL
        index = 1
        pid = taxa$parent_id[match(id, taxa$tax_id)]
        while (pid != -1) {
            name = taxa$name[match(id, taxa$tax_id)]
            level = taxa$level[match(id, taxa$tax_id)]
            pid = taxa$parent_id[match(id, taxa$tax_id)]
            identity = ttree$identity[match(id, ttree$tax_id)]
            coverage = ttree$coverage[match(id, ttree$tax_id)]
            gene.count = ttree$gene_count[match(id, ttree$tax_id)]
            has.genome = ttree$has.genome[match(id, ttree$tax_id)]
            frac = round(100 * gene.count / total_genes, 2)
            df = rbind(df,
                data.frame(anchor=anchor, index=index, tax_id=id, level=level, name=name, identity=identity, coverage=coverage,
                           gene.count=gene.count, has.genome=has.genome, frac=frac))
            id = pid
        }
        result = rbind(result, df)
    }

    cat(sprintf("saving result table: %s\n", ofn))
    write.table(result, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}
