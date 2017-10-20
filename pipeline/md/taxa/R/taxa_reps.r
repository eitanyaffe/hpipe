get.children.list=function(tax) {
    children.list = list()
    for (i in 1:dim(tax)[1]) {
        parent = tax$parent_id[i]
        child = tax$tax_id[i]
        if (is.na(parent))
            next
        children.list[[as.character(parent)]] = c(children.list[[as.character(parent)]], child)
    }
    return (children.list)
}

select.reps=function(ifn, ifn.order, ifn.taxa, ifn.summary, skip.ids, skip.names, ofn, min.ratio.1, min.ratio.12)
{
    table = load.table(ifn)
    anchor.table = load.table(ifn.order)

    cat(sprintf("reading tax table: %s\n", ifn.taxa))
    tax = read.delim(ifn.taxa, quote="")

    cat(sprintf("Skipping internal nodes: %s\n", paste(tax$name[is.element(tax$tax_id, skip.ids)], "(",
                                                       tax$tax_id[is.element(tax$tax_id, skip.ids)], ")", sep="", collapse=", ")))

    skip.names = gsub('\\.', ' ', skip.names)
    cat(sprintf("Skipping internal nodes which match names: %s\n", paste(skip.names, collapse=",", sep="")))
    united.skip.ids = sort(unique(c(skip.ids, tax$tax_id[is.element(tax$name, skip.names)])))

    anchors = unique(table$anchor)

    # set children list
    children.list = get.children.list(tax)
    base.root = 1

    select.node=function(tax.id=base.root, table=table, debug=F) {
        count = ifelse(is.na(match(tax.id,table$tax_id)), 0, table$weight[match(tax.id,table$tax_id)])
        if (debug)
            cat(sprintf("%d: %s W=%f\n", tax.id, tax$name[match(tax.id,tax$tax_id)], count))
        if (!(as.character(tax.id) %in% names(children.list)))
            return(tax.id)

        children = children.list[[as.character(tax.id)]]
        child.count = NULL
        for (i in seq_along(children))
            child.count[i] = ifelse(is.na(match(children[i],table$tax_id)), 0, table$weight[match(children[i],table$tax_id)])
        if (debug)
            cat(sprintf("children: %s\n", paste(children, "(", round(child.count), ")", sep="", collapse=",")))
        df = data.frame(id=children, count=child.count)
        if (dim(df)[1] == 1)
            return (select.node(tax.id=df$id[1], table=table, debug=debug))
        df = df[order(df$count, decreasing=T),]
        df$count = round(df$count,1)
        ratio.1 = df$count[1] / sum(df$count)
        ratio.12 = df$count[1] / df$count[2] - 1
        mask = is.element(tax.id, united.skip.ids)
        if ((ratio.1 < min.ratio.1 || ratio.12 < min.ratio.12) && !mask)
            return (tax.id)
        else
            return (select.node(tax.id=df$id[1], table=table, debug=debug))
    }

    debug = F
    rep = NULL
    for (i in seq_along(anchors)) {
        anchor = anchors[i]
        ttable = table[table$anchor == anchor,]
        rep.id = select.node(base.root, ttable, debug)
        ix = match(rep.id, ttable$tax_id)
        name = tax$name[match(rep.id, tax$tax_id)]
        level = tax$level[match(rep.id, tax$tax_id)]

        pid = tax$parent_id[match(rep.id, tax$tax_id)]

        # handle CAG genomes
        is.cag = grepl("CAG", name)
        cag.name = "none"
        cag.level = "none"
        if(is.cag) {
            genus.id = tax$parent_id[match(pid, tax$tax_id)]
            if (tax$name[match(pid, tax$tax_id)] != "environmental samples")
                stop("expecting 'environmental samples' as the parent of CAG")
            cag.level = tax$level[match(genus.id, tax$tax_id)]
            cag.name = tax$name[match(genus.id, tax$tax_id)]
        }

        # add strain level
        if (level == "no rank") {
            if (tax$level[match(pid, tax$tax_id)] == "species")
                level = "strain"
        }

        rep = rbind(rep, data.frame(
            anchor.id=anchor.id(anchor,anchor.table), anchor=anchor, tax_id=rep.id,
            gene_count=ttable$gene_count[ix],
            identity=round(ttable$identity[ix],2),
            coverage=round(100*ttable$coverage[ix],2),
            level=level, name=name,
            is.cag=is.cag, cag.level=cag.level, cag.name=cag.name))
    }

    # construct tree
    order.tree=function(tax.id, selected.ids) {
        result = NULL
        if ((as.character(tax.id) %in% names(children.list))) {
            children = children.list[[as.character(tax.id)]]
            for (i in seq_along(children))
                result = c(result, order.tree(children[i], selected.ids=selected.ids))
        }
        if (is.element(tax.id, selected.ids))
            result = c(result, tax.id)
        return (result)
    }

    cat(sprintf("reading summary table: %s\n", ifn.summary))
    summary = read.delim(ifn.summary, quote="")

    rep$total_genes = rowSums(summary[match(rep$anchor, summary$anchor),-1])
    rep$frac = round(100 * rep$gene_count / rep$total,2)
    rep$lineage = tax$lineage[match(rep$tax_id, tax$tax_id)]

    z = match(rep$tax_id, order.tree(base.root, selected.ids=rep$tax_id))
    rep = rep[order(z),]
    rep$order = 1:dim(rep)[1]
    rep = rep[match(anchor.table$set,rep$anchor),]

    cat(sprintf("saving result table: %s\n", ofn))
    write.table(rep, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}

select.species=function(ifn, ifn.taxa, ifn.gene.table, level, ofn)
{
    cat(sprintf("reading tax table: %s\n", ifn.taxa))
    tax = read.delim(ifn.taxa, quote="")

    cat(sprintf("reading gene table: %s\n", ifn.gene.table))
    gene.table = read.delim(ifn.gene.table, quote="")
    names(gene.table)[1] = "accession"

    table = load.table(ifn)
    table$level = tax$level[match(table$tax_id, tax$tax_id)]

    # ids = sort(unique(result$tax_id))
    ids = sort(unique(table$tax_id[table$identity > 70 & table$gene_count > 100 & table$level == level]))

    # unique genome list
    result.genomes = NULL
    for (id in ids) {
        xtable = gene.table[gene.table$taxid == id,]
        xtable = xtable[!grepl("CAG", xtable$organism_name),]
        if (dim(xtable)[1] == 0)
            next
        fields = c("accession", "taxid", "organism_name", "infraspecific_name", "assembly_level")
        if (any(xtable$assembly_level == "Complete Genome"))
            df = xtable[match("Complete Genome", xtable$assembly_level), fields]
        else {
            if (any(xtable$assembly_level == "Chromosome"))
                df = xtable[match("Chromosome", xtable$assembly_level), fields]
            else
                df = xtable[1, fields]
        }
        result.genomes = rbind(result.genomes, df)
    }
    save.table(result.genomes, ofn)
}
