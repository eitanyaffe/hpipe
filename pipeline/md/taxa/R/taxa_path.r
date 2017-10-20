taxa.path=function(ifn.order, ifn.tree, ifn.taxa, ifn.gene.table, ifn.summary, ofn.path, ofn.species)
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
    tree$level = taxa$level[match(tree$tax_id, taxa$tax_id)]
    anchor.ids = anchor.table$id

    df.path = NULL
    for (anchor.id in anchor.ids) {
        anchor = anchor.table$set[match(anchor.id, anchor.table$id)]
        ttree = tree[tree$anchor == anchor,]
        tleaf = ttree[(ttree$is_leaf | (ttree$level == "species")) & ttree$has.genome,]
        tleaf = tleaf[order(tleaf$gene_count, decreasing=T),]

        if(dim(tleaf)[1] < 1)
            stop(sprintf("no rep found (anchor=%d, id=%s)", anchor, anchor.id))
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

            # fix level
            if (level == "no rank" && name == "cellular organisms")
                level = "cellular.organisms"
            if (level == "no rank" && name == "environmental samples")
                level = paste("env.", taxa$level[match(pid, taxa$tax_id)], sep="")
            if (level == "no rank" && grepl("group", name))
                level = "group"

            df = rbind(df,
                data.frame(anchor.id=anchor.id, anchor=anchor, index=index, level=level, name=name, tax_id=id, coverage=coverage,
                           gene.count=gene.count, has.genome=has.genome, identity=identity, frac=frac, parent.tax.id=pid))
            id = pid
            index = index + 1
        }
        df.path = rbind(df.path, df)
    }

    # path
    save.table(df.path, ofn.path)

    # focus on species for single rep per anchor
    tt = df.path[df.path$level == "species",]
    tt$anchor.id = anchor.table$id[match(tt$anchor, anchor.table$set)]
    tt = tt[match(anchor.table$id, tt$anchor.id),]
    df.species = data.frame(
        anchor.id=tt$anchor.id,
        anchor=tt$anchor,
	tax_id=tt$tax_id,
        gene_count=tt$gene.count,
        identity=tt$identity,
	level=tt$level,
	name=tt$name,
        frac=tt$frac,
        parent.tax.id=tt$parent.tax.id,
        is.cag=F, cag.name="none", cag.level="none")
    for (i in 1:length(anchor.ids)) {
        anchor = df.species$anchor[i]
        name = df.species$name[i]
        pid = df.species$parent.tax.id[i]
        gene_count = df.species$gene_count[i]
        is.cag = grepl("CAG", name)
        cag.name = "none"
        cag.level = "none"
        if (is.cag) {
            genus.id = taxa$parent_id[match(pid, taxa$tax_id)]
            if (taxa$name[match(pid, taxa$tax_id)] != "environmental samples")
                stop("expecting 'environmental samples' as the parent of CAG")
            cag.level = taxa$level[match(genus.id, taxa$tax_id)]
            cag.name = taxa$name[match(genus.id, taxa$tax_id)]
        }
        df.species$is.cag[i] = is.cag
        df.species$cag.level[i] = cag.level
        df.species$cag.name[i] = cag.name

        # add coverage
        gene.total = sum(summary[match(anchor, summary$anchor),-1])
        df.species$coverage[i] = round(100 * gene_count / gene.total,2)
    }
    save.table(df.species, ofn.species)
}
