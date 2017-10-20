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

genome.table=function(ifn.species, ifn.taxa, ifn.gene.table, only.complete.threshold, ofn)
{
    species = load.table(ifn.species)
    tax = load.table(ifn.taxa, quote="")
    gene.table = load.table(ifn.gene.table, quote="")
    names(gene.table)[1] = "accession"

    children.list = get.children.list(tax)

    append.id=function(anchor, rid, id, only.complete) {
        result.count = 0
        ix = gene.table$taxid == id

        if (only.complete || sum(ix) > only.complete.threshold) {
            only.complete = T
            ix = ix & (gene.table$version_status == "latest" & gene.table$assembly_level == "Complete Genome"
                & gene.table$release_type == "Major" & gene.table$genome_rep == "Full")
        }
        if (any(ix)) {
            ttable = gene.table[ix,]
            df = data.frame(
                anchor=anchor, species.taxid=rid, taxid=id, accession=ttable$accession, name=ttable$organism_name,
                assembly_level=ttable$assembly_level, release=ttable$release_type,
                excluded=ttable$excluded_from_refseq, ftp.path=ttable$ftp_path)
#            if (!is.null(result))
#                df = df[!is.element(df$accession, result$accession),]
            # cat(sprintf("anchor: %d, rid: %d, add id: %d, name: %s, accession: %s\n", anchor, rid, id, ttable$organism_name, ttable$accession))
            result <<- rbind(result, df)
            result.count = result.count + dim(df)[1]
        }
        cids = children.list[[as.character(id)]]
        for (cid in cids)
            result.count = result.count + append.id(anchor=anchor, rid=rid, id=cid, only.complete=only.complete)
        result.count
    }

    result = NULL
    all.covered = T
    species = species[order(species$anchor),]
    for (i in 1:dim(species)[1]) {
        anchor = species$anchor[i]
        id = species$tax_id[i]
        name = species$name[i]

        count = append.id(anchor=anchor, rid=id, id=id, only.complete=F)
        cat(sprintf("anchor: %d, species: %s, taxa id: %d, number of strains: %d\n", anchor, name, id, count))
        if (count == 0)
            all.covered = F
    }
    if (!all.covered)
        cat(sprintf("warning: some anchors have no matching genomes\n"))

    save.table(result, ofn)
}
