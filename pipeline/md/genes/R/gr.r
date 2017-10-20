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

get.taxa.table=function(ifn.taxa, ifn.tree, ifn.genebank, ofn)
{
    taxa = load.table(ifn.taxa)
    tree = load.table(ifn.tree, quote="")
    genebank = load.table(ifn.genebank, quote="")
    names(genebank)[1] = "accession"

    ids = sort(unique(taxa$id))
    result = NULL
    for (i in 1:length(ids)) {
        rid = ids[i]
        rname = tree$name[tree$tax_id == rid]

        # species taxa ids
        sids = tree$tax_id[tree$root_id == rid & tree$level == "species"]

        # accession ids
        ix = is.element(genebank$taxid, sids) & genebank$excluded_from_refseq == ""
        if (any(ix)) {
            ttable = genebank[ix,]
            df = data.frame(
                root.taxid=rid, root.name=rname, taxid=ttable$taxid, species_taxid=ttable$species_taxid,
                accession=ttable$accession, name=ttable$organism_name,
                assembly_level=ttable$assembly_level, release=ttable$release_type,
                excluded=ttable$excluded_from_refseq, ftp.path=ttable$ftp_path, bioproject=ttable$bioproject, biosample=ttable$biosample)
            result = rbind(result, df)
        }
    }
    result = result[!duplicated(result$accession),]
    if (is.null(result))
        stop("non found")

    save.table(result, ofn)
}

get.gold.table=function(ifn.gold.all, ifn.gold.select, ifn.genebank, ofn)
{
    gold = load.table(ifn.gold.all)
    gold.select = load.table(ifn.gold.select)
    genebank = load.table(ifn.genebank, quote="")
    names(genebank)[1] = "accession"
    # genebank = genebank[genebank$excluded_from_refseq != "",]

    gold.ids = gold.select$GOLD.Project.ID[gold.select$Domain == "Bacteria"]
    gold = gold[is.element(gold$GOLDSTAMP, gold.ids),]
    gold$id = gold$NCBI.TAXON.ID
    gold = gold[!duplicated(gold$id),]

    ids = sort(unique(gold$id))
    x = genebank[is.element(genebank$taxid, ids),]

    result = NULL
    for (i in 1:length(ids)) {
        rid = ids[i]

        # accession ids
        ix = genebank$species_taxid == rid
        if (any(ix & genebank$assembly_level == "Complete Genome"))
            ix = ix & genebank$assembly_level == "Complete Genome"

        if (any(ix)) {
            ttable = genebank[ix,]

            # !!! for now take only first genome
            ttable = ttable[1,]

            df = data.frame(
                root.taxid=rid, name=ttable$organism_name, taxid=ttable$taxid, species_taxid=ttable$species_taxid,
                accession=ttable$accession, name=ttable$organism_name,
                assembly_level=ttable$assembly_level, release=ttable$release_type,
                excluded=ttable$excluded_from_refseq, ftp.path=ttable$ftp_path, bioproject=ttable$bioproject, biosample=ttable$biosample)
            result = rbind(result, df)
        }
    }
    result = result[!duplicated(result$accession),]
    if (is.null(result))
        stop("non found")

    save.table(result, ofn)
}

sample.taxa.table=function(ifn, seed, max.genomes, ofn)
{
    set.seed(seed)
    df = load.table(ifn)
    N = dim(df)[1]
    N.select = min(N, max.genomes)
    ix = sample(x=1:N, size=N.select, replace=F)
    result = df[ix,]
    save.table(result, ofn)
}

shared.distrib=function(ifn.genes, ifn.interference, ifn.genomes, ifn.taxa, ifn.lookup, ofn.table, ofn.distrib)
{
    genes = load.table(ifn.genes)
    inter = load.table(ifn.interference)
    genomes = load.table(ifn.genomes)
    taxa = load.table(ifn.taxa)
    lookup = load.table(ifn.lookup)

    genes$id = genomes$id[match(genes$genome, genomes$set)]

    # remove gene cluster if any gene is intererred
    genes$inter = is.element(genes$gene, inter$gene)
    s = sapply(split(genes$inter, genes$cgene), any)
    genes = genes[!s[match(genes$cgene, names(s))],]

    sf = sapply(split(genes$id, genes$cgene), unique)

    result = as.data.frame(table(sapply(sf, length)))
    names(result) = c("multiplicity", "count")
    save.table(result, ofn.distrib)

    return (NULL)

    s = sf[sapply(sf, length) > 2]
    sx = sapply(s, function(x) { paste(x, collapse="_")} )
    tt = sort(table(sx), decreasing=T)
    xx = data.frame(key=names(tt), count=tt)
    save.table(xx, ofn.table)

    ids = genomes$id
    ids.genomes = genomes$set[match(ids,genomes$id)]
    ids.accession = lookup$genome_id[match(ids.genomes,lookup$genome)]
    ids.bioproject = taxa$bioproject[match(ids.accession,taxa$accession)]
    ids.label = match(ids.bioproject, unique(ids.bioproject))

    N = length(ids)
    mm = matrix(rep(0,N^2), N, N)
    for (i in 1:length(tt)) {
        title = names(tt)[i]
        count = tt[i]
        indices = match(unlist(strsplit(title, "_")), ids)
        for (ix in indices)
        for (iy in indices) {
            if (ix < iy)
                mm[ix,iy] = mm[ix,iy] + count
        }
    }
    sm = matrix2smatrix(mm)
    sm = sm[sm$value > 0,]
    colors = c("gray", "blue", "red", "orange")
    breaks = c(1, 3, 10, 100)
    panel = make.color.panel(colors)
    sm$color = panel[vals.to.cols(vals=sm$value, breaks=breaks)]
    plot.new()
    plot.window(xlim=c(0,N), ylim=c(0,N))
    abline(h=0:N, col="gray")
    abline(v=0:N, col="gray")
    rect(sm$i-1, sm$j-1, sm$i, sm$j, col=sm$color, border=NA)
    rect(sm$j-1, sm$i-1, sm$j, sm$i, col=sm$color, border=NA)
    axis(1, at=1:N-0.5,ids.label, las=2)
    axis(2, at=1:N-0.5,ids.bioproject, las=2)
}
