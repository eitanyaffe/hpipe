info=function(response.ifn, clusters.ifn, coverage.ifn, contig.ifn, gc.ifn, gene.ifn, uniref.ifn, ca.ifn, ofn)
{
    my.apply=function(s, f, name) {
        ss = sapply(s, f)
        result = data.frame(id=names(ss), value=ss)
        names(result)[1] = name
        result
    }
    table = load.table(clusters.ifn)
    contigs = load.table(contig.ifn)
    gc = load.table(gc.ifn)
    coverage = load.table(coverage.ifn)
    response = load.table(response.ifn)
    gene = load.table(gene.ifn)
    uniref = load.table(uniref.ifn)
    ca = load.table(ca.ifn)
    anchors = sort(unique(ca$anchor))

    # contig count
    result = field.count(table, "cluster")
    result$cluster = as.numeric(result$cluster)
    result = result[order(result$cluster),]
    result = result[result$cluster != -1,]
    names(result)[2] = "contig.count"

    # contigs
    contigs = lookup.append(table=contigs, lookup.table=table, lookup.field="contig", value.field="cluster")
    s = split(contigs$length, contigs$cluster)
    df = my.apply(s, function(x) sum(x), "cluster")
    result$total.bp = lookup(table=result, lookup.table=df, lookup.field="cluster", value.field="value")
    df = my.apply(s, function(x) round(median(x)), "cluster")
    result$median.contig.length = lookup(table=result, lookup.table=df, lookup.field="cluster", value.field="value")

    # genes
    gene = lookup.append(table=gene, lookup.table=table, lookup.field="contig", value.field="cluster")
    s = split(gene$gene, gene$cluster)
    df = my.apply(s, function(x) length(x), "cluster")
    result$gene.count = lookup(table=result, lookup.table=df, lookup.field="cluster", value.field="value", na.value=0)

    # uniref gene
    uniref = lookup.append(table=uniref, lookup.table=gene, lookup.field="gene", value.field="contig")
    uniref = lookup.append(table=uniref, lookup.table=table, lookup.field="contig", value.field="cluster")

    s = split(uniref$identity, uniref$cluster)
    df = my.apply(s, function(x) median(x), "cluster")
    result$median.identity = lookup(table=result, lookup.table=df, lookup.field="cluster", value.field="value", na.value=0)

    # novel genes
    s = split(uniref$prot_desc, uniref$cluster)
    df = my.apply(s, function(x) {
        length(x) - sum(grepl("uncharacterized protein", x, ignore.case=T) | grepl("hypothetical protein", x, ignore.case=T))
    }, "cluster")
    result$gene.count.uniref = lookup(table=result, lookup.table=df, lookup.field="cluster", value.field="value", na.value=0)
    result$novel.score = ifelse(result$gene.count>0, 1 - result$gene.count.uniref /  result$gene.count, 0)
    result$novel.score = round(100 * result$novel.score,1)

    # gc-content
    gc = lookup.append(table=gc, lookup.table=table, lookup.field="contig", value.field="cluster")
    s = split(gc$gc, gc$cluster)
    df = my.apply(s, function(x) round(median(x),3), "cluster")
    result$gc = round(lookup(table=result, lookup.table=df, lookup.field="cluster", value.field="value", na.value=0),2)

    # pearson radius
    cat(sprintf("computing pearson over contigs: %d\n", dim(contigs)[1]))
    response = response[is.element(response$contig, contigs$contig),]
    rownames(response) = response$contig
    cc = cor(t(response[,-1]))
    cc[is.na(cc) | !is.finite(cc)] = -1
    s = split(contigs$contig, contigs$cluster)
    result$min.intra.pearson = 1
    for (i in 1:length(s)) {
        scontigs = s[[i]]
        cluster = names(s)[i]
        ix = match(scontigs, rownames(cc))
        scc = cc[ix,ix]
        result$min.intra.pearson[match(cluster,result$cluster)] = min(scc)
    }

    # abundance
    coverage = lookup.append(table=coverage, lookup.table=table, lookup.field="contig", value.field="cluster")
    s = split(coverage$abundance.enrichment, coverage$cluster)
    df = my.apply(s, function(x) {max(x) - min(x)}, "cluster")
    result$abundance.range = round(lookup(table=result, lookup.table=df, lookup.field="cluster", value.field="value"),2)
    df = my.apply(s, function(x) {max(x)}, "cluster")
    result$max.abundance = round(lookup(table=result, lookup.table=df, lookup.field="cluster", value.field="value"),2)
    df = my.apply(s, function(x) {min(x)}, "cluster")
    result$min.abundance = round(lookup(table=result, lookup.table=df, lookup.field="cluster", value.field="value"),2)
    df = my.apply(s, function(x) {median(x)}, "cluster")
    result$median.abundance = round(lookup(table=result, lookup.table=df, lookup.field="cluster", value.field="value"),2)

    result$is.anchor = is.element(result$cluster, anchors)
    save.table(result, ofn)
}


classify=function(clusters.ifn, contig.ifn, gene.ifn, uniref.ifn, ca.ifn, ofn)
{
    my.apply=function(s, f, name) {
        ss = sapply(s, f)
        result = data.frame(id=names(ss), value=ss)
        names(result)[1] = name
        result
    }

    options(stringsAsFactors=F)

    table = load.table(clusters.ifn)
    contigs = load.table(contig.ifn)
    gene = load.table(gene.ifn)
    uniref = load.table(uniref.ifn)
    ca = load.table(ca.ifn)

    # contig count
    result = field.count(table, "cluster")
    result$cluster = as.numeric(result$cluster)
    result = result[order(result$cluster),]
    result = result[result$cluster != -1,]
    names(result)[2] = "contig.count"

    # contigs
    contigs = lookup.append(table=contigs, lookup.table=table, lookup.field="contig", value.field="cluster")

    # genes
    gene = lookup.append(table=gene, lookup.table=table, lookup.field="contig", value.field="cluster")

    # uniref gene
    uniref = lookup.append(table=uniref, lookup.table=gene, lookup.field="gene", value.field="contig")
    uniref = lookup.append(table=uniref, lookup.table=table, lookup.field="contig", value.field="cluster")

    clusters = sort(unique(table$cluster))
    N = length(clusters)

    classes = c("phage", "plasmid", "transposase", "defense", "recombination")
    M = length(classes)

    classes.protein = list(
        phage=c("phage", "holin", "integrase", "excisionase"),
        plasmid=c("parb", "plasmid", "traa", "trse", "moba", "conju"),
        transposase="transposase",
        defense=c("nuclease", "methyltransferase", "DNA methylase"),
        recombination=c("recombinase", "resolvase"))
    if (!all(names(classes.protein) == classes))
        stop("internal error")

    classes.taxa = list(
        phage=c("phage"),
        plasmid=c("plasmid"),
        transposase=NULL,
        defense=NULL,
        recombination=NULL)

    if (!all(names(classes.taxa) == classes))
        stop("internal error")

    result = matrix(0, N, M)
    colnames(result) = classes

    s = split(uniref$prot_desc, uniref$cluster)
    for (i in 1:length(classes)) {

        # proteins
        for (pattern in classes.protein[[i]]) {
            x = sapply(s, function(x) { any(grepl(pattern=pattern, x=x, ignore.case=T)) })
            xx = data.frame(cluster=as.numeric(names(x)), value=x)
            result[match(xx$cluster, clusters),i] = result[match(xx$cluster, clusters),i] | xx$value
        }

        # taxa
        for (pattern in classes.taxa[[i]]) {
            x = sapply(s, function(x) { any(grepl(pattern=pattern, x=x, ignore.case=T)) })
            xx = data.frame(cluster=as.numeric(names(x)), value=x)
            result[match(xx$cluster, clusters),i] = result[match(xx$cluster, clusters),i] | xx$value
        }
    }
    save.table(data.frame(cluster=clusters, result), ofn)
}

print.info=function(ifn.info, odir)
{
    info = load.table(ifn.info)
    info = info[!info$is.anchor,]
    lines = NULL
    lines = c(lines, sprintf("number of elements: %d", dim(info)[1]))
    lines = c(lines, sprintf("total length: %f", sum(info$total.bp)))
    lines = c(lines, sprintf("median length: %f", median(info$total.bp)))

    save.lines(odir=odir, ofn=paste(odir,"/element_info.txt",sep=""), lines)
}
