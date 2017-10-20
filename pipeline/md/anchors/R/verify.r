verify.ga=function(ifn, genes.ifn, original.genome.table, original.genes.table, o2a.ifn, a2o.ifn, fdir)
{
    options(stringsAsFactors=F)
    ga = read.delim(ifn)
    genes = read.delim(genes.ifn)
    a2o = read.delim(a2o.ifn)
    o2a = read.delim(o2a.ifn)

    # assign max

    tt = table(o2a$hit)
    df = data.frame(gene=names(tt), count=as.numeric(tt))

    o2a$count = df$count[match(o2a$hit, df$gene)]
    o2a$count[is.na(o2a$count)] = 0
    o2a$anchor = ga$anchor[match(o2a$hit, ga$gene)]
    o2a$anchor[is.na(o2a$anchor)] = "none"

    ogenes = read.delim(original.genes.table)
    ogenomes = read.delim(original.genome.table)
    ogenes$genome = ogenomes$genome[match(ogenes$contig, ogenomes$contig)]

    s = split(ogenes$gene, ogenes$genome)
    p = sapply(s, function(x) table(o2a$anchor[match(x, o2a$query)]))
}

