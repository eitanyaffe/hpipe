field.count=function(x, field="gene")
{
    tt = table(x[,field])
    result = data.frame(x=names(tt), count=as.vector(tt))
    names(result)[1] = field
    result[order(result$count, decreasing=T),]
}

plot.ref.assembly=function(
    original.genes.table,
    assembly.genes.table,
    ifn.o2a,
    ifn.a2o,
    identity.threshold,
    coverage.threshold,
    fdir)
{
    options(stringsAsFactors=F)
    ogenes = read.delim(original.genes.table)

    agenes = read.delim(assembly.genes.table)

    o2a = read.delim(ifn.o2a)
    a2o = read.delim(ifn.a2o)
    o2a = o2a[o2a$identity > identity.threshold,]
    a2o = a2o[a2o$identity > identity.threshold,]
    o2a = o2a[o2a$coverage > coverage.threshold,]
    a2o = a2o[a2o$coverage > coverage.threshold,]

    # original gene stats
    ogenes$found = is.element(ogenes$gene, o2a$query)
    lost = 100*sum(!ogenes$found)/dim(ogenes)[1]

    # assembly gene stats
    s = sapply(split(o2a$query, o2a$hit), length)
    agenes$count = s[match(agenes$gene, names(s))]
    agenes$count[is.na(agenes$count)] = 0
    false = 100*sum(agenes$count==0)/dim(agenes)[1]
    collapsed = 100*sum(agenes$count>1)/dim(agenes)[1]

    df = data.frame(o_lost=lost, a_false=false, a_collapsed=collapsed)
    system(paste("mkdir -p", fdir))
    ofn = paste(fdir, "/original_genes.png", sep="")
    cat(sprintf("generating figure: %s\n", ofn))
    png(ofn, 300, 500)
    par(mai=c(2,1,1,0.5))
    barplot(as.matrix(df), main="gene breakdown", las=2, col="lightgrey")
    dev.off()
}

mark.contigs=function(
    original.genes.table,
    assembly.genes.table,
    assembly.contig.table,
    ifn.o2a,
    ifn.a2o,
    identity.threshold,
    coverage.threshold,
    ofn)
{
    options(stringsAsFactors=F)
    ogenes = read.delim(original.genes.table)
    agenes = read.delim(assembly.genes.table)
    acontigs = read.delim(assembly.contig.table)

    o2a = read.delim(ifn.o2a)
    a2o = read.delim(ifn.a2o)
    o2a = o2a[o2a$identity > identity.threshold & o2a$coverage > coverage.threshold,]
    a2o = a2o[a2o$identity > identity.threshold & a2o$coverage > coverage.threshold,]

    o2a$genome = ogenes$contig[match(o2a$source, ogenes$gene)]

    s = sapply(split(o2a$source, o2a$target), length)
    agenes$count = s[match(agenes$gene, names(s))]
    agenes$count[is.na(agenes$count)] = 0
    agenes$genome = ifelse(agenes$count == 1, o2a$genome[match(agenes$gene, o2a$target)], ifelse(agenes$count > 1, "multi", "false"))

    s = split(agenes$genome, agenes$contig)

    ss = sapply(s, function(x) {
        N = length(x)
        x = x[x != "false" & x != "multi"]
        t = table(x)
        if (length(t) == 1 && t>0.5*N)
            return (names(t)[1])
        else
            return ("multi")
    } )
    df = data.frame(contig=names(ss), genome=ss)
    df$contig_length = acontigs$length[match(df$contig, acontigs$contig)]

    cat(sprintf("saving result table: %s\n", ofn))
    write.table(df, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}
