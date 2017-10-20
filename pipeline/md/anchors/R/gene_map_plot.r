map.plot=function(
    cgenes.ifn1, cgenes.ifn2, set.field1,
    matrix.ifn1, matrix.ifn2, set.field2,
    map.1to2.ifn, map.2to1.ifn, fdir)
{
    options(stringsAsFactors=F)

    cgenes1 = read.delim(cgenes.ifn1)
    mat1 = read.delim(matrix.ifn1)

    cgenes2 = read.delim(cgenes.ifn2)
    mat2 = read.delim(matrix.ifn2)

    map.1to2 = read.delim(map.1to2.ifn)
    map.2to1 = read.delim(map.2to1.ifn)
}

xx=function()
{
    ogenes = read.delim(original.genes.table)
    ogenomes = read.delim(original.genome.table)

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
