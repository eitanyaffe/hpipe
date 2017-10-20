map.plot=function(
    genes.ifn1, genes.ifn2, set.field1,
    matrix.ifn1, matrix.ifn2, set.field2,
    map.ifn, fdir)
{
    options(stringsAsFactors=F)

    genes1 = read.delim(genes.ifn1)
    mat1 = read.delim(matrix.ifn1)
    mat1$set = mat1[,set.field1]
    set1 = mat1[,c("cgene", "set")]

    genes2 = read.delim(genes.ifn2)
    mat2 = read.delim(matrix.ifn2)
    mat2$set = mat2[,set.field2]
    set2 = mat2[,c("cgene", "set")]

    map = read.delim(map.ifn)

    map$cgene1 = genes1$cgene[match(map$source, genes1$gene)]
    map$cgene2 = genes2$cgene[match(map$target, genes2$gene)]

    s1 = split(set1$cgene, set1$set)
    s2 = split(set2$cgene, set2$set)

    result = NULL
    for (i1 in 1:length(s1))
    for (i2 in 1:length(s2)) {
        x1 = s1[[i1]]
        x2 = s2[[i2]]

        N1 = length(x1)
        N2 = length(x2)

        tmap = map[is.element(map$cgene1, x1),]
    }
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
