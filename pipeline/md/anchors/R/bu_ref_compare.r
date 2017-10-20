
cluster.matrix=function(mat.ifn, ofn1, ofn2)
{
    options(stringsAsFactors=F)
    smat = read.delim(mat.ifn)
    sets1 = sort(unique(smat$set1))
    sets2 = sort(unique(smat$set2))
    N1 = length(sets1)
    N2 = length(sets2)
    smat$set1.i = match(smat$set1, sets1)
    smat$set2.i = match(smat$set2, sets2)
    m = smatrix2matrix(smat, dim=c(N1,N2), i.field="set1.i", j.field="set2.i", value.field="identity", default.value=0)
    d = max(m) - m
    hc1 = hclust(dist(d))
    hc2 = hclust(dist(t(d)))

    df1 = data.frame(set=sets1[hc1$order])
    cat(sprintf("saving result: %s\n", ofn1))
    write.table(df1, ofn1, sep="\t", quote=F, col.names=T, row.names=F)

    df2 = data.frame(set=sets2[hc2$order])
    cat(sprintf("saving result: %s\n", ofn2))
    write.table(df2, ofn2, sep="\t", quote=F, col.names=T, row.names=F)
}
