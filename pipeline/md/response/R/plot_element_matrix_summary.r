plot.diameter.analysis=function(ifn.elements, mat.ifn, fdir)
{
    mat = load.table(mat.ifn)
    mat = mat[mat$identity > 0,]
    df = load.table(ifn.elements)
#    hist(df$diameter, col="darkgreen", border=1, las=2)
    N = 8
    breaks = c(40,50,60,70,80)
    bg = sapply(split(mat$identity, cut(mat$identity, breaks)), length)
    shared = sapply(split(df$median.identity, cut(df$identity.diameter, breaks)), length)
    bg = (bg / sum(bg)) * sum(shared)
    fig.start(fdir=fdir, ofn=paste(fdir, "/base_diameter.pdf", sep=""), type="pdf", width=2.5, height=4)
    barplot(rbind(shared, bg), beside=T, col=c("darkgreen", "gray"), border=NA, las=2)
    fig.end()
}

plot.element.matrix.summary=function(ifn.anchors, ifn.elements, ifn.majors, ifn.network, ifn.cluster.class, fdir)
{
    anchor.table = load.table(ifn.anchors)
    element.table = load.table(ifn.elements)
    network = load.table(ifn.network)
    df = load.table(ifn.cluster.class)
    df$any = rowSums(df[,c("phage", "plasmid")]) > 0
    types = names(df)[-1]
    df = df[is.element(df$cluster, element.table$cluster),]
    df$diameter = element.table$identity.diameter[match(df$cluster, element.table$cluster)]
    M = 3
    # breaks = quantile(df$diameter, 0:M/M)
    breaks = c(40,50,60,70,80)
    N = length(types)
    fig.start(fdir=fdir, ofn=paste(fdir, "/type_diameter_total.pdf", sep=""), type="pdf", width=8, height=4)
    layout(matrix(1:N, 1, N))
    for (type in types) {
        s = split(df[,type], cut(df$diameter, breaks=breaks))
        with = sapply(s, sum)
        without = sapply(s, length) - sapply(s, sum)
        barplot(rbind(with, without), col=(c("darkgreen", "gray")), main=type, las=2, border=NA)
    }
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/type_diameter_enrich.pdf", sep=""), type="pdf", width=8, height=4)
    layout(matrix(1:N, 1, N))
    for (type in types) {
        s = split(df[,type], cut(df$diameter, breaks=breaks))
        with = sapply(s, sum)
        without = sapply(s, length) - sapply(s, sum)
        barplot(100*with/(with+without), col="darkgreen", main=type, las=2)
    }
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/type_diameter_all_sum.pdf", sep=""), type="pdf", width=2, height=4)
    s = split(df$any, cut(df$diameter, breaks=breaks))
    with = sapply(s, sum)
    without = sapply(s, length) - sapply(s, sum)
    barplot(rbind(with, without), col=(c("red", "darkblue")), main=type, las=2, border=NA)
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/type_diameter_all_enrich.pdf", sep=""), type="pdf", width=2, height=4)
    s = split(df$any, cut(df$diameter, breaks=breaks))
    with = sapply(s, sum)
    without = sapply(s, length) - sapply(s, sum)
    barplot(100*with/(with+without), col="red", main="any", las=2, border=NA)
    fig.end()

}
