cluster.contigs=function(ifn, thresholds, odir)
{
    library(igraph)

    tt = load.table(ifn)
    contigs = tt$contig
    m = as.matrix(tt[,-1])
    m = t(t(m)/colSums(m))
    m = log10((1 + m[,-1]) / (1+m[,1]))

    cc = cor(t(m))
    system(paste("mkdir -p", odir))
    for (th in thresholds) {
        am = cc>th
        diag(am) = F
        rs = rowSums(am)
        ig = graph_from_adjacency_matrix(am, mode="undirected")
        co = components(ig)
    }
}
