
field.count=function(x, field="gene")
{
    tt = table(x[,field])
    result = data.frame(x=names(tt), count=as.vector(tt))
    names(result)[1] = field
    result[order(result$count, decreasing=T),]
}

smatrix2matrix=function(smatrix, dim, i.field="i", j.field="j", value.field="value", default.value=0)
{
  indices = smatrix[,i.field] + (smatrix[,j.field]-1) * dim[1]
  v = rep(default.value, dim[1]*dim[2])
  v[indices] = smatrix[,value.field]
  matrix(v, dim[1], dim[2])
}

plot.ga.matrix=function(ifn.ga, fdir)
{
    options(stringsAsFactors=F)

    cat(sprintf("reading table: %s\n", ifn.ga))
    table = read.delim(ifn.ga)
    df = field.count(table, "gene")
    table$count = df$count[match(table$gene, df$gene)]

    ga = table[table$count > 1, c("gene", "anchor")]
    anchors = sort(unique(table$anchor))
    genes = sort(unique(ga$gene))
    ga$anchor.i = match(ga$anchor, anchors)
    ga$gene.i = match(ga$gene, genes)
    ga$value = 1

    N = length(genes)
    M = length(anchors)
    m = smatrix2matrix(smatrix=ga, dim=c(N,M), i.field="gene.i", j.field="anchor.i", value.field="value")

    system(paste("mkdir -p", fdir))
    ofn = paste(fdir, "/ga.png", sep="")
    cat(sprintf("generating figure: %s\n", ofn))
    png(ofn, 800, 800)
    heatmap(m, labRow=F)
    dev.off()
}

plot.gene.graph=function(ifn.ga, ifn.genes, fdir)
{
    library(igraph)
    options(stringsAsFactors=F)

    cat(sprintf("reading table: %s\n", ifn.ga))
    table = read.delim(ifn.ga)

    cat(sprintf("reading table: %s\n", ifn.genes))
    gene.table = read.delim(ifn.genes)
    df = field.count(table, "gene")
    genes = gene.table$gene

    anchors = sort(unique(table$anchor))

    s = split(table$gene, table$anchor)
    N = length(s)
    edges = NULL
    for (i in 1:N) {
        t = expand.grid(s[[i]],names(s)[i])
        edges = rbind(edges, t)
    }
    vs = data.frame(
        id=c(anchors, genes),
        type=c(rep("anchor", length(anchors)), rep("gene", length(genes))),
        count=c(rep(0, length(anchors)), df$count[match(genes,df$gene)]))
    vs$gene.anchor = table$anchor[match(vs$id, table$gene)]
    vs$count[is.na(vs$count)] = 0
    vs$col = ifelse(vs$type == "anchor", NA, ifelse(vs$count == 1, vs$gene.anchor, ifelse(vs$count == 0, "white", "black")))
    vs$cex = ifelse(vs$type == "anchor", NA, 1)

    cat(sprintf("generating graph with %d vertices and %d edges\n", dim(vs)[1], dim(edges)[1]))
    gg = graph_from_data_frame(as.matrix(edges), directed=FALSE, vertices=vs$id)

    cat(sprintf("computing layout...\n"))
    lo = layout_nicely
    lo = layout_with_dh
    lo = layout_with_lgl
    lo = layout_with_kk
    cc = as.data.frame(lo(gg))

    vs$x = cc[,1]
    vs$y = cc[,2]
    es = data.frame(
        gene=edges[,1], anchor=edges[,2],
        x1=vs$x[match(edges[,1], vs$id)],
        x2=vs$x[match(edges[,2], vs$id)],
        y1=vs$y[match(edges[,1], vs$id)],
        y2=vs$y[match(edges[,2], vs$id)])
    es$col = "lightgrey"

    system(paste("mkdir -p", fdir))
    ofn = paste(fdir, "/gene_graph.png", sep="")
    cat(sprintf("plotting figure: %s\n", ofn))
    png(ofn, width=800, height=800)

    pvs = vs[vs$count > 0,]
    pes = es[is.element(es$gene, vs$id),]
    plot(pvs$x, pvs$y, type="p", col=pvs$col, pch=".", cex=pvs$cex)
    segments(pes$x1, pes$y1, pes$x2, pes$y2, col=pes$col)
    points(pvs$x, pvs$y, col=pvs$col, pch=".", cex=pvs$cex)
    # points(pvs$x, pvs$y, cex=pvs$cex)

    dev.off()

}

plot.gene.graph2=function(ifn, fdir)
{
    library(igraph)
    options(stringsAsFactors=F)

    table = read.delim(ifn)
    df = field.count(table, "gene")
    singles = df[df$count == 1, "gene"]
    shared = df[df$count > 1, "gene"]

    table.singles = table[is.element(table$gene, singles),]
    table.shared = table[is.element(table$gene, shared),]

    x = sapply(split(table.singles$gene, table.singles$anchor), length)
    agenes = data.frame(anchor=names(x), count=x)
    anchors = sort(unique(table$anchor))

    vertices = c(anchors, shared)
    type = c(rep("anchor", length(anchors)), rep("gene", length(shared)))
    edges = cbind(table.shared$gene, table.shared$anchor)

    gg = graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
    cc = as.data.frame(layout_nicely(gg))

    vs = data.frame(name=vertices, type=type, x=cc[,1], y=cc[,2])
    vs$col = ifelse(vs$type == "anchor", 1, 2)
    vs$cex = ifelse(vs$type == "anchor", 2, 0.5)

    es = data.frame(
        x1=vs$x[match(edges[,1], vertices)],
        x2=vs$x[match(edges[,2], vertices)],
        y1=vs$y[match(edges[,1], vertices)],
        y2=vs$y[match(edges[,2], vertices)])
    es$col = "lightgrey"

    system(paste("mkdir -p", fdir))
    ofn = paste(fdir, "/gene_graph.png", sep="")
    cat(sprintf("plotting figure: %s\n", ofn))
    png(ofn, width=800, height=800)

    plot.new()
    plot(vs$x, vs$y, type="p", col=NA, pch=19)
    segments(es$x1, es$y1, es$x2, es$y2, col=es$col)
    points(vs$x, vs$y, col=vs$col, pch=19, cex=vs$cex)

    dev.off()
}
