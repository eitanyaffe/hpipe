plot.shared.genes=function(ca.matrix.ifn, ifn.order, contigs.ifn, genes.ifn, uniref.ifn)
{
    gene.table = load.table(genes.ifn)
    uniref = load.table(uniref.ifn)

    mat = load.table(ca.matrix.ifn)
    ca = load.table(contigs.ifn)
    df = field.count(ca, "contig")
    contigs = df$contig[df$count > 2]

    genes = gene.table$gene[is.element(gene.table$contig, contigs)]
    u = uniref[is.element(uniref$gene, genes),]

    anchors = load.table(ifn.order)$set
    Na = length(anchors)
    Nc = length(contigs)

    mx = mat[is.element(mat$contig, contigs),]
    smat = data.frame(ci=match(mx$contig, contigs), ai=match(mx$anchor, anchors), value=log10(mx$contig_total_count/mx$contig_expected))

    dim = c(Nc, Na)
    m = smatrix2matrix(smat, dim, i.field="ci", j.field="ai", value.field="value")

    method = "average"
    d = dist(m)
    hc = hclust(d, method=method)
    o.contigs = hc$order

    method = "average"
    d = dist(t(m))
    hc = hclust(d, method=method)
    o.anchors = hc$order

    ncols = 256
    colors = c("white", "blue", "orange")
    panel = make.color.panel(colors, ncols)
    breaks = c(-1, 1, 4)

    mo = m[o.contigs,o.anchors]
    sm = matrix2smatrix(mo)
    sm$col = panel[vals.to.cols(sm$value, breaks, ncols)]

    system(paste("mkdir -p", fdir))
    fig.fn = paste(fdir, "/", title, "_matrix.png", sep="")
    cat(sprintf("plotting: %s\n", fig.fn))
    png(fig.fn, 1000, 1000)
    plot.new()
    plot.window(xlim=c(0,Nc), ylim=c(0,Na))
    rect(sm$i-1, sm$j-1, sm$i, sm$j, col=sm$col, border=NA)
    mtext(text=anchors[o.anchors], side=2, at=1:Na-0.5, las=2, line=0)
#    mtext(text=contigs[o.contigs], side=2, at=1:Nc-0.5, las=2, line=0)
    dev.off()
}
