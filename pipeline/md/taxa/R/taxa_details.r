plot.anchors=function(ifn.tree, ifn.rep, ifn.taxa, min.gene.count=100, fdir)
{
    options(stringsAsFactors=F)

    table = load.table(ifn.tree)
    tax = read.delim(ifn.taxa, quote="")
    rep = load.table(ifn.rep)

    anchors = sort(unique(table$anchor))
    root.id = 1

    # start with all nodes
    nodes = tree.node(id=tax$tax_id, root.id=root.id, parent.id=tax$parent_id, depth=1, name=tax$name, level=tax$level, lineage=tax$lineage)
    cat(sprintf("total number of nodes: %d\n", dim(nodes)[1]))

    # identity
    identity.colors = c("white", "blue")
    identity.breaks = c(50, 100)
    panel.identity = make.color.panel(identity.colors, ncols=256)
    wlegend(fdir=fdir, names=identity.breaks, cols=identity.colors, title="identity")

    # gene count
    genes.colors = c("white", "red")
    genes.breaks = c(50, 100)
    panel.genes = make.color.panel(genes.colors, ncols=256)
    wlegend(fdir=fdir, names=genes.breaks, cols=genes.colors, title="genes")

    pnode=function(depth, index, node, leaf, final=F) {
        label = paste(node$name)
        pch = 15
        pos = if(leaf||final) 4 else 2
        col = ifelse(node$is.rep, "red", "gray")
        lwd = 3
        cex = 2
        x = depth
        y = index

        rect(xleft=x-0.25, xright=x+0.25, ybottom=y-0.4, ytop=y, border=node$identity.col, col=node$identity.col)
        rect(xleft=x-0.25, xright=x+0.25, ybottom=y+0, ytop=y+0.4, border=node$genes.col, col=node$genes.col)
        segments(x0=x-0.25, x1=x+0.25, y0=y, y1=y)

        border.lwd = if (node$is.rep) 4 else 1
        rect(xleft=x-0.25, xright=x+0.25, ybottom=y-0.4, ytop=y+0.4, border=1, lwd=border.lwd)

        if (leaf || final)
            text(x=depth, y=index, labels=label, pos=pos, offset=1)
    }

    for (anchor in anchors) {
        rep.id = make.id(rep$tax_id[rep$anchor == anchor])
        ids = make.id(c(rep.id, table$tax_id[table$anchor == anchor & table$gene_count>min.gene.count]))
        snodes = get.spanning.tree(root.id=root.id, nodes=nodes, ids=ids)
        cat(sprintf("anchor: %d, number of nodes: %d\n", anchor, dim(snodes)[1]))
        table.anchor = table[table$anchor == anchor,]

        total.genes = table.anchor$gene_count[table.anchor$tax_id == root.id]

        ix = match(snodes$id, make.id(table.anchor$tax_id))
        snodes$genes = 100 * table.anchor$gene_count[ix] / total.genes
        snodes$identity = table.anchor$identity[ix]
        snodes$identity.col = panel.identity[vals.to.cols(snodes$identity, identity.breaks)]
        snodes$genes.col = panel.genes[vals.to.cols(snodes$genes, genes.breaks)]

        # mark rep
        snodes$is.rep = (snodes$id == rep.id)
        if (sum(snodes$is.rep) != 1)
            stop("rep.id not in tree")
        tree = make.tree(snodes)

        props = get.props(tree=tree)
        rep.node = snodes[snodes$is.rep,]
        title.text = paste(rep.node$name, ", G=", round(rep.node$genes,1), "%, I=", round(rep.node$identity,1), "%", sep="")

        fig.start(fdir=fdir, ofn=paste(fdir, "/", anchor, ".png", sep=""), height=(200+props$count*25), width=(400+props$depth*25))
        plot.tree(tree=tree, trim.f=NULL, plot.internal.node.f=pnode, plot.leaf.f=pnode, gap.left=10, gap.right=15)
        title(main=paste("#", anchor, ", N=", total.genes, ", ", title.text, sep=""))
        fig.end()
    }

}
