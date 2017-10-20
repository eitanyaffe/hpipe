short.level=function(levels) {
    levels = sub(" ", ".", levels)
    level.list = list(strain="sbS", class="C", family="F", forma="FO", genus="G", infraclass="iC", infraorder="iO", kingdom="K", no.rank="nr", order="O", parvorder="pO", phylum="P", species="S", species.group="Sg" , species.subgroup="Ssg" , subclass="sbC", subfamily="sbF", subgenus="sbG", subkingdom="sbK", suborder="sbO", subphylum="sbP", subspecies="sbS", subtribe="sbT", superclass="spC", superfamily="spF", superkingdom="spK", superorder="spO", tribe="T", varietas="V")
    level.map = data.frame(level=names(level.list), short=unlist(level.list))
    ix = match(levels, level.map$level)
    ifelse(!is.na(ix), level.map$short[ix], "X")
}

get.props.internal=function(nodes, root.id, ids)
{
    snodes = get.spanning.tree(root.id=root.id, nodes=nodes, ids=ids)
    tree = make.tree(snodes)
    get.props(tree=tree)
}

plot.tree.internal=function(nodes, root.id, ids, table, panel.identity, res.id, rep.id, main, fdir, ofn)
{
    x.gap = 0.2
    y.gap = 0.1
    pnode=function(depth, index, node, leaf, final=F) {
        label = paste(node$name, " I=", round(node$identity,1), " C=", round(node$coverage,1), sep="")
        pos = if(leaf||final) 4 else 2
        x = depth
        y = index
        ybottom = y-0.5+y.gap
        ydelta = 1 - 2*y.gap
        rect(xleft=x-0.5+x.gap, xright=x+0.5-x.gap, ybottom=ybottom, ytop=ybottom+ydelta, col="lightgray", border=NA)
        rect(xleft=x-0.5+x.gap, xright=x+0.5-x.gap, ybottom=ybottom, ytop=ybottom+ydelta*node$coverage/100, col=node$identity.col, border=NA)
        if (node$is.rep)
            rect(xleft=x-0.5+x.gap, xright=x+0.5-x.gap, ybottom=ybottom, ytop=ybottom+ydelta, col=NA, border=1, lty=2)
        if (node$is.res)
            rect(xleft=x-0.5+x.gap, xright=x+0.5-x.gap, ybottom=ybottom, ytop=ybottom+ydelta, col=NA, border=1)
        if (leaf || final)
            text(x=x, y=y, labels=label, pos=pos, offset=0.3, cex=0.8)
    }

    snodes = get.spanning.tree(root.id=root.id, nodes=nodes, ids=ids)

    ix = match(snodes$id, make.id(table$tax_id))
    snodes$parent.level = ifelse(snodes$parent.id > 0, snodes$level[match(snodes$parent.id, snodes$id)], "root")
    snodes$coverage = table$coverage[ix]
    snodes$identity = table$identity[ix]
    snodes$identity.col = table$identity.col[ix]
    snodes$is.res = (snodes$id == res.id)
    snodes$is.rep = (snodes$id == rep.id)
    tree = make.tree(snodes)
    props = get.props(tree=tree)
    plot.tree(tree=tree, plot.internal.node.f=pnode, plot.leaf.f=pnode, gap.left=1, gap.right=25)
    title(main=main)
}

plot.anchor.trees=function(ifn.tree, ifn.rep, ifn.taxa, ifn.res, ifn.summary, min.gene.count, fdir)
{
    library(plotrix)

    table.all = load.table(ifn.tree)
    tax = read.delim(ifn.taxa, quote="")
    rep = load.table(ifn.rep)
    res = load.table(ifn.res)
    summary = load.table(ifn.summary)

    root.id = 1

    # start with all nodes
    nodes = tree.node(id=tax$tax_id, root.id=root.id, parent.id=tax$parent_id, depth=1, name=tax$name, level=tax$level, lineage=tax$lineage)
    nodes$short.level = short.level(nodes$level)
    cat(sprintf("total number of nodes: %d\n", dim(nodes)[1]))

    anchor.ids = rep$anchor.id

    # identity
    identity.colors = c("white", "blue", "red", "orange", "yellow")
    identity.breaks = c(0, 60, 70, 80, 100)
    panel.identity = make.color.panel(identity.colors, ncols=256)
    wlegend(fdir=fdir, names=identity.breaks, cols=identity.colors, title="identity")
    table.all$identity.col = panel.identity[vals.to.cols(table.all$identity, identity.breaks)]

    for (anchor.id in anchor.ids) {
        anchor = rep$anchor[match(anchor.id, rep$anchor.id)]
        rep.id = make.id(rep$tax_id[rep$anchor == anchor])
        res.id = make.id(res$tax.id[res$anchor.id == anchor.id])
        rnode = res[res$anchor.id == anchor.id,]

        table = table.all[table.all$anchor == anchor,]
        gene.total = sum(summary[match(anchor, summary$anchor),-1])
        table$coverage = 100 * table$gene_count / gene.total

        main = paste(anchor.id, ", ", rnode$name, " ",  rnode$short.level, "\nC=",
            round(rnode$coverage,1), "%, I=", round(rnode$identity,1), "%", ", type=", rnode$type, sep="")
        ids = make.id(c(rep.id, res.id, table.all$tax_id[table.all$anchor == anchor & table.all$gene_count>min.gene.count]))

        props = get.props.internal(nodes=nodes, root.id=root.id, ids=ids)
        ofn = paste(fdir, "/", anchor.id, ".pdf", sep="")
        fig.start(fdir=fdir, ofn=ofn,
                  type="pdf", height=(2+props$count*0.2), width=(4+props$depth*0.25))
        plot.tree.internal(nodes=nodes, root.id=root.id, ids=ids, table=table, res.id=res.id, rep.id=rep.id, main=main, fdir=fdir, ofn=ofn)
        fig.end()
    }
}

plot.anchor.trees.limited=function(ifn.tree, ifn.rep, ifn.taxa, ifn.res, ifn.summary, min.gene.count, fdir)
{
    library(plotrix)

    table.all = load.table(ifn.tree)
    tax = read.delim(ifn.taxa, quote="")
    rep = load.table(ifn.rep)
    res = load.table(ifn.res)
    summary = load.table(ifn.summary)

    anchor.ids = rep$anchor.id

    # identity
    identity.colors = c("white", "blue", "red", "orange", "yellow")
    identity.breaks = c(0, 60, 70, 80, 100)
    panel.identity = make.color.panel(identity.colors, ncols=256)
    wlegend(fdir=fdir, names=identity.breaks, cols=identity.colors, title="identity")
    table.all$identity.col = panel.identity[vals.to.cols(table.all$identity, identity.breaks)]

    root.id = 1

    nodes = tree.node(id=tax$tax_id, root.id=root.id, parent.id=tax$parent_id, depth=1, name=tax$name, level=tax$level, lineage=tax$lineage)
    nodes$short.level = short.level(nodes$level)
    cat(sprintf("total number of nodes: %d\n", dim(nodes)[1]))

    for (anchor.id in anchor.ids) {
        anchor = rep$anchor[match(anchor.id, rep$anchor.id)]
        rnode = res[res$anchor.id == anchor.id,]
        table = table.all[table.all$anchor == anchor,]
        gene.total = sum(summary[match(anchor, summary$anchor),-1])
        table$coverage = 100 * table$gene_count / gene.total

        if (rnode$type == "ok")
            next
        rep.id = make.id(rep$tax_id[rep$anchor == anchor])
        res.id = rnode$tax.id

        tree = make.tree(nodes=nodes[is.element(nodes$id, table$tax_id),])
        ids = get.subtree(tree=tree, id=res.id, depth=rnode$depth)
        ids = table$tax_id[is.element(table$tax_id, ids) & table$gene_count>min.gene.count]

        main = paste(anchor.id, ", ", rnode$name, " ",  rnode$short.level, "\nC=",
            round(rnode$coverage,1), "%, I=", round(rnode$identity,1), "%", ", type=", rnode$type, sep="")

        ofn = paste(fdir, "/", anchor.id, "_", rnode$type, ".pdf", sep="")
        props = get.props.internal(nodes=nodes, root.id=root.id, ids=ids)
        fig.start(fdir=fdir, ofn=ofn,
                  type="pdf", height=(2+props$count*0.2), width=(4+props$depth*0.25))
        plot.tree.internal(nodes=nodes, root.id=root.id, ids=ids, table=table, res.id=res.id, rep.id=rep.id, main=main, fdir=fdir, ofn=ofn)
        fig.end()
    }

}
