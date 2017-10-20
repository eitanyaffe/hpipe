make.legend.tree=function(taxa, df, default.value)
{
    root.id = taxa$tax_id[taxa$name == "Bacteria"]
    all.nodes = tree.node(id=taxa$tax_id, root.id=root.id, parent.id=taxa$parent_id)

    ids = df$tax_id[!is.na(df$tax_id)]
    nodes = get.spanning.tree(root.id=root.id, nodes=all.nodes, ids=ids)
    nodes$is.legend = is.element(nodes$id, ids)
    ix = match(nodes$id, df$tax_id)
    nodes$value = ifelse(!is.na(ix), df$value[ix], default.value)

    names = taxa$name[match(nodes$id, taxa$tax_id)]
    nodes$name = ifelse(!is.na(ix), paste(names, " (", df$count[ix], ")", sep=""), names)
    make.tree(nodes)
}

pnode.text=function(depth, index, node, leaf, final=F)
{
    label = node$name
    pos = if(leaf||final) 4 else 1
    if(leaf||final)
        text(x=depth, y=index, labels=label, pos=pos, offset=1)
    else
        text(x=depth, y=index, labels=label, adj=1, srt=45)
}

plot.tree.wrapper=function(fdir, tree, title, height, width, pnode, pnode.text) {
    fig.start(fdir=fdir, type="pdf", ofn=paste(fdir, "/", title, ".pdf", sep=""), height=height, width=width)
    par(mai=c(1,0,0.5,0.5))
    par(xpd=NA)
    plot.tree(tree=tree, trim.f=NULL, plot.internal.node.f=pnode, plot.leaf.f=pnode, gap.right=1)
    plot.tree(tree=tree, trim.f=NULL, plot.internal.node.f=pnode.text, plot.leaf.f=pnode.text, add=T)
    title(main=title)
    fig.end()
}

plot.taxa.legend.color=function(ifn.taxa, ifn.legend.color, fdir)
{
    taxa = load.table(ifn.taxa)
    df = load.table(ifn.legend.color)
    tree = make.legend.tree(taxa=taxa, df=df, default.value="lightgray")

    pnode=function(depth, index, node, leaf, final=F) {
        rect(xleft=depth-0.25, xright=depth+0.1, ybottom=index-0.4, ytop=index+0.4, col=node$value, border=NA)
    }
    plot.tree.wrapper(fdir=fdir, tree=tree, pnode=pnode, pnode.text=pnode.text, title="color", height=6, width=10)

    lo = leaf.order(tree)
    df$index = lo$index[match(df$tax_id, lo$id)]
    df$index[is.na(df$name)] = 0
    df = df[order(df$index, decreasing=T),]
    df$name[is.na(df$name)] = "unknown"
    names = paste(df$name, " (", df$count, ")", sep="")
    wlegend(fdir, names=names, cols=df$value, title="color_short")

}

plot.taxa.legend.letter=function(ifn.taxa, ifn.legend.color, ifn.legend.letter, fdir)
{
    taxa = load.table(ifn.taxa)
    df = load.table(ifn.legend.letter)
    df.color = load.table(ifn.legend.color)
    tree = make.legend.tree(taxa=taxa, df=df, default.value="")

    # add color
    root.id = taxa$tax_id[taxa$name == "Bacteria"]
    ix = match(tree$nodes$id, df.color$tax_id)
    tree$nodes$color = ifelse(!is.na(ix), df.color$value[ix], "lightgray")
    for (id in tree$nodes$id) {
        tid = id
        color = "lightgray"
        while (tid != root.id && color == "lightgray") {
            node = tree$nodes[tree$nodes$id == tid,]
            if (node$color != "lightgray")
                color = node$color
            tid = node$parent.id
        }
        tree$nodes$color[tree$nodes$id == id] = color
    }

    pnode=function(depth, index, node, leaf, final=F) {
        rect(xleft=depth-0.1, xright=depth+0.1, ybottom=index-0.4, ytop=index+0.4, col=node$color, border=NA)
        if (node$is.legend)
            text(x=depth, y=index, labels=node$value)
    }
    plot.tree.wrapper(fdir=fdir, tree=tree, pnode=pnode, pnode.text=pnode.text, title="letter", height=8, width=10)

    lo = leaf.order(tree)
    df$index = lo$index[match(df$tax_id, lo$id)]
    df$index[is.na(df$tax_id)] = dim(df)[1]
    df = df[order(df$index, decreasing=T),]
    df$name[is.na(df$tax_id)] = "unknown"

    ix = match(df$tax_id, tree$nodes$id)
    df$color = ifelse(!is.na(ix), tree$nodes$color[ix], "lightgray")

    names = paste(df$name, " (", df$count, ")", sep="")
    letters = df$value
    colors = df$color

    N = length(names)
    fig.start(fdir=fdir, width=4, height=0.5 + N*0.18, ofn=paste(fdir, "/letter_legend.pdf", sep=""), type="pdf")
    par(mai=c(0,0,0,3.7))
    par(xpd=NA)
    plot.new()
    plot.window(xlim=c(0,1), ylim=c(0,N))
    rect(xleft=rep(0.25, N), xright=rep(0.75, N), ybottom=1:N-0.9, ytop=1:N-0.1, col=colors, border=NA)
    text(x=rep(0.5,N), y=1:N-0.5, labels=letters)
    text(x=rep(0.7,N), y=1:N-0.5, labels=names, pos=4)
    fig.end()
}

