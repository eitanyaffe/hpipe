#######################################
# tree basic functions
#######################################

tree.node=function(
    id,
    parent.id=NA, # NULL for root
    depth=1,      # plot depth of node
    ...           # user data
    )
{
    data.frame(id=id, is.root=is.na(parent.id), parent.id=parent.id, depth=depth, stringsAsFactors=F, ...)
}

make.tree=function(root.id, ...)
{
    nodes = NULL
    for (node in list(...))
        nodes = rbind(nodes, node)

    # compute children
    children = list()
    for (i in 1:dim(nodes)[1]) {
        node = nodes[i,]
        if (!node$is.root)
            children[[node$parent.id]] = c(children[[node$parent.id]], node$id)
    }

    # set leaf bit
    for (id in nodes$id)
        nodes$is.leaf[nodes$id == id] = (length(children[[id]])  == 0)

    list(root.id=root.id, nodes=nodes, children=children)
}

#######################################
# tree plot functions
#######################################

plot.node.default=function(depth, index, node, leaf, final=F)
{
    points(x=depth, y=index, pch=ifelse(!final,19,1), col=ifelse(leaf,1,2))
    text(x=depth, y=index, pos=4, labels=node$id)
}

plot.tree=function(tree, trim.f=NULL, plot.internal.node.f=plot.node.default, plot.leaf.f=plot.node.default)
{
    # first go over tree and compute total depth and number of nodes
    max.depth = 0
    count = 0
    nodes = tree$nodes
    children = tree$children
    root.id = tree$root.id

    # get depth and count
    get.props=function(id=root.id) {
        node = nodes[nodes$id == id,]
        if (node$is.leaf)
            return (list(count=1, depth=node$depth))
        trimmed = !is.null(trim.f) && !trim.f(node)
        result = list(count=0, depth=0)
        if (trimmed)
            return (result)
        children.count = 0
        for (child.id in children[[node$id]]) {
            rchild = get.props(child.id)
            children.count = children.count + rchild$count
            if (rchild$count > 0) {
                result$count = result$count + rchild$count
                result$depth = max(result$depth, node$depth+rchild$depth)
            }
        }
        if (children.count == 0) {
            result = list(count=1, depth=node$depth)
        } else {
            return (result)
        }
    }
    props = get.props(tree$root.id)
    cat(sprintf("max depth: %d\n", props$depth))
    cat(sprintf("number of leaves: %d\n", props$count))

    # actually do the plot
    plot.node=function(id, depth=0, index=0) {
        node = nodes[nodes$id == id,]
        node.depth = depth+node$depth
        if (node$is.leaf) {
            segments(x0=depth, y0=index, x1=node.depth, y1=index)
            plot.leaf.f(depth=node.depth, index=index, node=node, leaf=T)
            return (list(index=index, count=1))
        }
        trimmed = !is.null(trim.f) && !trim.f(node)
        if (trimmed)
            return (list(index=index, count=0))
        children.range = NULL
        children.count = 0
        for (child.id in children[[node$id]]) {
            rchild = plot.node(id=child.id, depth=node.depth, index=index)
            children.count = children.count + rchild$count
            if (rchild$count > 0) {
                index = index + rchild$count
                children.range = range(c(children.range, rchild$index))
            }
        }
        if (children.count > 0) {
            node.index = mean(children.range)
            count = children.count

            # connect children stems
            segments(x0=node.depth, y0=children.range[1], x1=node.depth, y1=children.range[2])
        } else {
            node.index = index
            count = 1
        }

        # node stem
        segments(x0=depth, y0=node.index, x1=node.depth, y1=node.index)
        plot.internal.node.f(depth=node.depth, index=node.index, leaf=F, node=node, final=(children.count==0))

        return (list(index=node.index, count=count))
    }

    plot.new()
    plot.window(xlim=c(0,props$depth+15), ylim=c(-2,props$count-.5))
    axis(1)
    axis(2)
    grid()
    plot.node(tree$root.id)
    return (0)
}

#######################################
# test functions
#######################################

test.tree=function()
{
    tree = make.tree(
        root.id = "n1",
        tree.node("n1"),
        tree.node("n2", "n1"),
        tree.node("n3", "n1"),
        tree.node("n4", "n3"),
        tree.node("n5", "n3"))

    plot.tree(tree)
}

test2=function()
{
    xx = read.delim("/relman02/users/eitany/bcc_output/Relman/assembly/pre_minia_contig/fold_0/datasets/pre_hic/anchors/pre_hic/uniref/uniref100_2015_12/taxa_trees")
    xt = read.delim("/relman02/users/eitany/bcc_output/Relman/assembly/pre_minia_contig/fold_0/datasets/pre_hic/anchors/pre_hic/uniref/uniref100_2015_12/taxa_table")
    mm = match(xx$tax_id, xt$tax_id)
    for (field in c("parent_id", "name", "level"))
        xx[,field]= xt[mm,field]

    anchors = sort(unique(xx$anchor))

    ptrim = function(node) {
        node$gene.count > 200
    }
    pnode=function(depth, index, node, leaf, final=F) {
        labels = if(leaf||final) paste(node$name, node$gene.count, sep="_", collapse="") else node$name
        srt = ifelse(leaf||final, 0, 45)
        pos = if(leaf||final)  4 else NULL
        adj = c(1,1.2)

        points(x=depth, y=index, pch=ifelse(!final,19,1), col=ifelse(leaf,1,2))
        text(x=depth, y=index, pos=pos, adj=adj, labels=labels, srt=srt)
    }

    xx$parent_id = ifelse(xx$parent_id != -1, paste("t", xx$parent_id, sep=""), NA)
    xx$tax_id = paste("t", xx$tax_id, sep="")

    for (anchor in anchors) {
        xxa = xx[xx$anchor == anchor,]
        tree = make.tree(root.id="t1", tree.node(id=xxa$tax_id, parent.id=xxa$parent_id, depth=1, name=xxa$name, gene.count=xxa$gene_count))
        plot.tree(tree=tree, trim.f=ptrim, plot.internal.node.f=pnode, plot.leaf.f=pnode)

    }
}
