#######################################
# tree basic functions
#######################################

# convert ids to string for internal purposes of lists
make.id=function(id)
{
    as.character(id)
}

tree.node=function(
    id,
    root.id,
    parent.id, # NULL for root
    depth=1,      # plot depth of node
    ...           # user data
    )
{
    if (sum(id==root.id) != 1) stop(paste("expecting a single root with id:", root.id))
    data.frame(id=make.id(id), is.root=(id==root.id), parent.id=make.id(parent.id), depth=depth, stringsAsFactors=F, ...)
}

make.tree=function(nodes)
{
    nodes$id = make.id(nodes$id)
    nodes$parent.id = make.id(nodes$parent.id)

    if (sum(nodes$is.root) != 1) stop(paste("expecting a single root with id:", root.id))
    root.id = nodes$id[nodes$is.root]

    # compute children
    children = list()
    for (i in 1:dim(nodes)[1]) {
        node = nodes[i,]
        if (node$is.root || node$parent.id < 0)
            next
        if (sum(nodes$id == node$parent.id) == 0)
            stop(paste("parent node not found:", node$parent.id))
        if (sum(nodes$id == node$parent.id) > 1)
            stop(paste("more than 2 parent nodes found:", node$parent.id))
        children[[node$parent.id]] = c(children[[node$parent.id]], node$id)
    }

    # set leaf bit
    for (id in nodes$id)
        nodes$is.leaf[nodes$id == id] = (length(children[[id]])  == 0)

    list(root.id=root.id, nodes=nodes, children=children)
}

# returns subtree
get.subtree=function(tree, id, depth)
{
    result = id
    if (depth > 0) {
        cids = get.children(tree=tree, id=id)
        for (cid in cids)
            result = c(result, get.subtree(tree=tree, id=cid, depth=depth-1))
    }
    result
}

# returns nodes, limited to the spanning tree over a set of ids
get.spanning.tree=function(root.id, nodes, ids)
{
    ids = make.id(ids)

    result.ids = ids
    for (id in ids) {
        tid = id
        while (tid != root.id) {
            tid = nodes$parent.id[nodes$id == tid]
            result.ids = c(result.ids, tid)
        }
    }
    result.ids = c(result.ids, root.id)
    sort(unique(result.ids))
    nodes[is.element(nodes$id, result.ids),]
}

# trim according to function
trim.tree=function(tree, trim.f)
{
    nodes = tree$nodes
    children = tree$children
    root.id = tree$root.id
    get.ids=function(id) {
        node =  nodes[nodes$id == id,]
        if (trim.f(node)) {
            return (NULL)
        }
        result = node$id
        for (child.id in children[[node$id]]) {
            result = c(result, get.ids(child.id))
        }
        result
    }
    ids = get.ids(root.id)
    sort(unique(ids))
}

get.children=function(tree, id) {
    tree$children[[make.id(id)]]
}

# compute map projecting from all ids to a limited set of ids
map.ids=function(nodes, root.id, ids)
{
    result = NULL
    for (id in nodes$id) {
        tid = id
        rid = NULL
        while (tid != root.id) {
            mm = match(tid, ids)
            if (!is.na(mm)) {
                rid = ids[mm]
                break
            }
            tid = nodes$parent.id[nodes$id == tid]
        }
        if (is.null(rid))
            rid = root.id
        result = rbind(result, data.frame(source.id=id, target.id=rid))
    }
    result$source.index = match(result$source.id, nodes$id)
    result$target.index = match(result$target.id, nodes$id)
    result
}

get.lineage.ids=function(tree, id)
{
    nodes = tree$nodes
    root.id = tree$root.id

    result = NULL
    tid = id
    level = 1
    while (tid != root.id) {
        result = rbind(tid, result)
        level = level + 1
        tid = nodes$parent.id[nodes$id == tid]
    }
    c(root.id, result)
}

#######################################
# tree plot functions
#######################################

plot.node.default=function(depth, index, node, leaf, final=F)
{
    points(x=depth, y=index, pch=ifelse(!final,19,1), col=ifelse(leaf,1,2))
    text(x=depth, y=index, pos=4, labels=node$id)
}

# get depth and count, recursive function
get.props=function(tree, id=tree$root.id, trim.f=NULL)
{
    nodes = tree$nodes
    children = tree$children
    node = nodes[nodes$id == id,]
    if (node$is.leaf)
        return (list(count=1, depth=node$depth))
    trimmed = !is.null(trim.f) && trim.f(node)
    result = list(count=0, depth=0)
    if (trimmed)
        return (result)
    children.count = 0
    for (child.id in children[[node$id]]) {
        rchild = get.props(tree=tree, id=child.id, trim.f=trim.f)
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

plot.tree=function(tree, add=F, trim.f=NULL, plot.internal.node.f=plot.node.default, plot.leaf.f=plot.node.default, gap.left=0, gap.right=5)
{
    # first go over tree and compute total depth and number of nodes
    max.depth = 0
    count = 0
    nodes = tree$nodes
    children = tree$children
    root.id = tree$root.id

    props = get.props(tree=tree, id=tree$root.id)
    cat(sprintf("max depth: %d\n", props$depth))
    cat(sprintf("number of leaves: %d\n", props$count))

    # actually do the plot
    plot.node=function(id, depth=0, index=0) {
        node = nodes[nodes$id == id,]
        node.depth = depth+node$depth
        if (node$is.leaf) {
            if (!add) segments(x0=depth, y0=index, x1=node.depth, y1=index)
            plot.leaf.f(depth=node.depth, index=index, node=node, leaf=T)
            return (list(index=index, count=1))
        }
        trimmed = !is.null(trim.f) && trim.f(node)
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
            if (!add) segments(x0=node.depth, y0=children.range[1], x1=node.depth, y1=children.range[2])
        } else {
            node.index = index
            count = 1
        }

        # node stem
        if (!add) segments(x0=depth, y0=node.index, x1=node.depth, y1=node.index)
        plot.internal.node.f(depth=node.depth, index=node.index, leaf=F, node=node, final=(children.count==0))

        return (list(index=node.index, count=count))
    }

    if (!add) {
        plot.new()
        plot.window(xlim=c(-gap.left,props$depth+gap.right), ylim=c(-2,props$count-.5))
    }

    # axis(1)
    # axis(2)
    # grid()
    plot.node(tree$root.id)
    return (0)
}

leaf.order=function(tree, trim.f=NULL)
{
    nodes = tree$nodes
    children = tree$children
    root.id = tree$root.id

    result = NULL
    visit.node=function(id, index=1) {
        node = nodes[nodes$id == id,]
        if (node$is.leaf) {
            result <<- rbind(result, data.frame(index=index, id=id))
            return (1)
        }
        trimmed = !is.null(trim.f) && trim.f(node)
        if (trimmed)
            return (list(index=index, count=0))
        children.count = 0
        for (child.id in children[[node$id]]) {
            rchild = visit.node(id=child.id, index=index)
            children.count = children.count + rchild
            if (rchild > 0) {
                index = index + rchild
            }
        }
        if (children.count > 0) {
            count = children.count
        } else {
            count = 1
        }
        return (count)
    }

    visit.node(tree$root.id)
    return (result)
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
