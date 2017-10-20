# converts NCBI taxa to short 3 letter format
short.rank=function(ranks) {
    ranks = sub(" ", ".", ranks)
    rank.list = list(strain="sbS", class="C", family="F", forma="FO", genus="G", infraclass="iC", infraorder="iO", kingdom="K", no.rank="nr", order="O", parvorder="pO", phylum="P", species="S", species.group="Sg" , species.subgroup="Ssg" , subclass="sbC", subfamily="sbF", subgenus="sbG", subkingdom="sbK", suborder="sbO", subphylum="sbP", subspecies="sbS", subtribe="sbT", superclass="spC", superfamily="spF", superkingdom="spK", superorder="spO", tribe="T", varietas="V")
    rank.map = data.frame(rank=names(rank.list), short=unlist(rank.list))
    ix = match(ranks, rank.map$rank)
    ifelse(!is.na(ix), rank.map$short[ix], "X")
}

plot.rep.summary=function(ifn.rep, ifn.legend, ifn.order, fdir)
{
    rep = load.table(ifn.rep)
    table.legend = load.table(ifn.legend)
    rep = merge(rep, table.legend, by="anchor", sort=F)

    anchor.table = load.table(ifn.order)
    anchors = anchor.table$set
    rep = rep[match(anchors, rep$anchor),]

    # rep$short.rank = short.rank(rep$level)
    # rep$divergence = log10(1 + 100 - rep$identity)
    # rep$ugenes = log10(1 + 100 - rep$frac)

    rep$divergence = log10(1 + 100 - rep$anchor.identity)
    rep$ugenes = log10(1 + 100 - rep$anchor.coverage)
    N = dim(rep)[1]

#    rep$label = paste(rep$name, " | ",make.anchor.id(rep$anchor,anchor.table), " | ", rep$short.rank, sep="")
    rep$label = paste(rep$ref.name, " | ",make.anchor.id(rep$anchor,anchor.table), sep="")

    # color levels
    # levels = c("phylum", "family", "genus", "species", "strain", "other")
    # level.colors = c("darkgray", "darkblue", "yellow", "red", "orange", "gray")
    # rep$level = ifelse(is.element(rep$level, levels), rep$level, "other")
    # rep$level.color = level.colors[match(rep$level, levels)]

    fig.start(fdir=fdir, type="pdf", ofn=paste(fdir, "/summary.pdf", sep=""), width=5, height=9)
    layout(matrix(1:3, 1, 3), widths=c(1,1,4))

    axis.coords = c(0, 2, 10, 50)
    axis.at = log10(1 + axis.coords)
    pplot=function(field) {
        par(mai=c(1,0.1,0.3,0.1))
        par(yaxs="i")
        par(xaxs="i")
        plot.new()
        plot.window(xlim=c(0,2.2), ylim=c(0, N))
        axis(1, at=axis.at, labels=axis.coords, las=2)
        abline(v=axis.at, lty=1, col="gray")
        # abline(h=0:N, lty=1, col="gray")
        title(main=field)
        rect(xleft=0, ybottom=1:N-1+0.1, xright=rep[,field], ytop=1:N-0.1, col="darkblue", border=NA)
        box()
    }
    pplot("divergence")
    pplot("ugenes")

    # names, taxa, level
    par(mai=c(1,0,0.3,3.1))
    plot.new()
    plot.window(xlim=c(0,1.5), ylim=c(0, N))
    axis(4, at=1:N-0.5, labels=rep$label, las=2)
    rect(xleft=0, ybottom=1:N-1, xright=1.5, ytop=1:N, col=rep$group.color, border=NA)
    # rect(xleft=1.1, ybottom=1:N-1, xright=2, ytop=1:N, col=rep$level.color, border=NA)
    # rect(xleft=2.1, ybottom=1:N-1, xright=3, ytop=1:N, col=ifelse(rep$is.cag, "black", "gray"), border=NA)
    text(x=0.75, y=1:N-0.5, paste(rep$group.sign, rep$group.index, sep=""))

    fig.end()
}

plot.rep.types=function(ifn.rep, fdir)
{
    rep = load.table(ifn.rep)
    tt = as.matrix(table(rep$level))
    df = data.frame(level=rownames(tt), count=tt)
    rorder = c("phylum", "family", "genus", "species", "strain")
    ix = match(df$level, rorder)
    df$order = ifelse(!is.na(ix), ix, 0)
    df = df[order(df$order),]

    pp=function(clean) {
        fig.start(fdir=fdir, ofn=paste(fdir, "/types", ifelse(clean, "_clean", ""), ".png", sep=""), height=400, width=100+dim(df)[1]*20)
        names = if(!clean) df$level else NA
        mx = barplot(df$count, ylab="count", names.arg=names, border=NA, ylim=c(0, max(df$count*1.1)), las=2)
        if (!clean) text(x=mx, y=df$count, labels=df$count, pos=3)
        fig.end()
    }
    pp(T)
    pp(F)
}

plot.rep.tree=function(ifn.order, ifn.taxa, ifn.rep, ofn.legend,
    trim.levels=NULL, trim.names=NULL, trim.grep.names=NULL, collapse.depth.ids=NULL,
    level.names, level.cols, no.level.color="grey",
    branch.ids, branch.color.indices, no.branch.color="grey",
   fdir)
{
    anchor.table = load.table(ifn.order)
    anchors = anchor.table$set

    cat(sprintf("reading tax table: %s\n", ifn.taxa))
    tax = read.delim(ifn.taxa, quote="")
    tax$level.short = short.rank(tax$level)

    cat(sprintf("reading rep table: %s\n", ifn.rep))
    rep = read.delim(ifn.rep)

    root.id = 1

    # start with all nodes
    nodes = tree.node(id=tax$tax_id, root.id=root.id, parent.id=tax$parent_id, depth=1,
        name=tax$name, level=tax$level, level.short=tax$level.short, lineage=tax$lineage)
    cat(sprintf("total number of nodes: %d\n", dim(nodes)[1]))

    # spanning tree limited to tax_id of reps
    snodes = get.spanning.tree(root.id=root.id, nodes=nodes, ids=rep$tax_id)
    cat(sprintf("rep spanning tree nodes: %d\n", dim(snodes)[1]))

    # add rep counts to nodes
    tt = table(rep$tax_id)
    snodes$rep.count = ifelse(is.na(match(snodes$id, names(tt))), 0, tt[match(snodes$id, names(tt))])
    snodes$index = 1:dim(snodes)[1]
    snodes = snodes[,c(dim(snodes)[2], 1:(dim(snodes)[2]-1))]
    stree = make.tree(snodes)

    # trim out branches
    trim.names = gsub("_", " ", trim.names)
    trim.f=function(node) {
        result = is.element(node$name, trim.names) || is.element(node$level, trim.levels)
        for (p in trim.grep.names)
            result = result || grepl(p, node$name)
        result
    }
    trim.ids = trim.tree(tree=stree, trim.f=trim.f)

    # map from tax_id to selected trimmed tax ids
    map = map.ids(nodes=snodes, root.id=root.id, ids=trim.ids)

    # save full tree
    for (i in 1:dim(snodes)[1])
        snodes$lineage.ids[i] = paste(get.lineage.ids(stree, snodes$id[i]), collapse=" ")

    fig.dir(dir=fdir)

    # save tax tree
    tax = tax[order(tax$tax_id),]
    save.table(tax, paste(fdir, "/taxa.txt", sep=""))

    # fold counts into trimmed tree
    umap = map[map$source.id != map$target.id,]
    for (i in 1:dim(umap)[1])
        snodes$rep.count[umap$target.index[i]] = snodes$rep.count[umap$source.index[i]] + snodes$rep.count[umap$target.index[i]]

    save.table(snodes, paste(fdir, "/tree.txt", sep=""))

    # make trimmed final tree
    fnodes = snodes[is.element(snodes$id, trim.ids),]
    cat(sprintf("final number of nodes: %d\n", dim(fnodes)[1]))

    #####################################################################################
    # plots
    #####################################################################################

    if (length(branch.ids) != length(branch.color.indices))
        stop("set branch ids and colors to be same length")
    branch.cols = colors()[branch.color.indices]

    # branch color legend
    mx = match(branch.ids, tax$tax_id)
    branch.names = ifelse(!is.na(mx), paste(tax$name[mx], tax$level[mx]), "NA")
    wlegend(fdir=fdir, title="branch", cols=branch.cols, names=branch.names)

    # branch color legend
    wlegend(fdir=fdir, title="rank", cols=level.cols, names=level.names)

    # color by branch
    branch.color=function(nodes, root.id) {
        bmap =  map.ids(nodes=nodes, root.id=root.id, ids=branch.ids)
        bmap$match = match(bmap$target.id, branch.ids)
        bmap$col = ifelse(is.na(bmap$match), no.branch.color, branch.cols[bmap$match])
        bmap$col
    }
    fnodes$branch.color = branch.color(nodes=fnodes, root.id=root.id)
    stree$nodes$branch.color = branch.color(nodes=snodes, root.id=root.id)

    # color by level
    level.names = gsub("_", " ", level.names)
    level.color=function(nodes) {
        mm = match(nodes$level, level.names)
        ifelse(!is.na(mm), level.cols[mm], no.level.color)
    }
    fnodes$level.color = level.color(fnodes)
    stree$nodes$level.color = level.color(snodes)

    # assign signs
    ix = fnodes$rep.count > 0
    N = sum(ix)
    fnodes$sign = ""
    fnodes$sign[ix] = letters[1:N]

    # assign to reps the group color and sign
    legend.ids = fnodes$id[ix]
    map2 = map.ids(nodes=snodes, root.id=root.id, ids=legend.ids)
    rep$group.id = map2$target.id[match(rep$tax_id, map2$source.id)]
    rep$group.sign = fnodes$sign[match(rep$group.id, fnodes$id)]
    rep$group.color = fnodes$branch.color[match(rep$group.id, fnodes$id)]

    # assign group index
    ss = split(rep$anchor, rep$group.id)
    df = NULL
    for (i in 1:length(ss)) {
        x = anchors[is.element(anchors,ss[[i]])]
        df = rbind(df, data.frame(anchor=x, index=seq_along(x)))
    }
    rep$group.index = df$index[match(rep$anchor, df$anchor)]
    rep$anchor.id = make.anchor.id(rep$anchor, anchor.table)
    # save reps
    save.table(rep[c("anchor.id", "anchor", "group.id", "group.sign", "group.index", "group.color")], ofn.legend)

    # collapse some node depths
    for (id in collapse.depth.ids)
        fnodes$depth[fnodes$id == id] = 0

    # trimmed tree
    trimmed.tree = make.tree(fnodes)

    pnode.raw=function(depth, index, node, leaf, final=F) {
        hrep = node$rep.count > 0
        label = paste(node$index, ":", node$id, "(", node$rep.count, ")", sep="")

        pch = if (hrep) 15 else 1
        pos = if(leaf||final) 4 else 1
        col = node$level.color
        lwd = 3
        cex = 2
        points(x=depth, y=index, pch=pch, col=col, cex=cex, lwd=lwd)
        if(leaf||final)
            text(x=depth, y=index, labels=label, pos=pos)
        else
            text(x=depth, y=index, labels=label, adj=1.1, srt=45)
    }

    pnode.simple=function(depth, index, node, leaf, final=F) {
        hrep = node$rep.count > 0
#        if (!hrep)
#            return (NULL)
        if (hrep)
            label = paste(node$level.short, " ", node$name, "(", node$rep.count, ")", sep="")
        else
            label = paste(node$level.short, " ", node$name, sep="")

        pch = if (hrep) 15 else 1
        pos = if(leaf||final) 4 else 1
        col = node$branch.color
        lwd = 3
        cex = 4
        points(x=depth, y=index, pch=pch, col=col, cex=cex, lwd=lwd)

        if(leaf||final)
            text(x=depth, y=index, labels=label, pos=pos, offset=1)
        else
            text(x=depth, y=index, labels=label, adj=1.1, srt=45)

        text(x=depth, y=index, labels=node$sign, cex=1.2)
    }

    # only nodes with reps
    pnode.simple2=function(depth, index, node, leaf, final=F) {
        hrep = node$rep.count > 0
        if (!hrep)
            return (NULL)
        label = paste(node$level.short, " ", node$name, "(", node$rep.count, ")", sep="")

        pos = if(leaf||final) 4 else 1
        col = node$branch.color
        lwd = 3
        cex = 4
        rect(xleft=depth-0.25, xright=depth+0.25, ybottom=index-0.4, ytop=index+0.4, col=col, border=NA)

        if(leaf||final)
            text(x=depth, y=index, labels=label, pos=pos, offset=1)
        else
            text(x=depth, y=index, labels=label, adj=1.1, srt=45)

        text(x=depth, y=index, labels=node$sign, cex=1.2)
    }

    plot.tree.wrapper=function(tree, title, height, width, pnode, pnode.title) {
        main = paste(title, "_", pnode.title, sep="")
        fig.start(fdir=fdir, type="pdf", ofn=paste(fdir, "/", main, ".pdf", sep=""), height=height, width=width)
        par(mai=c(2,1,0.5,0.5))
        par(xpd=NA)
        plot.tree(tree=tree, trim.f=NULL, plot.internal.node.f=pnode, plot.leaf.f=pnode)
        title(main=main)
        fig.end()
    }

    plot.tree.wrapper(tree=stree, title="full", height=12, width=14, pnode=pnode.raw, pnode.title="raw")
    plot.tree.wrapper(tree=trimmed.tree, title="trimmed", height=12, width=14, pnode=pnode.raw, pnode.title="raw")

    plot.tree.wrapper(tree=trimmed.tree, title="trimmed", height=10, width=14, pnode=pnode.simple, pnode.title="simple_all")
    plot.tree.wrapper(tree=trimmed.tree, title="trimmed", height=8, width=8, pnode=pnode.simple2, pnode.title="simple_reps")
}
