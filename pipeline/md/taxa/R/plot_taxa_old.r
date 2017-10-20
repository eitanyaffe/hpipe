
    ## ranks = c("Superkingdom", "Kingdom", "Subkingdom",
    ##     "Superphylum", "Phylum", "Subphylum",
    ##     "Superclass", "Class", "Subclass", "Infraclass",
    ##     "Superorder", "Order", "Suborder", "Infraorder", "Parvorder",
    ##     "Superfamily", "Family", "Subfamily",
    ##     "Tribe", "Subtribe",
    ##     "Genus", "Subgenus",
    ##     "Species group", "Species subgroup", "Species", "Subspecies",
    ##     "Varietas", "Forma")

add.uniref=function(ifn, ofn, uniref.ifn)
{
  options(stringsAsFactors=F)

  cat(sprintf("reading table: %s\n", ifn))
  table = read.delim(ifn)

  cat(sprintf("reading table: %s\n", uniref.ifn))
  uniref = read.delim(uniref.ifn)

  m = merge(table, uniref)

  cat(sprintf("input: %d, uniref hits: %d\n", dim(table)[1], dim(m)[1]))
  cat(sprintf("saving result table: %s\n", ofn))
  write.table(m, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}

get.children.list=function(tax)
{
    children.list = list()
    for (i in 1:dim(tax)[1]) {
        parent = tax$parent_id[i]
        child = tax$tax_id[i]
        if (is.na(parent))
            next
        children.list[[as.character(parent)]] = c(children.list[[as.character(parent)]], child)
    }
    return (children.list)
}

plot.tree=function(tree, table, tax, main, ofn, highlight.id=-1, color.tip=T, adj.node=0.5, tip.label.f, node.label.f, rep, step=16)
{
    tip.ids = as.numeric(tree$tip.label)
    node.ids = as.numeric(tree$node.label)
    node.ids[1] = 0

    tree$tip.label = tip.label.f(ids=tip.ids, table=table)
    tree$node.label = node.label.f(ids=node.ids, table=table)

    # cat(sprintf("plotting: %s\n", ofn))
    png(ofn, width=1400, height=(150 + length(tip.ids) * step))
    par(mai=c(0.4, 0.2, 0.5, 0.2))

    v = plot(tree, main=main, cex=1, plot=F)
    xlim = v$x.lim
    ylim = v$y.lim
    dd = diff(xlim)
    xlim[1] = xlim[1] - 0.1*dd
    xlim[2] = xlim[2] + 0.1*dd

    plot(tree, main=main, cex=1, x.lim=xlim, y.lim=ylim, show.tip.label=F)
    rep.col = ifelse(tip.ids == highlight.id, "orange", "white")
    tiplabels(text=tree$tip.label, cex=1, bg=rep.col, adj=0)

    color.tip.f=function(ids) {
        cols = colorpanel(101, low="white", high="red")
        v = table[match(ids, table$tax_id),"frac"]
        idx = log10(pmax(1,100*v))*50+1
        idx[idx==0] = 1
        cols[idx]
    }

    if (color.tip) {
        cols = color.tip.f(tip.ids)
        tiplabels(text="     ", adj=1.5, bg=cols)
    }

    node.col = ifelse(node.ids == highlight.id, "orange", "white")
    nodelabels(text=tree$node.label, adj=adj.node, cex=0.9, bg=node.col)

    # tiplabels(text=tree$tip.label, adj=0, bg=NA)
    dev.off()
}

plot.taxa=function(ifn, ifn.taxa, ifn.summary, ifn.rep, fdir, skip.ids, select.ratio, ifn.order)
{
    options(stringsAsFactors=F)
    library(ape)
    library(gplots)

    cat(sprintf("reading table: %s\n", ifn))
    table = read.delim(ifn)

    cat(sprintf("reading rep table: %s\n", ifn.rep))
    rep = read.delim(ifn.rep)

    cat(sprintf("reading tax table: %s\n", ifn.taxa))
    tax = read.delim(ifn.taxa, quote="")

    anchors = unique(table$anchor)

    # set children list
    children.list = get.children.list(tax)
    base.root = 1

    cat(sprintf("reading summary table: %s\n", ifn.summary))
    summary = read.delim(ifn.summary, quote="")

    # use the total gene count for collapsing
    collapse = "lca"
    table$collapse.count = switch(collapse,
        total=table$gene_count,
        lca=table$lca_count,
        weight=table$weight)

    # construct tree
    build.tree=function(tax.id, table, min.count, depth=2) {
        count = ifelse(is.na(match(tax.id,table$tax_id)), 0, table$collapse.count[match(tax.id,table$tax_id)])
        if (count < min.count)
            return ("")

        if (!(as.character(tax.id) %in% names(children.list)))
            return (sprintf("%d:%d", tax.id, depth))
        children = children.list[[as.character(tax.id)]]
        child.result = NULL
        child.count = NULL
        for (i in seq_along(children)) {
            child.result[i] = build.tree(tax.id=children[i], table=table, min.count=min.count, depth=depth+2)
            child.count[i] = ifelse(is.na(match(children[i],table$tax_id)), 0, table$collapse.count[match(children[i],table$tax_id)])
        }

        # children with low coverage, stop here
        if (all(child.count < min.count))
            return (sprintf("%d:%d", tax.id, depth))

        child.result = child.result[gsub("[\\(|\\)|,]", "", child.result) != ""]
        if (length(child.result) > 1)
            return (sprintf("(%s)%s:%d", paste(child.result, sep=",", collapse=","), tax.id, depth))
        if (length(child.result) == 1)
            return (sprintf("%s", paste(child.result, sep=",", collapse=",")))
        if (length(child.result) == 0)
            return ("")
    }

    fdir.anchors = paste(fdir, "/by_anchors", sep="")
    system(paste("mkdir -p", fdir, fdir.anchors))

    df.main = data.frame(
        anchor=summary$anchor, strong=summary$hit, weak=summary$weak, masked=summary$masked,
        no_taxa=summary$no_taxa, no_hit=summary$no_hit, shared=summary$shared)

    cat(sprintf("reading order table: %s\n", ifn.order))
    otable = read.delim(ifn.order)

    df = df.main[match(otable$set, df.main$anchor),]
    scol = c("red", "darkred", "blue", "darkgrey", "black", "orange")

    # genes
    ofn = paste(fdir, "/gene_summary.png", sep="")
    cat(sprintf("plotting: %s\n", ofn))
    png(ofn, width=1200, height=400)
    barplot(t(as.matrix(df[,-1])), col=scol, ylab="#genes", xlab="anchors", las=2, border=NA, names.arg=df$anchor)
    dev.off()
    wlegend(fdir=fdir, title="gene_summary", cols=scol, names=names(df)[-1])

    # fraction of genes
    ofn = paste(fdir, "/gene_percent_summary.png", sep="")
    cat(sprintf("plotting: %s\n", ofn))
    png(ofn, width=1200, height=400)
    barplot(t(as.matrix(100*df[,-1]/rowSums(df[,-1]))), col=scol, ylab="%genes", xlab="anchors", las=2, border=NA, names.arg=df$anchor)
    dev.off()

    table$total = rowSums(summary[match(table$anchor, summary$anchor),-1])
    table$frac = table$weight / table$total

    tip.label.f = function(ids, table) {
        tab = table[match(ids, table$tax_id),]
        tax = tax[match(ids, tax$tax_id),]
        paste(tab$lca_count, "/", round(tab$weight), "/", tab$gene_count, " | ",
              "I=", round(tab$identity,2), "%", " | ",
              "C=", round(100*tab$coverage,1), "%", " | ",
              "S=", round(100*tab$shared,1), "%", " | ",
              tax$name, " | ",
              ids, sep="")
    }
    node.label.f = function(ids, table) {
        tab = table[match(ids, table$tax_id),]
        tax = tax[match(ids, tax$tax_id),]
        paste(tab$lca_count, "/", round(tab$weight), "/", tab$gene_count, " | ", tax$name, sep="")
    }

    cat(sprintf("plotting anchors in: %s\n", fdir.anchors))
    for (anchor in anchors) {
        ixx = match(anchor, summary$anchor)
        shared = summary$shared[ixx]
        no.hit = summary$no_hit[ixx]
        no.taxa = summary$no_taxa[ixx]
        mask = summary$masked[ixx]
        weak = summary$weak_hit[ixx]
        strong = summary$hit[ixx]

        ttable = table[table$anchor == anchor,]
        cat(sprintf("plotting anchor: %s\n", anchor))

        rep.id = rep$tax_id[match(anchor, rep$anchor)]
        found = F
        min.count = 11
        while (!found) {
            root.result = NULL

            root.result = build.tree(tax.id=base.root, table=ttable, min.count=min.count)
            str = sprintf("%s;", root.result)

            if (gsub("[\\(|\\)|,]", "", str) != str) {
                found = T
            } else {
                cat(sprintf("single node, decreasing trimming\n"))
                min.count = min.count - 2
            }
            if (min.count < 0) {
                cat(sprintf("single node, skipping anchor\n"))
                break
            }
        }
        if (!found)
            next

        tree = read.tree(text=str)
        main = paste("anchor:", anchor, ", shared=", shared, ", no_hit=", no.hit, ", no_taxa=", no.taxa,
            ", masked=", mask, ", weak_hit=", weak, ", strong_hit=", strong, sep="")
        ofn = paste(fdir.anchors, "/", anchor, ".png", sep="")

        plot.tree(tree=tree, table=ttable, tax=tax, main=main, ofn=ofn, highlight.id=rep.id, tip.label.f=tip.label.f, node.label.f=node.label.f)
    }
}

plot.taxa.all=function(ifn, ifn.taxa, ifn.summary, fdir, skip.ids, select.ratio)
{
    options(stringsAsFactors=F)
    library(ape)
    library(gplots)

    cat(sprintf("reading table: %s\n", ifn))
    table = read.delim(ifn)

    cat(sprintf("reading tax table: %s\n", ifn.taxa))
    tax = read.delim(ifn.taxa, quote="")

    # set children list
    children.list = get.children.list(tax)
    base.root = 1

    cat(sprintf("reading summary table: %s\n", ifn.summary))
    summary = read.delim(ifn.summary, quote="")

    # use the total gene count for collapsing
    collapse = "lca"
    table$collapse.count = switch(collapse,
        total=table$gene_count,
        lca=table$lca_count,
        weight=table$weight)

    # construct tree
    build.tree=function(tax.id, table, min.count, depth=2) {
        count = ifelse(is.na(match(tax.id,table$tax_id)), 0, table$collapse.count[match(tax.id,table$tax_id)])
        if (count < min.count)
            return ("")

        if (!(as.character(tax.id) %in% names(children.list)))
            return (sprintf("%d:%d", tax.id, depth))
        children = children.list[[as.character(tax.id)]]
        child.result = NULL
        child.count = NULL
        for (i in seq_along(children)) {
            child.result[i] = build.tree(tax.id=children[i], table=table, min.count=min.count, depth=depth+2)
            child.count[i] = ifelse(is.na(match(children[i],table$tax_id)), 0, table$collapse.count[match(children[i],table$tax_id)])
        }

        # children with low coverage, stop here
        if (all(child.count < min.count))
            return (sprintf("%d:%d", tax.id, depth))

        child.result = child.result[gsub("[\\(|\\)|,]", "", child.result) != ""]
        if (length(child.result) > 1)
            return (sprintf("(%s)%s:%d", paste(child.result, sep=",", collapse=","), tax.id, depth))
        if (length(child.result) == 1)
            return (sprintf("%s", paste(child.result, sep=",", collapse=",")))
        if (length(child.result) == 0)
            return ("")
    }

    df = data.frame(strong=summary$hit, weak=summary$weak, masked=summary$masked, no_taxa=summary$no_taxa,  no_hit=summary$no_hit)
    scol = c("red", "orange", "blue", "darkgrey", "black")

    # genes
    ofn = paste(fdir, "/all_gene_summary.png", sep="")
    cat(sprintf("plotting: %s\n", ofn))
    png(ofn, width=400, height=400)
    barplot(unlist(df), col=scol, ylab="#genes", xlab="", las=2, border=NA, names.arg=names(df))
    legend("topright", fill=scol, legend=names(df)[-1])
    dev.off()

    tip.label.f = function(ids, table) {
        tab = table[match(ids, table$tax_id),]
        tax = tax[match(ids, tax$tax_id),]
        paste(tab$lca_count, "/", round(tab$weight), "/", tab$gene_count, " | ",
              "I=", round(tab$identity,2), "%", " | ",
              "C=", round(100*tab$coverage,1), "%", " | ",
              "S=", round(100*tab$shared,1), "%", " | ",
              tax$name, " | ",
              ids, sep="")
    }
    node.label.f = function(ids, table) {
        tab = table[match(ids, table$tax_id),]
        tax = tax[match(ids, tax$tax_id),]
        paste(tab$lca_count, "/", round(tab$weight), "/", tab$gene_count, " | ", tax$name, sep="")
    }

    no.hit = summary$no_hit
    no.taxa = summary$no_taxa
    mask = summary$masked
    weak = summary$weak_hit
    strong = summary$hit

    found = F
    for (min.count in c(200, 500, 1000, 1500)) {
        root.result = build.tree(tax.id=base.root, table=table, min.count=min.count)
        str = sprintf("%s;", root.result)
        tree = read.tree(text=str)
        main = paste("no_hit=", no.hit, ", no_taxa=", no.taxa, ", masked=", mask, ", weak_hit=", weak, ", strong_hit=", strong, sep="")
        ofn = sprintf("%s/all_tree_%d.png", fdir, min.count)
        cat(sprintf("plotting tree: %s\n", ofn))
        plot.tree(tree=tree, table=table, tax=tax, main=main, ofn=ofn, tip.label.f=tip.label.f, node.label.f=node.label.f)
    }
}

gpath=function(tax, id=1262964)
{
    if (id == 1)
        return (id)
    else
        return (paste(gpath(tax, tax$parent_id[match(id, tax$tax_id)]), "_", id, sep=""))
}

plot.reps=function(ifn, ifn.taxa, ifn.rep, ifn.summary, fdir)
{
    options(stringsAsFactors=F)
    library(ape)
    library(gplots)

    cat(sprintf("reading table: %s\n", ifn))
    table = read.delim(ifn)

    cat(sprintf("reading summary table: %s\n", ifn.summary))
    summary = read.delim(ifn.summary, quote="")

    table$total = rowSums(summary[match(table$anchor, summary$anchor),-1])
    table$frac = table$weight / table$total

    cat(sprintf("reading tax table: %s\n", ifn.taxa))
    tax = read.delim(ifn.taxa, quote="")

    cat(sprintf("reading rep table: %s\n", ifn.rep))
    rep = read.delim(ifn.rep)

    # set children list
    children.list = get.children.list(tax)
    base.root = 1

    # construct tree
    build.rep.tree=function(tax.id, selected.ids, depth=2) {
        if (!(as.character(tax.id) %in% names(children.list)))
            return (ifelse(is.element(tax.id, selected.ids), sprintf("%d:%d", tax.id, depth), ""))

        children = children.list[[as.character(tax.id)]]
        child.result = NULL
        for (i in seq_along(children))
            child.result[i] = build.rep.tree(children[i], selected.ids=selected.ids, depth=depth+2)

        child.result = child.result[child.result != ""]
        if (length(child.result) > 1)
            return (sprintf("(%s)%s:%d", paste(child.result, sep=",", collapse=","), tax.id, depth))
        if (length(child.result) == 1)
            return (sprintf("%s", paste(child.result, sep=",", collapse=",")))
        if (length(child.result) == 0)
            return (ifelse(is.element(tax.id, selected.ids), sprintf("%d:%d", tax.id, depth), ""))
    }

    system(paste("mkdir -p", fdir))

    label.f = function(ids, table) {
        name = rep[match(ids, rep$tax_id),"name"]
        anchor = sapply(split(rep$anchor, rep$tax_id), function(x) paste(x, collapse=","))
        anchor = anchor[match(ids, as.numeric(names(anchor)))]
        sname = tax$name[match(ids, tax$tax_id)]
        ifelse(!is.na(name), paste(anchor, name, sep=" | "), sname)
    }

    ids = unique(rep$tax_id)
    all.tree = read.tree(text=paste(build.rep.tree(tax.id=base.root, selected.ids=ids), ";", sep=""))
    main = ""
    ofn = paste(fdir, "/rep_tree.png", sep="")
    cat(sprintf("plotting rep tree: %s\n", ofn))

    plot.tree(tree=all.tree, table=table, tax=tax, rep=rep, main=main, ofn=ofn, adj.node=0.5, col=F,
              tip.label.f=label.f, node.label.f=label.f, step=16)

    # plot rep species details

    ofn = paste(fdir, "/rep_details.png", sep="")
    cat(sprintf("plotting: %s\n", ofn))
    png(ofn, width=800, height=1200)
    par(mai=c(0.5, 7, 0.5, 0.2))
    rcols = c("red", "orange", "blue")
    rep$legend = paste(rep$name, rep$tax_id, rep$anchor, sep=" | ")

    rep = rep[order(rep$frac),]
    v = barplot(t(as.matrix(rep[,c("frac", "identity", "coverage")])), beside=T, col=rcols, ylab="%", las=2, border=NA, plot=F, horiz=T)
    xlim = c(0.5, max(v)*1.1)
    barplot(t(as.matrix(rep[,c("frac", "identity", "coverage")])), beside=T, col=rcols, ylab="", xlab="%", las=2, border=NA, names.arg=rep$legend, ylim=xlim, horiz=T)
    legend("topright", fill=rcols, legend=c("fraction", "identity", "coverage"))
    dev.off()
}
