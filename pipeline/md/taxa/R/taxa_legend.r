taxa.levels=function(ifn.species, ifn.path, levels, ofn.level.prefix, ofn.level.summary.prefix)
{
    df.species = load.table(ifn.species)
    df.path = load.table(ifn.path)
    ids = df.species$anchor.id
    for (level in levels) {
        # all anchors
        ofn = paste(ofn.level.prefix, ".", level, sep="")
        tlevel = df.path[df.path$level == level,]
        df.level = data.frame(anchor.id=ids)
        df.level$anchor = lookup(table=df.level, lookup.table=df.species, lookup.field="anchor.id", value.field="anchor", na.value="NA")
        df.level$is.cag = lookup(table=df.level, lookup.table=df.species, lookup.field="anchor.id", value.field="is.cag", na.value="NA")
        df.level$tax_id = lookup(table=df.level, lookup.table=tlevel, lookup.field="anchor.id", value.field="tax_id", na.value="NA")
        df.level$name = lookup(table=df.level, lookup.table=tlevel, lookup.field="anchor.id", value.field="name", na.value="NA")
        save.table(df.level, ofn)

        # level summary
        ofn = paste(ofn.level.summary.prefix, ".", level, sep="")
        tt = table(df.level$tax_id)
        tt = tt[order(tt, decreasing=T)]
        ix = match("NA", names(tt))
        if (!is.na(ix))
            tt = c(tt[-ix], tt[ix])
        df.summary = data.frame(tax_id=names(tt))
        df.summary$name = lookup(table=df.summary, lookup.table=df.level, lookup.field="tax_id", value.field="name", na.value="NA")
        df.summary$count = tt
        save.table(df.summary, ofn)
    }
}

taxa.legend=function(
    ifn.taxa, ifn.species, ifn.level.prefix, ifn.level.summary.prefix,
    letter.level, letter.order,
    color.level, force.ids, force.colors,
    ofn.legend.letter, ofn.legend.color, ofn)
{
    tax = read.delim(ifn.taxa, quote="")
    result = load.table(ifn.species)
    result = result[,c("anchor.id", "anchor", "is.cag")]

    # letter
    level = letter.level
    df = load.table(paste(ifn.level.prefix, ".", level, sep=""))
    df.letter = load.table(paste(ifn.level.summary.prefix, ".", level, sep=""))

    # place into tree to have letters in order
    all.nodes = tree.node(id=tax$tax_id, root.id=1, parent.id=tax$parent_id)
    ids = df.letter$tax_id[!is.na(df.letter$tax_id)]
    nodes = get.spanning.tree(root.id=1, nodes=all.nodes, ids=ids)
    tree=make.tree(nodes)
    lo = leaf.order(tree)

    if (letter.order == "by.tree") {
        # by tree order
        ix = match(df.letter$tax_id, lo$id)
        df.letter$value = ifelse(!is.na(ix), letters[lo$index[ix]], "-")
    } else {
        # by number of anchors
        df.letter$value = letters[1:dim(df.letter)[1]]
        df.letter$value[is.na(df.letter$name)] = "-"
    }

    df$letter = df.letter$value[match(df$tax_id, df.letter$tax_id)]
    result$letter = df$letter[match(result$anchor.id, df$anchor.id)]
    result$sub.group.id = df$tax_id[match(result$anchor.id, df$anchor.id)]

    # color
    level = color.level
    df = load.table(paste(ifn.level.prefix, ".", level, sep=""))
    df.color = load.table(paste(ifn.level.summary.prefix, ".", level, sep=""))
    df.color$value = rainbow(dim(df.color)[1])

    # force colors
    for (i in seq_along(force.ids)) {
        id = force.ids[i]
        color.index = force.colors[i]
        df.color$value[match(id, df.color$tax_id)] = colors()[color.index]
    }

    df.color$value[is.na(df.color$name)] = "gray"
    df$color = df.color$value[match(df$tax_id, df.color$tax_id)]
    result$color = df$color[match(result$anchor.id, df$anchor.id)]
    result$group.id = df$tax_id[match(result$anchor.id, df$anchor.id)]

    save.table(df.letter, ofn.legend.letter)
    save.table(df.color, ofn.legend.color)
    save.table(result[,c("anchor.id", "anchor", "group.id", "sub.group.id", "is.cag", "color", "letter")], ofn)
}
