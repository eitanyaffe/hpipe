
merge.GO=function(ifn.genes, ifn.genes.ctrl, ifn.prefix, ifn.prefix.ctrl, ofn)
{
    genes = load.table(ifn.genes)
    genes.ctrl = load.table(ifn.genes.ctrl)

    ll = list(
        func = load.table(paste(ifn.prefix, "_function", sep="")),
        process = load.table(paste(ifn.prefix, "_process", sep="")),
        component = load.table(paste(ifn.prefix, "_component", sep="")))

    ll.ctrl = list(
        func = load.table(paste(ifn.prefix.ctrl, "_function", sep="")),
        process = load.table(paste(ifn.prefix.ctrl, "_process", sep="")),
        component = load.table(paste(ifn.prefix.ctrl, "_component", sep="")))

    result = NULL
    for (i in 1:length(ll)) {
        name = names(ll)[i]
        table = ll[[i]]
        table.ctrl = ll.ctrl[[i]]
        m = merge(table, table.ctrl, by="id")
        df = data.frame(id=m$id, desc=m$desc.x, count=m$count.x, count.ctrl=m$count.y)
        df$percent = df$count / dim(genes)[1]
        df$percent.ctrl = df$count.ctrl / dim(genes.ctrl)[1]
        df$enrichment = df$percent / df$percent.ctrl
        result = rbind(result, data.frame(type=name, df))
    }
    save.table(result, ofn)
}
