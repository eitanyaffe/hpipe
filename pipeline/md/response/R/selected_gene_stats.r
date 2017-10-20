selected.gene.stats=function(ifn.genes, ifn.bg, poor.annotation.desc, ofn)
{
    gene = load.table(ifn.genes)
    counts.bg = load.table(ifn.bg)

    for (i in 1:length(poor.annotation.desc))
        poor.annotation.desc[i] = sub("_", " ", poor.annotation.desc[i])
    gene$prot_desc = tolower(gene$prot_desc)
    poor.annotation.desc = tolower(poor.annotation.desc)
    identity.below.70 = sum(!is.na(gene$identity) & gene$identity < 70)
    counts = data.frame(total=dim(gene)[1], no.uniref=sum(is.na(gene$uniref)), poor.annotation=sum(is.element(gene$prot_desc, poor.annotation.desc)), identity.below.70=identity.below.70)
    percent = 100 * counts[,-1] / counts[,1]
    percent.bg = 100 * counts.bg[,-1] / counts.bg[,1]

    df = as.data.frame(t(rbind(counts[,-1], percent, counts.bg[,-1], percent.bg, percent / percent.bg)))
    colnames(df) = c("counts", "percent", "counts.bg", "percent.bg", "enrichment")
    df = data.frame(type=rownames(df), df)
    rownames(df) = NULL
    df = rbind(df, data.frame(type="no.annotation.sum", df[df$type == "no.uniref",-1] + df[df$type == "poor.annotation",-1]))

    # fix enrichment
    ix = df$type == "no.annotation.sum"
    df$enrichment[ix] = df$percent[ix] / df$percent.bg[ix]
    save.table(df, ofn)
}
