plot.shared.analysis=function(ifn.uniref, ifn.genes, ifn.legend, fdir)
{
    options(stringsAsFactors=F)
    uniref = load.table(ifn.uniref)
    classes = list(phage="phage", plasmid=c("parb", "plasmid", "traa"), conjugation=c("moba"), transposase="transposase", toxin="toxin")
    uniref$class = "none"
    for (i in 1:length(classes)) {
        class = names(classes)[i]
        for (pattern in classes[[i]])
            uniref$class = ifelse(uniref$class == "none" & grepl(pattern=pattern, x=uniref$prot_desc, ignore.case=T), class, uniref$class)
    }

    # count genes per anchor
    genes = load.table(ifn.genes)
    df = field.count(genes, "gene")

    mm = match(genes$gene, uniref$gene)
    genes$identity = ifelse(!is.na(mm), uniref$identity[mm], 0)
    genes$class = ifelse(!is.na(mm), uniref$class[mm], "none")
    genes$is.mobile = genes$class != "none"
    genes$multi = df$count[match(genes$gene, df$gene)]


    # multiplicity counts
    tt = table(df$count)
    tt = tt[-1]
    fig.start(fdir=fdir, ofn=paste(fdir, "/multi.png", sep=""), height=300, width=200)
    barplot(tt, col="blue")
    title(main="anchors per gene", ylab="#genes")
    fig.end()

    # hist of identity
    breaks = 0:10*10
    s = split(genes$identity, ifelse(genes$multi == 1, "normal", "shared"))
    sx = sapply(s, function(x) { 100 * sapply(split(x, cut(x, breaks)), length) / length(x) })
    fig.start(fdir=fdir, ofn=paste(fdir, "/identity.png", sep=""), height=600, width=600)
    barplot(sx, beside=T, border=NA, col="gray")
    abline(h=0)
    title(main="identity", ylab="%")
    fig.end()

    s = sapply(split(genes$is.mobile, ifelse(genes$multi == 1, "normal", "shared")), function(x) { 100 * sum(x) / length(x) } )
    fig.start(fdir=fdir, ofn=paste(fdir, "/mobile_precentage.png", sep=""), height=300, width=200)
    barplot(s, border=NA, col="red")
    title(main="%mobile", ylab="%")
    fig.end()

    s = sapply(split(genes$is.mobile, genes$multi), function(x) { sum(x) } )
    fig.start(fdir=fdir, ofn=paste(fdir, "/mobile_counts.png", sep=""), height=300, width=200)
    barplot(s, border=NA, col="red")
    title(main="%mobile", ylab="%")
    fig.end()
}
