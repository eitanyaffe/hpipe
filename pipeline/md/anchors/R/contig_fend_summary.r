contig.fend.summary=function(ifn.fends, ifn.contigs, ofn)
{
    df = load.table(ifn.fends)
    contigs = load.table(ifn.contigs)$contig

    s = sapply(split(df$fend, df$contig), length)
    ix = match(contigs, names(s))
    count = ifelse(!is.na(ix),s[ix],0)
    result = data.frame(contig=contigs, count=count)

    save.table(result, ofn)
}
