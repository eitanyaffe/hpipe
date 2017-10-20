
multi.summary=function(ifn, assembly.dir, ids, map.tag, ofn)
{
    options(stringsAsFactors=F)
    contigs = load.table(ifn)
    ifns = paste(assembly.dir, "/datasets/", ids, "/map_", map.tag, "/var_contig_summary", sep="")

    for (i in 1:length(ids)) {
        ifn = ifns[i]
        id = ids[i]
        x = load.table(ifn)
        if (i == 1)
            result = x[,c("contig", "length")]
        result[[id]] = x$total
    }
    save.table(result, ofn)
}
