##################################################################################################
# generate read coords
##################################################################################################

generate.abundance.table=function(ifn, abundance.min, abundance.max, abundance.min.hic, abundance.max.hic, ofn)
{
    options(stringsAsFactors=F)
    table = read.delim(ifn)
    table$abundance = runif(dim(table)[1], min=abundance.min, max=abundance.max)
    table$abundance.hic = runif(dim(table)[1], min=abundance.min.hic, max=abundance.max.hic)

    cat(sprintf("saving table: %s\n", ofn))
    write.table(table, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}
