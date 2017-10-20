top = function(table, ofn, min.length)
{
    cat(sprintf("reading table: %s\n", table))
    table = read.delim(table)
    N = dim(table)[1]
    table = table[order(table$length, decreasing=T),]
    result = table[table$length >= min.length,]
    cat(sprintf("selected contig count: %d/%d (%.1f%%)\n", dim(result)[1], N, 100*dim(result)[1]/N))
    cat(sprintf("selected length bp: %d/%d (%.1f%%)\n", sum(result$length), sum(table$length), 100*sum(result$length)/sum(table$length)))
    cat(sprintf("contig length stats: min=%d ,max=%d, total=%d\n", min(result$length) ,max(result$length), sum(result$length)))
    cat(sprintf("generating file: %s\n", ofn))
    write.table(result, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}
