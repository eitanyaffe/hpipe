get.temporal.table=function(tables, ids, contigs, field="sum", baseline.idx=1)
{
    N = length(tables)

    cat(sprintf("uniting %d profiles\n", length(tables)))
    table = data.frame(contig=contigs)

    for (i in 1:N) {
        p = load.table(tables[i])
        p$value = p[,field]
        df = data.frame(contig=p$contig, value=p$value)
        id = ids[i]
        table[[id]] = p$value[match(contigs, p$contig)]
    }
    table
}

compute.contig.response=function(ifn, assembly.dir, ids, ofn.observed, ofn.expected, ofn.norm)
{
    options(stringsAsFactors=F)
    contigs = load.table(ifn)
    ttables = paste(assembly.dir, "/datasets/", ids, "/coverage/table", sep="")

    observed.table = get.temporal.table(ttables, ids=ids, contigs=contigs$contig, field="sum")
    save.table(observed.table, ofn.observed)

    # uniform expected values
    expected.table = data.frame(contig=contigs$contig)
    total.length = sum(contigs$length)
    for (i in 2:dim(observed.table)[2]) {
        total.reads = sum(observed.table[,i])
        expected = total.reads * (contigs$length / total.length)
        expected.table = cbind(expected.table, expected)
    }
    names(expected.table) = names(observed.table)
    save.table(observed.table, ofn.expected)

    # uniform expected values
    norm.table = data.frame(contig=contigs$contig, log10((observed.table[,-1]+1) / (expected.table[,-1]+1)))
    save.table(norm.table, ofn.norm)
}

