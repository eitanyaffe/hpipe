seq.summary=function(order.ifn, table.ifn, gene.table.ifn, anchor2ref.dir, ref2anchor.dir, length, ofn, ofn.unique)
{
    anchor.table = load.table(order.ifn)
    table = load.table(table.ifn)
    gtable = load.table(gene.table.ifn)
    gtable$id = gtable[,1]

    get.coverage=function(ifn) {
        tt = read.delim(ifn)
        round(100*sum(tt$count[tt$edit != -1]) / sum(tt$count),1)
    }
    get.identity=function(ifn) {
        tt = read.delim(ifn)
        tt = tt[tt$edit != -1,]
        tt$identity = 100 - 100*tt$edit/length
        tt$weight = tt$count / sum(tt$count)
        round(sum(tt$identity*tt$weight),3)
    }
    get.length=function(ifn) {
        tt = read.delim(ifn)
        sum(tt$count)
    }
    result = NULL
    for (i  in 1:dim(table)[1]) {
        anchor = table$anchor[i]
        accession = table$accession[i]
        src.anchor.ifn = paste(anchor2ref.dir, "/", anchor, "_", accession, "/src_summary", sep="")
        tgt.ref.ifn = paste(anchor2ref.dir, "/", anchor, "_", accession, "/tgt_summary", sep="")
        src.ref.ifn = paste(ref2anchor.dir, "/", anchor, "_", accession, "/src_summary", sep="")
        tgt.anchor.ifn = paste(ref2anchor.dir, "/", anchor, "_", accession, "/tgt_summary", sep="")
        df = data.frame(
            anchor=anchor, ref=accession,
            anchor.length=get.length(src.anchor.ifn),
            ref.length=get.length(src.ref.ifn),
            src.anchor.coverage=get.coverage(src.anchor.ifn),
            tgt.anchor.coverage=get.coverage(tgt.anchor.ifn),
            src.ref.coverage=get.coverage(src.ref.ifn),
            tgt.ref.coverage=get.coverage(tgt.ref.ifn),
            src.anchor.identity=get.identity(src.anchor.ifn),
            tgt.anchor.identity=get.identity(tgt.anchor.ifn),
            src.ref.identity=get.identity(src.ref.ifn),
            tgt.ref.identity=get.identity(tgt.ref.ifn))
        result = rbind(result, df)
    }
    result$anchor.coverage = round((result$src.anchor.coverage + result$tgt.anchor.coverage)/2,1)
    result$ref.coverage = round((result$src.ref.coverage + result$tgt.ref.coverage)/2,1)
    result$anchor.identity = round((result$src.anchor.identity + result$tgt.anchor.identity)/2,1)
    result$ref.identity = round((result$src.ref.identity + result$tgt.ref.identity)/2,1)
    result$anchor.id = anchor.table$id[match(result$anchor, anchor.table$set)]
    save.table(result, ofn)

    anchor.ids = anchor.table$id
    result.u = NULL
    for (anchor.id in anchor.ids) {
        xx = result[result$anchor.id == anchor.id,]
        ix = which.max(xx$anchor.coverage)
        result.u = rbind(result.u, xx[ix,
            c("anchor.id", "anchor", "ref", "anchor.length", "ref.length", "anchor.coverage", "ref.coverage", "anchor.identity", "ref.identity")])
    }
    result.u$ref.name = gtable$organism_name[match(result.u$ref, gtable$id)]
    result.u$ref.strain.name = gtable$infraspecific_name[match(result.u$ref, gtable$id)]
    result.u$is.cag = grepl("CAG", result.u$ref.name)
    save.table(result.u, ofn.unique)
}
