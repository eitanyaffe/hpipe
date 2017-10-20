anchor.stats=function(ifn, ifn.contigs, odir)
{
    ca = load.table(ifn)
    ctable = load.table(ifn.contigs)
    ca$length = lookup(table=ca, lookup.table=ctable, lookup.field="contig", value.field="length", na.value=NA)
    anchors = sort(unique(ca$anchor))
    df = field.count(ca, "contig")
    df$length = lookup(table=df, lookup.table=ctable, lookup.field="contig", value.field="length", na.value=NA)
    caa = ca[ca$anchor == ca$contig_anchor,]

    median.anchor.bp= median(sapply(split(caa$length, caa$anchor), sum))
    median.cloud.bp= median(sapply(split(ca$length, ca$anchor), sum))
    shared.bp = sum(df$length[df$count>1])
    shared.count = sum(df$count>1)
    shared3.count = sum(df$count>2)
    total.count = dim(df)[1]

    lines = NULL

    lines = c(lines, sprintf("number of anchors: %d", length(anchors)))
    lines = c(lines, sprintf("number of anchor contigs: %d", dim(caa)[1]))
    lines = c(lines, sprintf("number of cloud contigs: %d", dim(df)[1]))

    lines = c(lines, sprintf("median anchor: %.2fMbp", median.anchor.bp/10^6))
    lines = c(lines, sprintf("median cloud: %.2fMbp", median.cloud.bp/10^6))

    lines = c(lines, sprintf("anchors: %.2fMbp (%.1f%%)", sum(caa$length)/10^6, 100*sum(caa$length)/sum(ctable$length)))
    lines = c(lines, sprintf("clouds: %.2fMbp (%.1f%%)", sum(df$length)/10^6, 100*sum(df$length)/sum(ctable$length)))


    lines = c(lines, sprintf("shared: %.2fKbp (%.1f%% of assembly)", shared.bp/10^3, 100*shared.bp/sum(df$length)))
    lines = c(lines, sprintf("shared count(>1): %d (%.1f%% of contigs)", shared.count, 100*shared.count/total.count))
    lines = c(lines, sprintf("shared count(>2): %d (%.1f%% of contigs)", shared3.count, 100*shared3.count/total.count))

    system(paste("mkdir -p", odir))
    ofn = paste(odir, "/anchor_info.txt", sep="")
    cat(sprintf("generating file: %s\n", ofn))
    fc = file(ofn)
    writeLines(lines, fc)
    close(fc)
}

anchor.stats.table=function(ifn, ifn.contigs, ofn)
{
    ca = load.table(ifn)
    ctable = load.table(ifn.contigs)
    ca$length = lookup(table=ca, lookup.table=ctable, lookup.field="contig", value.field="length", na.value=NA)
    anchors = sort(unique(ca$anchor))
    caa = ca[ca$anchor == ca$contig_anchor,]

    get.values=function(table) {
        s = sapply(split(table$length, table$anchor), sum)
        ix = match(anchors,names(s))
        ifelse(!is.na(ix), s[ix], 0)
    }

    result = data.frame(anchor=anchors, anchor.size=get.values(caa), union.size=get.values(ca))
    save.table(result, ofn)
}
