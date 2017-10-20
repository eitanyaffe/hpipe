# deprecated, divides contigs into bins
anchor.coverage.rpkb=function(ca.ifn, idir, binsize, ofn)
{
    ca = load.table(ca.ifn)
    anchors = sort(unique(ca$anchor))

    result = NULL
    for (anchor in anchors) {
        contigs = ca$contig[ca$anchor == anchor]
        cat(sprintf("processing anchor: %d\n", anchor))
        ra = NULL
        length = 0
        count = 0
        for (contig in contigs) {
            x = read.delim(paste(idir, "/", contig, sep=""))
            x$bin = floor(x$coord / binsize)
            ll = max(x$bin) * binsize
            if (ll <= 0)
                next
            x = x[x$coord < ll,]
            ra = c(ra, as.vector(sapply(split(x$count, x$bin), sum)))
            length = length + max(x$bin+1)*binsize
            count = count +  max(x$bin+1)
        }
        rpkb = ra / (length/1000)
        df = data.frame(
            anchor=anchor,
            total.reads=sum(ra),
            count=count,
            bottom05=quantile(rpkb, 0.05),
            bottom25=quantile(rpkb, 0.25),
            median=median(rpkb),
            top75=quantile(rpkb, 0.75),
            top95=quantile(rpkb, 0.95))

        result = rbind(result, df)
    }
    save.table(result, ofn)
}

anchor.coverage=function(ca.ifn, cov.ifn, ofn)
{
    cov = load.table(cov.ifn)
    ca = load.table(ca.ifn)
    anchors = sort(unique(ca$anchor))

    result = NULL
    for (anchor in anchors) {
        contigs = ca$contig[ca$anchor == anchor]
        ra = cov$abundance.enrichment[match(contigs,cov$contig)]
        df = data.frame(
            anchor=anchor,
            count=length(ra),
            bottom05=round(quantile(ra, 0.05),3),
            bottom25=round(quantile(ra, 0.25),3),
            median=round(median(ra),3),
            top75=round(quantile(ra, 0.75),3),
            top95=round(quantile(ra, 0.95),3))
        result = rbind(result, df)
    }
    save.table(result, ofn)
}

anchor.gc=function(ca.ifn, ifn.gc, binsize, ofn)
{
    ca = load.table(ca.ifn)
    anchors = sort(unique(ca$anchor))
    gc.table = load.table(ifn.gc)

    result = NULL
    for (anchor in anchors) {
        contigs = ca$contig[ca$anchor == anchor]
        gc = gc.table[is.element(gc.table$contig, contigs),"gc"]
        df = data.frame(anchor=anchor, count=length(gc),
            bottom05=quantile(gc, 0.05),
            bottom25=quantile(gc, 0.25),
            median=median(gc),
            top75=quantile(gc, 0.75),
            top95=quantile(gc, 0.95))

        result = rbind(result, df)
    }
    save.table(result, ofn)
}

anchor.info=function(
    ga.ifn, ca.ifn,
    contig.ifn, gc.ifn, coverage.ifn,
    ofn.summary, ofn.genes, ofn.contigs)
{
    options(stringsAsFactors=F)

    ga = read.delim(ga.ifn)
    gc = read.delim(gc.ifn)
    coverage = read.delim(coverage.ifn)

    df = field.count(ga, "gene")
    ga$count = df$count[match(ga$gene, df$gene)]

    gene.fields = c("anchor", "enrichment", "count")
    s = split(ga, ga$anchor)
    genes.info = data.frame(
        anchor=as.numeric(names(s)),
        enrichment=sapply(s, function(x) median(x$enrichment)),
        gene.count=sapply(s, function(x) dim(x)[1]),
        shared.genes=sapply(s, function(x) sum(x$count > 1)))
    save.table(ga[,gene.fields], ofn.genes)

    contigs = read.delim(contig.ifn)
    contigs$anchor = ga$anchor[match(contigs$contig, ga$contig)]
    contigs = contigs[!is.na(contigs$anchor),]
    contigs$gc = gc$gc[match(contigs$contig, gc$contig)]
    contigs$reads = coverage$sum[match(contigs$contig, coverage$contig)]

    s = split(contigs, contigs$anchor)
    contigs.info = data.frame(
        anchor=as.numeric(names(s)),
        length=sapply(s, function(x) sum(x$length)),
        median.contig.length=sapply(s, function(x) median(x$length)),
        gc=sapply(s, function(x) median(x$gc)),
        contig.count=sapply(s, function(x) dim(x)[1]),
        reads=sapply(s, function(x) sum(x$reads)))
    contigs.info$reads.per.bp = contigs.info$reads / contigs.info$length
    contigs.info$coverage = contigs.info$reads.per.bp / sum(contigs.info$reads.per.bp)

    contig.fields = c("anchor", "gc", "reads", "length")
    save.table(contigs[,contig.fields], ofn.contigs)

    info = merge(genes.info, contigs.info)
    cat(sprintf("saving result table: %s\n", ofn.summary))
    write.table(info, ofn.summary, quote=F, col.names=T, row.names=F, sep="\t")
}

plot.anchor.info.summary=function(ifn, fdir)
{
    library(gplots)
    options(stringsAsFactors=F)

    table = read.delim(ifn)
    fields = c("shared.genes", "length", "median.contig.length", "coverage")

    table$shared.genes = table$shared.genes + 1
    dta = log10(as.matrix(table[,fields]))
    system(paste("mkdir -p", fdir))

    # complete
    ofn = paste(fdir, "/complete_matrix.png", sep="")
    cat(sprintf("plotting figure: %s\n", ofn))
    png(ofn, width=800, height=800)
    plot.new()
    panel.cor = function(x, y, digits = 2, prefix = "", cex.cor, ...)
    {
        points(x,y,...)
        usr = par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r = abs(cor(x, y))
        txt = format(c(r, 0.123456789), digits = digits)[1]
        txt = paste0(prefix, txt)
        text(0.5, 0.5, txt, cex=2)
    }
    pairs(dta, fields, lower.panel=panel.cor, upper.panel=panel.cor, pch="+")
    dev.off()

    make.color = function(x,
        ncols=256, panel=colorpanel(ncols, "blue", "white", "red"), start=min(x), end=max(x))
    {
        v = (x-start) / (end - start)
        vi = v*(ncols-1) + 1
        panel[vi]
    }

    ofn = paste(fdir, "/genes.png", sep="")
    cat(sprintf("plotting figure: %s\n", ofn))
    png(ofn, width=600, height=600)
    barplot(table$gene.count, names.arg=table$anchor, las=2)
    dev.off()

    ofn = paste(fdir, "/shared_genes.png", sep="")
    cat(sprintf("plotting figure: %s\n", ofn))
    png(ofn, width=600, height=600)
    barplot(table$shared.genes, names.arg=table$anchor, las=2)
    dev.off()

    # fragmented vs. abundance
    ofn = paste(fdir, "/fragment_vs_shared.png", sep="")
    cat(sprintf("plotting figure: %s\n", ofn))
    png(ofn, width=600, height=600)
    plot.new()
    plot(table$median.contig.length, table$shared.genes+1,
         type="p", pch=19, col=make.color(table$length), log="xy",
         xlab="median contig length", ylab="#shared genes", main="color=length")
    points(table$median.contig.length, table$shared.genes+1, pch=1)
    dev.off()

    # enrichment vs. abundance
    ofn = paste(fdir, "/enrichment_vs_abundance.png", sep="")
    cat(sprintf("plotting figure: %s\n", ofn))
    png(ofn, width=600, height=600)
    plot.new()
    plot(table$coverage, table$enrichment,
         type="p", pch=19, col=1, log="x",
         xlab="coverage", ylab="enrichment")
    dev.off()

    ofn = paste(fdir, "/fragment_vs_abundance.png", sep="")
    cat(sprintf("plotting figure: %s\n", ofn))
    png(ofn, width=600, height=600)
    plot.new()
    plot(table$median.contig.length, table$coverage,
         type="p", pch=19, col=make.color(table$length), log="xy",
         xlab="median contig length", ylab="abundance", main="color=length")
    points(table$median.contig.length, table$coverage)
    dev.off()

    ofn = paste(fdir, "/abundance_length.png", sep="")
    cat(sprintf("plotting figure: %s\n", ofn))
    png(ofn, width=600, height=600)
    plot.new()
    plot(table$coverage, table$length,
         type="p", pch=19, col=make.color(table$median.contig.length), log="xy",
         xlab="coverage", ylab="length", main="color=contig length")
    points(table$coverage, table$length)
    dev.off()

    ofn = paste(fdir, "/fragment_length.png", sep="")
    cat(sprintf("plotting figure: %s\n", ofn))
    png(ofn, width=600, height=600)
    plot.new()
    col = make.color(log10(table$shared.genes+1), panel=colorpanel(256, "white", "red"),
        start=0, end=max(log10(table$shared.genes+1)))
    plot(log10(table$median.contig.length), table$length,
         type="p", pch=19, col=col, log="x",
         xlab="log10(median contig length)", ylab="length", main="color=shared")
    points(log10(table$median.contig.length), table$length)
    dev.off()
}

plot.anchor.info=function(ifn, ifn.genes, ifn.contigs, fdir)
{
    library(gplots)
    options(stringsAsFactors=F)

    genes = read.delim(ifn.genes)
    contigs = read.delim(ifn.contigs)
    contigs$reads.per.bp = contigs$reads / contigs$length

    pp=function(field, table, ofn, sort=F, log.y=F) {
        values = split(table[,field], table$anchor)
        if (sort)
            values = values[order(sapply(values, median))]
        cat(sprintf("plotting figure: %s\n", ofn))
        png(ofn, width=800, height=600)
        boxplot(values, main=field, outline=F, las=2, log=ifelse(log.y,"y",""))
        dev.off()
    }
    system(paste("mkdir -p", fdir))

    pp(table=genes, field="enrichment", ofn = paste(fdir,"/full_contact_enrichment.png", sep=""))
    pp(table=contigs, field="reads.per.bp", ofn = paste(fdir,"/full_coverage.png", sep=""), log.y=T)
    pp(table=contigs, field="length", ofn = paste(fdir,"/full_contig_length.png", sep=""))
    pp(table=contigs, field="gc", ofn = paste(fdir,"/full_gc.png", sep=""))

    pp(table=genes, field="enrichment", ofn = paste(fdir,"/full_sort_contact_enrichment.png", sep=""), sort=T)
    pp(table=contigs, field="reads.per.bp", ofn = paste(fdir,"/full_sort_coverage.png", sep=""), sort=T, log.y=T)
    pp(table=contigs, field="length", ofn = paste(fdir,"/full_sort_contig_length.png", sep=""), sort=T)
    pp(table=contigs, field="gc", ofn = paste(fdir,"/full_sort_gc.png", sep=""), sort=T)
}

plot.report=function(
    ifn, contig.ifn, gc.ifn, coverage.ifn1, coverage.ifn2, gene.ifn, uniref.ifn, response.ifn, contig.cluster.ifn, fdir)
{
    classes = list(phage="phage", plasmid=c("parb", "plasmid", "traa"), conjugation=c("moba"), transposase="transposase", toxin="toxin")
    class.color = c("red", "blue", "purple", "green", "orange")
    wlegend(fdir=fdir, names=names(classes), cols=class.color, title="mobile", width=300)

    options(stringsAsFactors=F)
    tab = load.table(ifn)

    gc = read.delim(gc.ifn)
    contigs = read.delim(contig.ifn)
    coverage1 = read.delim(coverage.ifn1)
    coverage2 = read.delim(coverage.ifn2)

    rt = load.table(response.ifn)
    M = dim(rt)[2] - 1

    rcolors = c("black", "blue", "white", "red", "orange")
    panel = make.color.panel(colors=rcolors)
    breaks = c(-3, -1, 0, 1, 2)

    # contig clusters
    ccluster = load.table(contig.cluster.ifn)
    ccount = field.count(ccluster, "cluster")
    ccluster$size = ccount$count[match(ccluster$cluster, ccount$cluster)]
    tab$ccluster = ccluster$cluster[match(tab$contig, ccluster$contig)]
    tab$csize = ccluster$size[match(tab$contig, ccluster$contig)]

    genes = load.table(gene.ifn)
    uniref = load.table(uniref.ifn)
    uniref$contig = genes$contig[match(uniref$gene, genes$gene)]
    uniref$class = "none"
    for (i in 1:length(classes)) {
        class = names(classes)[i]
        for (pattern in classes[[i]])
            uniref$class = ifelse(uniref$class == "none" & grepl(pattern=pattern, x=uniref$prot_desc, ignore.case=T), class, uniref$class)
    }
    xx = uniref[uniref$class != "none",]
    s = split(xx, xx$contig)
    s = sapply(s, function(x) {
               cs = unique(x$class)
               if (length(cs) > 1)
                   "multi"
               else
                   cs
           } )

    tab$class = ifelse(is.na(match(tab$contig, names(s))), "none", s[match(tab$contig, names(s))])
    tab$class.color = ifelse(tab$class == "none", "lightgray", class.color[match(tab$class, names(classes))])

    tab$length = contigs$length[match(tab$contig, contigs$contig)]
    tab$gc = gc$gc[match(tab$contig, gc$contig)]

    tab$reads1 = coverage1$sum[match(tab$contig, coverage1$contig)]
    tab$reads.per.bp1 = tab$reads1 / tab$length

    tab$reads2 = coverage2$sum[match(tab$contig, coverage2$contig)]
    tab$reads.per.bp2 = tab$reads2 / tab$length

    df = field.count(tab, "contig")
    tab$anchor.count = df$count[match(tab$contig, df$contig)]

    tab$basic.color = ifelse(tab$anchor.count>1, "orange", ifelse(tab$contig_anchor == tab$anchor, "red", "black"))

    anchors = sort(unique(tab$anchor))
    for (anchor in anchors) {
        fig.start(fdir=fdir, ofn=paste(fdir, "/", anchor, ".png", sep=""), width=1000, height=600)

        x = tab[tab$anchor == anchor,]
        x = x[order(x$reads.per.bp1),]
        x$index = 1:dim(x)[1]

        xlim = c(0, max(x$index))
        heights = c(8,8,8,8,4,1,1,1)
        N = length(heights)

        layout(matrix(1:N, N, 1), heights=heights)
        par(mai=c(0.1,1,0.1,0.2))

        # 1 reads per bp
        plot.init(xlim=xlim, ylim=range(1000*tab$reads.per.bp1), x.axis=F, ylab="reads/kbp", log="y")
        lines(x$index, 1000*x$reads.per.bp1)
        lines(x$index, 1000*x$reads.per.bp2, lty=2)

        # 2 contig length
        plot.init(xlim=xlim, ylim=range(tab$length/1000), x.axis=F, ylab="length(kb)", log="y")
        lines(x$index, x$length/1000)

        # 3 enrichment
        plot.init(xlim=xlim, ylim=c(0, max(tab$enrichment)), x.axis=F, ylab="cc score")
        lines(x$index, x$enrichment)

        # 4 gc
        plot.init(xlim=xlim, ylim=range(tab$gc), x.axis=F, ylab="gc")
        lines(x$index, x$gc)

        # 5 response
        rr = rt[is.element(rt$contig, x$contig),]
        rr = rr[match(x$contig, rr$contig),]
        par(mai=c(0.01,1,0.02,0.2))
        plot.new()
        plot.window(xlim=xlim, ylim=c(0,M))
        for (i in 1:M) {
            vals = rr[,1+i]
            cols = ifelse(is.finite(vals), panel[vals.to.cols(vals=vals,breaks=breaks)], "black")
            rect(xleft=x$index-1, xright=x$index, ybottom=i-1, ytop=i, col=cols, border=NA)
        }
        rect(xleft=0, xright=max(x$index), ybottom=0, ytop=M+1, col=NA, border=1)

        # 6 shared
        par(mai=c(0.01,1,0.02,0.2))
        plot.new()
        plot.window(xlim=xlim, ylim=c(0,1))
        rect(xleft=x$index-1, xright=x$index, ybottom=0, ytop=1, col=x$basic.color, border=NA)

        # 7 class
        par(mai=c(0.01,1,0.02,0.2))
        plot.new()
        plot.window(xlim=xlim, ylim=c(0,1))
        rect(xleft=x$index-1, xright=x$index, ybottom=0, ytop=1, col=x$class.color, border=NA)

        # 8 cluster
        par(mai=c(0.01,1,0.02,0.2))
        plot.new()
        plot.window(xlim=xlim, ylim=c(0,1))
        if (any(x$csize>1)) {
            clusters = setdiff(sort(unique(x$ccluster[x$csize>1])), -1)
            clusters.colors = rainbow(length(clusters))
            x$ccolor = ifelse(x$ccluster == -1, "white", clusters.colors[match(x$ccluster, clusters)])
            wlegend(fdir=fdir,
                    names=paste(clusters, " : ", ccount$count[match(clusters, ccount$cluster)], sep=""),
                    cols=clusters.colors, title=paste(anchor, "_clusters", sep=""), width=300)
            rect(xleft=x$index-1, xright=x$index, ybottom=0, ytop=1, col=x$ccolor, border=NA)
        }

        fig.end()
    }
}
