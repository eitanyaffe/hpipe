anchor.genes=function(ifn, ifn.cgenes, ofn, min.contacts=10, min.contig.coverage=0.5, min.anchor.contigs=2, q.threshold=0.97)
{
  library(gplots)
  options(stringsAsFactors=F)

  cat(sprintf("reading cgenes table: %s\n", ifn.cgenes))
  cgenes = read.delim(ifn.cgenes)

  cat(sprintf("reading table: %s\n", ifn))
  table = read.delim(ifn)
  table$type = ifelse(table$gene_anchor == 0, "extended", ifelse(table$gene_anchor == table$anchor, "intra", "inter"))
  table$cgene = cgenes$cgene[match(table$gene, cgenes$gene)]

  # contact enrichment
  table$enc = ifelse(table$contig_total_count == 0, -Inf, log10(table$contig_total_count/table$contig_expected))
  table$eng = ifelse(table$gene_count == 0, -Inf, log10(table$gene_count/table$gene_expected))

  # contact enrichment threshold
  threshold = quantile(table$enc[table$type == "inter" & is.finite(table$enc)], q.threshold)
  cat(sprintf("input threshold q: %f\n", q.threshold))
  cat(sprintf("selected o/e enrichment threshold: %f\n", 10^threshold))

  # selected genes
  gene.selected = (table$gene_count>=min.contacts & table$eng>threshold &
                   table$anchor_contig_count>=min.anchor.contigs)

  # selected contigs
  contig.selected =
      (table$contig_total_count>=min.contacts & table$enc>threshold &
       table$contig_coverage>=min.contig.coverage & table$anchor_contig_count>=min.anchor.contigs)

  N = sum(contig.selected)

  table$support = ifelse(gene.selected & contig.selected, "both", ifelse(contig.selected, "contig", "none"))
  table$gene.selected = gene.selected
  table$contig.selected = contig.selected
  cat(sprintf("number of associated genes: %d (by contig=%d)\n",
              N, sum(contig.selected&!gene.selected)))

  # TBD select only by contig
  table = table[contig.selected,]

  table$enrichment = ifelse(table$gene.selected, table$eng, table$enc)

  cat(sprintf("saving result table: %s\n", ofn))
  write.table(table, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}

plot.ga.distrib=function(ifn, fdir, min.contacts=10, min.contig.coverage=0.5, min.anchor.contigs=2, q.threshold=0.97)
{
  library(gplots)
  options(stringsAsFactors=F)

  cat(sprintf("reading table: %s\n", ifn))
  table = read.delim(ifn)
  table$type = ifelse(table$gene_anchor == 0, "extended", ifelse(table$gene_anchor == table$anchor, "intra", "inter"))

  # contact enrichment
  table$enc = ifelse(table$contig_total_count == 0, -Inf, log10(table$contig_total_count/table$contig_expected))
  table$eng = ifelse(table$gene_count == 0, -Inf, log10(table$gene_count/table$gene_expected))

  # contact enrichment threshold
  threshold = quantile(table$enc[table$type == "inter" & is.finite(table$enc)], q.threshold)
  cat(sprintf("input threshold q: %f\n", q.threshold))
  cat(sprintf("selected o/e enrichment threshold: %f\n", 10^threshold))

  # selected genes
  gene.selected = (table$gene_count>=min.contacts & table$eng>threshold &
                   table$anchor_contig_count>=min.anchor.contigs)

  # selected contigs
  contig.selected =
      (table$contig_total_count>=min.contacts & table$enc>threshold &
       table$contig_coverage>=min.contig.coverage & table$anchor_contig_count>=min.anchor.contigs)

  table$gene.selected = gene.selected
  table$contig.selected = contig.selected
  table$selected = contig.selected

  scatter.fdir = paste(fdir, "/anchor_scatter", sep="")
  system(paste("mkdir -p", fdir, scatter.fdir))
  cat(sprintf("plotting in directory: %s\n", fdir))

  ################################################################
  # density
  ################################################################

  fin = is.finite(table$enc)
  inter.ind = table$type == "inter" & fin & table$contig_total_count>=min.contacts
  intra.ind = table$type == "intra" & fin & table$contig_total_count>=min.contacts

  pdensity=function(field, title) {
      d.inter = density(table[inter.ind, field])
      d.intra = density(table[intra.ind, field])

      xlim = range(c(d.inter$x, d.intra$x))
      ylim = range(c(d.inter$y, d.intra$y))

      ofn = paste(fdir, "/", title, "_density.png", sep="")
      cat(sprintf("plotting figure: %s\n", ofn))
      png(ofn, width=400, height=400)
      main = paste(title, "-anchor contact enrichment density", sep="")
      plot(d.inter, xlab="log10(o/e)", lwd=2, xlim=xlim, ylim=ylim, col="black", main=main)
      lines(d.intra, col="red", lwd=2)
      abline(v=threshold, lty=2)
      legend("topright", fill=c("black", "red"), legend=c("inter", "intra"))
      dev.off()
  }
  pdensity(field="enc", title="contig")
  # pdensity(field="eng", title="gene")

  ################################################################
  # scatters obs vs exp
  ################################################################

  xlim = c(1, max(table$contig_total_count))
  ylim = range(c(table$gene_expected, table$contig_expected))

  x = seq(from=min.contacts, to=xlim[2], by=1)
  y = x / 10^threshold
  x = c(min.contacts, x)
  y = c(ylim[1], y)

  table$col =
      ifelse(table$type == "inter", "blue",
             ifelse(table$type=="intra", "black",
                    ifelse(table$gene.selected, "red",
                           ifelse(table$contig.selected, "orange", "grey"))))

  ofn = paste(fdir, "/contig_scatter.png", sep="")
  cat(sprintf("plotting figure: %s\n", ofn))
  png(ofn, width=600, height=600)
  plot(table$contig_total_count, table$contig_expected, pch=".", type="p", log="xy", cex=4,
       xlab="observed", ylab="expected", col=table$col, main="contig-anchor contact enrichment")
  lines(x, y, lty=2, lwd=2)
  dev.off()

  ofn = paste(fdir, "/gene_scatter.png", sep="")
  png(ofn, width=600, height=600)
  plot(table$gene_count, table$gene_expected, pch=".", type="p", log="xy", cex=4,
       xlab="observed", ylab="expected", col=table$col, main="gene-anchor contact enrichment")
  lines(x, y, lty=2, lwd=2)
  dev.off()

  pscatter=function(ttable, ofn, title, field.o, field.e)
  {
      ttable = ttable[ttable[,field.o] > 0 & ttable[,field.e] > 0,]
      png(ofn, width=600, height=600)

      main = title
      plot(ttable[,field.o], ttable[,field.e], pch=".", type="p", log="xy", cex=4, xlim=xlim, ylim=ylim,
           xlab="observed", ylab="expected", col=ttable$col, main=main)
      lines(x, y, lty=2, lwd=2)
      dev.off()
  }

  anchors = unique(table$anchor)
  for (anchor in anchors) {
      ttable = table[table$anchor == anchor,]
      if (dim(ttable)[1] == 0)
          next
      cat(sprintf("plotting anchor scatters: %d\n", anchor))

      ofn = paste(scatter.fdir, "/", anchor, "_contig_total.png", sep="")
      pscatter(ttable=ttable, ofn=ofn, title="contig total count", field.o="contig_total_count", field.e="contig_expected")

      ofn = paste(scatter.fdir, "/", anchor, "_gene.png", sep="")
      pscatter(ttable=ttable, ofn=ofn, title="gene count", field.o="gene_count", field.e="gene_expected")

      ttable = table[table$anchor == anchor,]
      ofn = paste(scatter.fdir, "/", anchor, "_total_vs_anchor_contigs.png", sep="")
      png(ofn, width=600, height=600)
      ttable = ttable[ttable$contig_total_count > 0 & ttable$anchor_contig_count > 0,]
      plot(ttable$contig_total_count, ttable$anchor_contig_count, pch=".", type="p", log="xy", cex=4, xlim=xlim, ylim=xlim,
           xlab="$contacts", ylab="#contigs", col=ttable$col, main=paste(anchor))
      dev.off()

      ttable = table[table$anchor == anchor & table$contig_total_count > 0 & table$anchor_contig_count > 0 ,]
      ofn = paste(scatter.fdir, "/", anchor, "_contig_total_vs_coverage.png", sep="")
      png(ofn, width=600, height=600)
      plot(ttable$contig_total_count, ttable$contig_coverage, pch=".", type="p", log="x", cex=4, xlim=xlim, ylim=c(0,1),
           xlab="$contacts", ylab="%", col=ttable$col, main=paste(anchor))
      dev.off()
  }
}

anchor.coverage=function(ifn, ifn.cov, fdir)
{
  cat(sprintf("reading table: %s\n", ifn))
  table = read.delim(ifn)

  cat(sprintf("reading table: %s\n", ifn.cov))
  table.cov = read.delim(ifn.cov)

  table$median.cov = table.cov$median[match(table$gene_contig, table.cov$contig)]
  table$enrichment = ifelse(table$gene.selected, table$eng, table$enc)

  system(paste("mkdir -p", fdir))

  ofn = paste(fdir, "/enrichment_vs_abundance_all.png", sep="")
  cat(sprintf("plotting figure: %s\n", ofn))
  png(ofn, width=600, height=600)
  plot(table$enrichment, table$median.cov, col=table$anchor, type="p", pch="+", log="y")
  dev.off()

  en = sapply(split(table$enrichment, table$anchor), median)
  ab = sapply(split(table$median.cov, table$anchor), median)
  col = as.numeric(names(split(table$enrichment, table$anchor)))
  ofn = paste(fdir, "/enrichment_vs_abundance_anchors.png", sep="")
  cat(sprintf("plotting figure: %s\n", ofn))
  png(ofn, width=600, height=600)
  plot(en, ab, col=col, type="p", pch="+", log="y")
  dev.off()
}

field.count=function(x, field="gene")
{
    tt = table(x[,field])
    result = data.frame(x=names(tt), count=as.vector(tt))
    names(result)[1] = field
    result[order(result$count, decreasing=T),]
}

anchor.info=function(ga.ifn, gene.ifn, contig.ifn, gc.ifn, coverage.ifn,
    ofn.summary, ofn.genes, ofn.contigs)
{
    options(stringsAsFactors=F)

    ga = read.delim(ga.ifn)
    contigs = read.delim(contig.ifn)
    genes = read.delim(gene.ifn)
    gc = read.delim(gc.ifn)
    coverage = read.delim(coverage.ifn)

    df = field.count(ga, "gene")
    ga$count = df$count[match(ga$gene, df$gene)]
    ga$enrichment = ifelse(ga$support == "gene", ga$eng, ga$enc)

    gene.fields = c("anchor", "enrichment", "count")
    s = split(ga, ga$anchor)
    genes.info = data.frame(
        anchor=as.numeric(names(s)),
        enrichment=sapply(s, function(x) median(x$enrichment)),
        gene.count=sapply(s, function(x) dim(x)[1]),
        shared.genes=sapply(s, function(x) sum(x$count > 1)))

    cat(sprintf("saving result table: %s\n", ofn.genes))
    write.table(ga[,gene.fields], ofn.genes, quote=F, col.names=T, row.names=F, sep="\t")

    contigs$anchor = ga$anchor[match(contigs$contig, ga$gene_contig)]
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

    cat(sprintf("saving result table: %s\n", ofn.contigs))
    contig.fields = c("anchor", "gc", "reads", "length")
    write.table(contigs[,contig.fields], ofn.contigs, quote=F, col.names=T, row.names=F, sep="\t")

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

    pp=function(field, table, ofn, sort=F) {
        values = split(table[,field], table$anchor)
        if (sort)
            values = values[order(sapply(values, median))]
        cat(sprintf("plotting figure: %s\n", ofn))
        png(ofn, width=800, height=600)
        boxplot(values, main=field, outline=F, las=2)
        dev.off()
    }
    system(paste("mkdir -p", fdir))

    pp(table=genes, field="enrichment", ofn = paste(fdir,"/full_contact_enrichment.png", sep=""))
    pp(table=contigs, field="reads.per.bp", ofn = paste(fdir,"/full_coverage.png", sep=""))
    pp(table=contigs, field="length", ofn = paste(fdir,"/full_contig_length.png", sep=""))
    pp(table=contigs, field="gc", ofn = paste(fdir,"/full_gc.png", sep=""))

    pp(table=genes, field="enrichment", ofn = paste(fdir,"/full_sort_contact_enrichment.png", sep=""), sort=T)
    pp(table=contigs, field="reads.per.bp", ofn = paste(fdir,"/full_sort_coverage.png", sep=""), sort=T)
    pp(table=contigs, field="length", ofn = paste(fdir,"/full_sort_contig_length.png", sep=""), sort=T)
    pp(table=contigs, field="gc", ofn = paste(fdir,"/full_sort_gc.png", sep=""), sort=T)
}

anchor.matrix=function(ga.ifn, map.ifn, ofn)
{
    options(stringsAsFactors=F)

    cat(sprintf("reading ga file: %s\n", ga.ifn))
    ga = read.delim(ga.ifn)[,c("cgene", "anchor")]
    ga = unique(ga)

    cat(sprintf("reading map file: %s\n", map.ifn))
    map = read.delim(map.ifn)

    s = split(ga$cgene, ga$anchor)
    N = length(s)
    result = NULL
    for (i1 in 1:N) {
        for (i2 in 1:N) {
            if (i1 == i2)
                next
            anchor1 = names(s)[i1]
            anchor2 = names(s)[i2]
            cgenes1 = s[[i1]]
            cgenes2 = s[[i2]]

            n1 = length(cgenes1)
            n2 = length(cgenes2)

            # shared cgenes
            shared = unique(intersect(cgenes1, cgenes2))
            n.shared = length(shared)

            cgenes1 = setdiff(cgenes1, shared)
            cgenes2 = setdiff(cgenes2, shared)
            tmap = map[is.element(map$cgene1, cgenes1) & is.element(map$cgene2, cgenes2),]

            # unmapped cgenes
            n.unmapped1 = length(cgenes1) - length(unique(tmap$cgene1))
            n.unmapped2 = length(cgenes2) - length(unique(tmap$cgene2))

            n.strong = sum(tmap$identity >= 70 & tmap$coverage >= 0.7)
            n.weak = sum(tmap$identity < 70 & tmap$coverage < 0.7)

            # compute average identity score
            tmap = tmap[tmap$identity >= 60 & tmap$coverage >= 0.7,]

            identity = mean(c(tmap$identity, rep(100, n.shared)))

            df = data.frame(anchor1=anchor1, anchor2=anchor2, identity=identity, shared=n.shared, strong=n.strong, weak=n.weak)
            result = rbind(result, df)
        }
    }

    cat(sprintf("saving result table: %s\n", ofn))
    write.table(result, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}
