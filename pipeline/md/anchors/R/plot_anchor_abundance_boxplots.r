plot.box.i=function(table, field.group, field.value, groups, col="blue", add.axis=T, mai, add.grid=F, label, plot.top95.line=T, axis.side=2)
{
    fields = c("bottom05", "bottom25", "median", "top75", "top95")

    s = split(table[,field.value], table[,field.group])
    x = as.data.frame(t(as.matrix(sapply(s, function(x) { quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)) }))))
    colnames(x) = fields
    x$group = rownames(x)
    x$coord = match(x$group, groups) - 0.5
    N = length(groups)
    coords = 0:N

    xlim = c(0,N)
    ylim = range(c(x$bottom05, 1.1*x$top95))

    par(mai=mai)
    plot.new()
    plot.window(xlim=xlim, ylim=ylim)
    grid()
    title(main=label)
    if (add.grid) grid(col="darkgray", lty=1)
    if (plot.top95.line)  segments(x0=x$coord, x1=x$coord, y0=x$bottom05, y1=x$top95)
    rect(xleft=x$coord-0.25, ybottom=x$bottom25, xright=x$coord+0.25, ytop=x$top75, col=col, border=col)
    segments(x0=x$coord-0.25, x1=x$coord+0.25, y0=x$median, y1=x$median, lwd=1)
    if (add.axis)
        axis(axis.side, las=2, cex.axis=0.75)
    box()
}

plot.anchor.abundance.boxplots=function(ifn.coverage, ifn, fdir)
{
    cov = load.table(ifn.coverage)
    table = load.table(ifn)
    table$abun = cov$abundance.enrichment[match(table$contig,cov$contig)]
    field.group = "anchor"
    field.value = "abun"
    groups = sort(unique(table$anchor))
    mai = c(0.1,0.5,0.5,0.1)

    fig.start(fdir=fdir, ofn=paste(fdir, "/anchor_abundance.pdf", sep=""), type="pdf", width=6, height=3)
    plot.box.i(table=table, field.group=field.group, field.value=field.value, groups=groups, col="orange", add.axis=T, mai=mai, add.grid=F,
               label="abundance", plot.top95.line=T, axis.side=2)
    # boxplot(s, outline=F, las=2, col="darkblue", border=NA)
    fig.end()
}
