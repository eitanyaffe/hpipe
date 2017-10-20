short.levels=function() {
    rank.list = list(
        superkingdom="spK", kingdom="K", subkingdom="sbK",
        phylum="P", subphylum="sbP",
        superclass="spC", class="C", subclass="sbC", infraclass="iC",
        superorder="spO", order="O", suborder="sbO", infraorder="iO", parvorder="pO",
        superfamily="spF", family="F", subfamily="sbF",
        tribe="T", subtribe="sbT",
        subgenus="sbG", genus="G",
        species.group="Sg" , species.subgroup="Ssg", species="S", subspecies="sbS", strain="sbS",
        varietas="V", forma="FO", no.rank="nr", unknown="X")
    return (rank.list)
}

plot.barplot=function(ifn, fdir)
{
    res = load.table(ifn)
    res$flevel = factor(res$short.level, levels=unique(unlist(short.levels())))
    tt = table(res$flevel)
    tt = tt[tt>0]
    df = data.frame(level=names(tt))
    df$count = tt
    fig.start(fdir=fdir, ofn=paste(fdir, "/resolve_barplot.pdf", sep=""), type="pdf", height=6, width=4)
    barplot(height=df$count, names.arg=df$flevel, border=NA, col="darkgray", axes=F, cex.axis=1, las=2, ylim=c(0, 1.1*max(df$count)))
    axis(2, las=2)
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/resolve_barplot_labels.pdf", sep=""), type="pdf", height=6, width=4)
    bb = barplot(height=df$count, names.arg=df$flevel, border=NA, col="darkgray", axes=F, cex.axis=1, las=2, ylim=c(0, 1.1*max(df$count)))
    text(bb, df$count, pos=3, labels=df$count)
    axis(2, las=2)
    fig.end()
}

plot.details=function(df, fdir)
{
    N = dim(df)[1]
    df$ycoord = 1:N-0.5
    size = 0.4

    colors = c("black", "darkblue", "blue", "red", "orange", "yellow")
    breaks = c(0, 85, 90, 95, 98, 100)
    panel = make.color.panel(colors)
    wlegend2(fdir=fdir, panel=panel, breaks=breaks, title="identity")
    df$color = panel[vals.to.cols(df$identity, breaks)]
    df$coord = df$coverage

    fig.start(fdir=fdir, ofn=paste(fdir, "/summary.pdf", sep=""), type="pdf", width=8, height=2 + N*0.16)

    par(mai=c(1, 1, 0.2, 6))
    par(yaxs="i")
    plot.new()
    plot.window(xlim=c(0,100), ylim=c(0,N))
    title(xlab="genes (%)")

    # all genome
    rect(xleft=0, xright=1, ybottom=df$ycoord-size, ytop=df$ycoord+size, border=NA, col="lightgray")

    # identity region
    rect(xleft=0, xright=df$coord, ybottom=df$ycoord-size, ytop=df$ycoord+size, border=NA, col=df$color)

    # axis
    axis(side=2, at=1:N-0.5, labels=df$anchor.id, las=2)
    axis(side=4, at=1:N-0.5, labels=df$title, las=2)

    # coord axis
    axis(side=1)
    fig.end()

    levels = unlist(short.levels())
    levels = levels[is.element(levels, df$short.level)]
    M = length(levels)
    fig.start(fdir=fdir, ofn=paste(fdir, "/taxa_legend.pdf", sep=""), type="pdf", width=5, height=0.1 + M*0.2)

    par(mai=c(0, 0, 0, 0))
    plot.new()
    plot.window(xlim=c(0,1), ylim=c(0,M))
    text(x=0, y=rev(1:M-0.5), pos=4, labels=paste(levels, ":", names(levels)))
    fig.end()
}

plot.details.resolve=function(ifn, fdir)
{
    df = load.table(ifn)
    df$title = paste(df$short.level, " | ", df$name, " | ", df$type)
    plot.details(df=df, fdir=fdir)
}

plot.details.ref=function(ifn, fdir)
{
    df = load.table(ifn)
    df$title = df$name
    plot.details(df=df, fdir=fdir)
}
