
plot.path=function(ifn.order, ifn.path, ifn.species, fdir)
{
    table = load.table(ifn.path)
    species = load.table(ifn.species)
    anchor.table = load.table(ifn.order)
    table$anchor.index = match(table$anchor, anchor.table$set)
    ids = anchor.table$id
    N = length(ids)

    levels = c("species", "genus", "family", "order", "class", "phylum", "group", "superkingdom")
    ix = match(table$level, levels)
    table$level.index = ifelse(is.na(ix), -1, ix)
    M = length(levels)
    table = table[table$level.index > 0,]
    gap = 0.5

    xlim = c(0, M+0.5)
    ylim = c(0, N+0.5+gap)

    plot.f=function(title, field) {
        fig.start(fdir=fdir, ofn=paste(fdir, "/", field, ".pdf", sep=""), type="pdf", height=3+N*0.2, width=5+M*0.2)

        plot.new()
        par(mai=c(2,0.5,0.5,4))
        par(xaxs="i", yaxs="i")
        plot.window(xlim=xlim, ylim=ylim)
        rect(ybottom=table$anchor.index-0.5, ytop=table$anchor.index+0.5,
             xleft=table$level.index-0.5, xright=table$level.index+0.5, border=1, col=table[,field])
        axis(2, at=1:N, labels=ids, las=2)
        axis(1, at=1:M, labels=levels, las=2)
        axis(4, at=1:N, labels=species$name, las=2)

        fig.end()
    }

    identity.colors = c("white", "blue", "red", "orange", "yellow")
    identity.breaks = c(0, 60, 70, 80, 100)
    panel.identity = make.color.panel(identity.colors, ncols=256)
    wlegend2(fdir=fdir, panel=panel.identity, breaks=identity.breaks, title="identity")
    table$identity.col = panel.identity[vals.to.cols(table$identity, identity.breaks)]
    plot.f(title="identity", field="identity.col")

    frac.colors = c("white", "blue", "red", "orange", "yellow")
    frac.breaks = c(0, 60, 70, 80, 100)
    panel.frac = make.color.panel(frac.colors, ncols=256)
    wlegend2(fdir=fdir, panel=panel.frac, breaks=frac.breaks, title="frac")
    table$frac.col = panel.frac[vals.to.cols(table$frac, frac.breaks)]
    plot.f(title="frac", field="frac.col")
}
