plot.info=function(adir, mreads, fdir, aid, ids, titles, intra.min=30)
{
    library(graphics)
    intra.close.range = NULL
    inter.range = NULL

    ll = list()
    for (i in 1:length(ids)) {
        id = ids[i]
        title = titles[i]

        xname = if(mreads == -1) "ccontig_analysis" else sprintf("ccontig_analysis_S%d", mreads)
        ifn = paste(adir, "/datasets/", id, "/", xname, sep="")
        x = load.table(ifn)
        ll[[title]] = x
        intra.close.range = range(c(intra.close.range, x$intra_close_count))
        inter.range = range(c(inter.range, x$inter_count))
    }

    cols = rainbow(length(ids))

    # inter / intra percentage
    result = NULL
    for (i in 1:length(ids)) {
        title = titles[i]
        x = ll[[title]]
        r = 100 * sum(x$inter_count) / sum(x$intra_close_count)
        r2 = 100 * sum(x$intra_far_count) / sum(x$intra_close_count)
        result = rbind(result, data.frame(title=title, inter=r, far.intra=r2))
    }
    fig.start(fdir=fdir, ofn=paste(fdir, "/", aid, "_enrichment.png", sep=""), width=800, height=400)
    layout(matrix(1:2, 1, 2))
    barplot(result$inter, names.arg=result$title, col=cols, main="inter")
    barplot(result$far.intra, names.arg=result$title, col=cols, main="intra >2k")
    fig.end()

    # scatters
    fdir = paste(fdir, "/", aid, sep="")

    xlim=c(1, max(intra.close.range))
    ylim=c(1, max(inter.range))
    for (i in 1:length(ids)) {
        title = titles[i]
        x = ll[[title]]
        x = x[x$intra_close_count > 0 & x$inter_count > 0,]

        # points
        fig.start(fdir=fdir, ofn=paste(fdir, "/scatter_", title, ".png", sep=""))
        plot.new()
        plot.window(log="xy", xlim=xlim, ylim=ylim)
        title(xlab="#intra", ylab="#inter", main=paste(title, ", intra_n=", sum(x$intra_close_count), sep=""))
        axis(side=1)
        axis(side=2)
        grid()
        points(x$intra_close_count, x$inter_count, pch=".")
        fig.end()

        # smooth
        fig.start(fdir=fdir, ofn=paste(fdir, "/smooth_scatter_", title, ".png", sep=""))
        smoothScatter(log10(x$intra_close_count), log10(x$inter_count), xlim=log10(xlim), ylim=log10(ylim), main=title)
        fig.end()
    }
}
