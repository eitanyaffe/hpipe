plot.timeline=function(ids, disturb.days, days, zoom.days, labels, fdir)
{
    xlim = c(min(days)-1, max(days)+1)

    N = length(days)
    bin.step = 7
    dbin.step = bin.step*2
    bin.start = floor(days[1]/dbin.step)
    bin.end = ceiling(days[N]/dbin.step)
    bins = bin.start:bin.end * dbin.step
    ylim = c(-1,1.3)

    plot.internal=function(drange, clean) {
        ix = days>=drange[1] & days<=drange[2]
        ddays = days[ix]
        dlabels = labels[ix]

        par(mai=c(0.1,0.3,0.8,0.3))
        par(xaxs="i")
        plot.new()
        plot.window(xlim=drange, ylim=ylim)
        segments(x0=ddays, y0=0.8, x1=ddays, y1=1.3)
        rect(bins, 0, bins+bin.step, 0.8, col="darkgray", border="darkgray")
        rect(bins+bin.step, 0, bins+bin.step*2, 0.8, col="lightgray", border="lightgray")
        rect(xleft=disturb.days[1], ybottom=0, xright=disturb.days[2], ytop=0.8, col="darkgreen", border=NA, density=40, angle=0)
        # abline(h=0)

        if (!clean)
            axis(3, at=ddays, labels=dlabels, las=2, tick=F)
    }

    height = 1.2

    fig.start(fdir=fdir, ofn=paste(fdir, "/clean.pdf", sep=""), type="pdf", width=8, height=height)
    plot.internal(drange=xlim, clean=T)
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/labels.pdf", sep=""), type="pdf", width=8, height=height)
    plot.internal(drange=xlim, clean=F)
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/clean_zoom.pdf", sep=""), type="pdf", width=8, height=height)
    plot.internal(drange=zoom.days, clean=T)
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/labels_zoom.pdf", sep=""), type="pdf", width=8, height=height)
    plot.internal(drange=zoom.days, clean=F)
    fig.end()
}
