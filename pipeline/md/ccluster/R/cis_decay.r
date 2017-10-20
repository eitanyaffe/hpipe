
plot.cis.decay=function(ifn, id, fdir, distal.threshold, ylim.n)
{
    lwd = 2
    x = load.table(ifn)
    x = x[x$potential > 0,]
#    for (f in c("face", "back", "pos", "neg", "inter"))
#        x[,f] = x[,f] + 1

    x$total = x$face + x$back + x$pos + x$neg + x$inter
    x$cs = 1 - cumsum(x$total) / sum(x$total)

    x$face.n = x$face / x$potential
    x$back.n = x$back / x$potential
    x$pos.n = x$pos / x$potential
    x$neg.n = x$neg / x$potential
    ss = sum(c(x$face.n, x$back.n, x$pos.n, x$neg.n))
    x$face.n = x$face.n / ss
    x$back.n = x$back.n / ss
    x$pos.n = x$pos.n / ss
    x$neg.n = x$neg.n / ss

    x = x[x$potential > 0 & x$contig_count > 100,]

#    ylim = range(1, max(c(x$inter, x$face)))
    ylim = range(1, max(x$face))

    fig.start(fdir=fdir, ofn=paste(fdir, "/raw_", id, ".png", sep=""))
    plot(x$log_start, x$face, log="y", type="l", main=paste("raw",id),lwd=lwd, ylim=ylim)
    grid()
    lines(x$log_start, x$back, col=2,lwd=lwd)
    lines(x$log_start, x$pos, col="gray",lwd=lwd)
    lines(x$log_start, x$neg, col="gray",lwd=lwd)
#    lines(x$log_start, x$inter, col="blue",lwd=lwd)
    abline(v=log10(distal.threshold), lty=2)
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/norm_", id, ".png", sep=""))
    plot(x$log_start, x$face.n, log="y", type="l", main=paste("norm",id), lwd=lwd, ylim=ylim.n)
    grid()
    lines(x$log_start, x$back.n, col=2,lwd=lwd)
    lines(x$log_start, x$pos.n, col="gray",lwd=lwd)
    lines(x$log_start, x$neg.n, col="gray",lwd=lwd)
    abline(v=log10(distal.threshold), lty=2)
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/contigs_", id, ".png", sep=""))
    plot(x$log_start, x$contig_count, log="y", type="l", main=paste("contigs",id),lwd=lwd)
    fig.end()

    # distrib
    fig.start(fdir=fdir, ofn=paste(fdir, "/distrib_", id, ".png", sep=""))
    plot.new()
    ylim = c(0.00001, 10)
    xlim = range(x$log_start)
    plot.window(xlim=xlim, ylim=ylim, log="y")
    title(main=paste("distrib",id))
    grid()
    abline(h=1, lty=2)
    lines(x$log_start, x$cs, lwd=lwd)
    for (i in 1:2) axis(i)

    fig.end()

}

plot.summary=function(distal.threshold, adir, mreads, aid, ids, titles, fdir, ymax)
{
    ymax.v = ifelse(ymax==0, NA, ymax)
    table = NULL
    for (i in 1:length(ids)) {
        id = ids[i]
        title = titles[i]

        xname = if(mreads == -1) "decay" else sprintf("cis_decay_S%d", mreads)
        ifn = paste(adir, "/datasets/", id, "/", xname, "/cis_decay", sep="")
        x = load.table(ifn)
        x$total = x$face + x$back + x$pos + x$neg

        ind = 10^x$log_start < distal.threshold
        close = sum(x$total[ind]) / sum(x$potential[ind])
        far = sum(x$total[!ind]) / sum(x$potential[!ind])
        ss = close + far
        table = rbind(table, data.frame(id=id, close=close/ss, far=far/ss, far.count=sum(x$total[!ind])))
    }
    N = dim(table)[1]

    # distrib
    fig.start(fdir=fdir, ofn=paste(fdir, "/", aid, "_summary.png", sep=""), height=400, width=120+N*35)
    par(mai=c(2,1,1,1))
    if (ymax > 0)
        ylim = c(0, ymax)
    else
        ylim = c(0, max(table$far*1.1))
    barplot(table$far, names.arg=table$id, las=2, main=sprintf("%s, P(>%dk)", aid, distal.threshold/1000), ylim=ylim)
    fig.end()

    fig.start(fdir=fdir, ofn=paste(fdir, "/", aid, "_summary_count.png", sep=""), height=400, width=120+N*35)
    par(mai=c(2,1,1,1))
    barplot(table$far.count, names.arg=table$id, las=2, main=sprintf("%s, N(>%dk)", aid, distal.threshold/1000))
    fig.end()

}
