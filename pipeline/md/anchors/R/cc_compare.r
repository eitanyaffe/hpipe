plot.cc.compare=function(ifn1, ifn2, fdir)
{
    t1 = load.table(ifn1)
    t2 = load.table(ifn2)

    t1$score1 = log10(t1$contig_total_count / t1$contig_expected)
    t2$score2 = log10(t2$contig_total_count / t2$contig_expected)

    m = merge(t1, t2, by=c("contig", "anchor"), all=T)
    m$score1[is.na(m$score1)] = 0
    m$score2[is.na(m$score2)] = 0
    m = m[m$score1 > 1 | m$score2 > 1,]

    anchors = sort(unique(m$anchor))
    for (anchor in anchors) {
        fig.start(fdir=fdir, ofn=paste(fdir, "/", anchor, ".png", sep=""))
        x = m[m$anchor == anchor,]
        lim = c(0,7)
        plot.init(xlim=lim, ylim=lim, main=anchor)
        points(x$score1, x$score2, pch=".")
        fig.end()
    }
}
