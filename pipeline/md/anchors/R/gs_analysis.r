bin.gs.by.aai=function(ifn.prefix, min.k, max.k, aai.breaks, odir)
{
    for (k in min.k:max.k) {
        ifn = paste(ifn.prefix, ".", k, sep="")
        df = load.table(ifn)
        s.min = split(df$gene_count, cut(df$min_aai,breaks=aai.breaks))
        s.max = split(df$gene_count, cut(df$max_aai,breaks=aai.breaks))

        process=function(s, oprefix) {
            names = names(s)
            total = sapply(s, length)
            observed = sapply(s, function(x) { sum(x>0)} )
            Pshared = sapply(s, function(x) {
                if (length(x)>0)
                    round(sum(x>0) / length(x),5)
                else
                    0
            } )
            counts = sapply(s, function(x) { x[x>0] }, simplify=F)
            median.count = sapply(counts, median)

            median.count.result = NULL
            for (i in 1:length(counts)) {
                if (any(counts[[i]] > 0))
                    median.count.result = rbind(median.count.result, data.frame(bin=i, count=counts[[i]]))
            }

            save.table(x=as.data.frame(t(observed)), ofn=paste(oprefix, ".observed", sep=""))
            save.table(x=as.data.frame(t(total)), ofn=paste(oprefix, ".total", sep=""))
            save.table(x=as.data.frame(t(Pshared)), ofn=paste(oprefix, ".Pshared", sep=""))
            save.table(x=as.data.frame(t(median.count)), ofn=paste(oprefix, ".median_count", sep=""))
            save.table(x=median.count.result, ofn=paste(oprefix, ".detailed_count", sep=""))
        }

        if (k == 2) {
            process(s.min, oprefix=paste(odir, "/k", k, ".aai", sep=""))
        } else {
            process(s.min, oprefix=paste(odir, "/k", k, ".min_aai", sep=""))
            process(s.max, oprefix=paste(odir, "/k", k, ".max_aai", sep=""))
        }
    }

}

plot.gs=function(ref.dir, anchor.dir, min.k, max.k, aai.breaks, fdir)
{
    labels = paste(aai.breaks[-length(aai.breaks)], "-", aai.breaks[-1], sep="")
    N = length(labels)
    col = "darkgreen"
    col1 = "darkgreen"
    col2 = "gray"
    min.logp = -6

    make.ifn=function(idir, k, aai.type, ext) {
        sprintf("%s/k%d.%s.%s", idir, k, aai.type, ext)
    }
    make.odir=function(k, aai.type, title) {
        sprintf("%s/%s/k%d_%s", fdir, title, k, aai.type)
    }
    make.ofn=function(odir, ext, add.labels) {
        sprintf("%s/%s_%s.pdf", odir, ext, if (add.labels) "labels" else "clean")
    }
    get.values=function(ifn) {
        df = load.table(ifn)
        vals = unlist(df)
        vals[is.na(vals)] = 0
        vals
    }
    get.binned.values=function(ifn) {
        df = load.table(ifn)
        df$label = labels[df$bin]
        df$f = factor(df$label, levels=labels)
        s = split(df$count, df$f)
        m = boxplot(s, plot=F)
        list(split=s, max=1.1*max(m$stats, na.rm=T))
    }
    mlog=function(vals) {
        ifelse(vals>0, log10(vals), min.logp)
    }
    get.p.values=function(idir, k, aai.type, ext.t, ext.o)
    {
        vals.o = get.values(make.ifn(idir=idir, k=k, aai.type=aai.type, ext=ext.o))
        vals.t = get.values(make.ifn(idir=idir, k=k, aai.type=aai.type, ext=ext.t))
        vals.p = ifelse(vals.t>0, vals.o/vals.t, 0)

        vals.o.high = vals.o + sqrt(vals.o)
        vals.o.low = vals.o - sqrt(vals.o)
        vals.o.low[vals.o.low<0] = 0
        vals.p.high = ifelse(vals.t>0, vals.o.high/vals.t, 0)
        vals.p.low = ifelse(vals.t>0, vals.o.low/vals.t, 0)
        vals.p.high[vals.p.high>1] = 1
        list(vals=mlog(vals.p), low=mlog(vals.p.low), high=mlog(vals.p.high))
    }

    # barplots for real values
    plot.bars=function(title, idir, k, aai.type, ext, add.labels=T) {
        odir = make.odir(k=k, title=title, aai.type=aai.type)
        ofn = make.ofn(odir=odir, ext=ext, add.labels=add.labels)
        vals = get.values(make.ifn(idir=idir, k=k, aai.type=aai.type, ext=ext))
        ylim = c(0, 1.1*max(vals))
        fig.start(ofn=ofn, fdir=odir, width=1+N*0.3, height=4, type="pdf")
        m = barplot(vals, names.arg=labels, las=2, border=NA, ylim=ylim)
        if (add.labels)
            text(x=m, y=vals, pos=3, labels=vals, cex=0.75)
        title(main=paste(ext, k, aai.type), cex=0.75)
        fig.end()
    }
   plot.bars.compare=function(title, idir1, idir2, k, aai.type, ext, add.labels=T) {
        odir = make.odir(k=k, title=title, aai.type=aai.type)
        ofn = make.ofn(odir=odir, ext=ext, add.labels=add.labels)
        vals1 = get.values(make.ifn(idir=idir1, k=k, aai.type=aai.type, ext=ext))
        vals2 = get.values(make.ifn(idir=idir2, k=k, aai.type=aai.type, ext=ext))
        ylim = c(0, 1.1*max(c(vals1,vals2)))
        df = rbind(vals1, vals2)
        fig.start(ofn=ofn, fdir=odir, width=1+N*0.6, height=4, type="pdf")
        m = barplot(df, beside=T, col=c("darkgreen", "gray"), names.arg=labels, las=2, border=NA, ylim=ylim)
        if (add.labels) {
            text(x=m[1,], y=vals1, pos=3, labels=vals1, cex=0.75)
            text(x=m[2,], y=vals2, pos=3, labels=vals2, cex=0.75)
        }
        title(main=paste(ext, k, aai.type), cex=0.75)
        fig.end()
    }

    # P values
    get.plabel=function(vals) {
        ifelse(vals == min.logp, "NA", round(vals, 1))
    }
    plot.bars.p=function(title, idir, k, aai.type, ext.t, ext.o, add.labels=T) {
        odir = make.odir(k=k, title=title, aai.type=aai.type)
        ofn = make.ofn(odir=odir, ext="P", add.labels=add.labels)
        pvals = get.p.values(idir=idir, k=k, aai.type=aai.type, ext.t=ext.t, ext.o=ext.o)

        fig.start(ofn=ofn, fdir=odir, width=1+N*0.3, height=4, type="pdf")

        ylim = extendrange(c(pvals$low, pvals$high), f=0.1)
        coord = 1:N-0.5
        plot.init(xlim=c(0, N), ylim=ylim, ylab="P", x.axis=F, axis.las=2)
        rect(xleft=coord-0.25, xright=coord+0.25, ybottom=pvals$low, ytop=pvals$high, col=col, border=col)
        segments(x0=coord-0.25, x1=coord+0.25, y0=pvals$vals, y1=pvals$vals)
        axis(1, at=coord, labels=labels, lty=0, las=2)
        if (add.labels)
            text(x=coord, y=pvals$high, pos=3, labels=get.plabel(pvals$vals), cex=0.75)
        title(main=paste("P", k, aai.type), cex=0.75)

        fig.end()
    }
    plot.bars.p.compare=function(title, idir1, idir2, k, aai.type, ext.t, ext.o, add.labels=T) {
        odir = make.odir(k=k, title=title, aai.type=aai.type)
        ofn = make.ofn(odir=odir, ext="P", add.labels=add.labels)
        pvals1 = get.p.values(idir=idir1, k=k, aai.type=aai.type, ext.t=ext.t, ext.o=ext.o)
        pvals2 = get.p.values(idir=idir2, k=k, aai.type=aai.type, ext.t=ext.t, ext.o=ext.o)


        fig.start(ofn=ofn, fdir=odir, width=1+N*0.6, height=4, type="pdf")

        ylim = extendrange(c(pvals1$low, pvals1$high, pvals2$low, pvals2$high), f=0.1)
        coord = 2*(1:N)-1
        coord1 = 2*(1:N)-1.5
        coord2 = 2*(1:N)-0.5

        plot.init(xlim=c(0, N*2), ylim=ylim, ylab="P", x.axis=F, axis.las=2)
        rect(xleft=coord1-0.25, xright=coord1+0.25, ybottom=pvals1$low, ytop=pvals1$high, col=col1, border=col1)
        segments(x0=coord1-0.25, x1=coord1+0.25, y0=pvals1$vals, y1=pvals1$vals)
        rect(xleft=coord2-0.25, xright=coord2+0.25, ybottom=pvals2$low, ytop=pvals2$high, col=col2, border=col2)
        segments(x0=coord2-0.25, x1=coord2+0.25, y0=pvals2$vals, y1=pvals2$vals)

        axis(1, at=coord, labels=labels, lty=0, las=2)
        if (add.labels) {
            text(x=coord1, y=pvals1$high, pos=3, labels=get.plabel(pvals1$vals), cex=0.75)
            text(x=coord2, y=pvals2$high, pos=3, labels=get.plabel(pvals2$vals), cex=0.75)
        }
        title(main=paste("P", k, aai.type), cex=0.75)

        fig.end()
    }

    plot.box=function(title, idir, k, aai.type, ext, add.labels=T) {
        odir = make.odir(k=k, title=title, aai.type=aai.type)
        ofn = make.ofn(odir=odir, ext=ext, add.labels=add.labels)
        vals = get.binned.values(make.ifn(idir=idir, k=k, aai.type=aai.type, ext=ext))
        ylim = c(0, vals$max)
        fig.start(ofn=ofn, fdir=odir, width=1+N*0.5, height=4, type="pdf")
        coord = 1:N
        xlim = c(0.5,N+0.5)
        m = boxplot(vals$split, at=coord, names=labels, las=2, ylim=ylim, outline=F, col=col, xlim=xlim)
        if (add.labels)
            text(x=coord, y=m$stats[2,], pos=1, labels=m$stats[3,], cex=0.75)
        title(main=paste(ext, k, aai.type), cex=0.75)
        fig.end()
    }
    plot.box.compare=function(title, idir1, idir2, k, aai.type, ext, add.labels=T) {
        odir = make.odir(k=k, title=title, aai.type=aai.type)
        ofn = make.ofn(odir=odir, ext=ext, add.labels=add.labels)
        vals1 = get.binned.values(make.ifn(idir=idir1, k=k, aai.type=aai.type, ext=ext))
        vals2 = get.binned.values(make.ifn(idir=idir2, k=k, aai.type=aai.type, ext=ext))

        fig.start(ofn=ofn, fdir=odir, width=1+N*0.5, height=4, type="pdf")
        ylim = c(0, max(c(vals1$max, vals2$max)))
        xlim = c(0.5,2*N+0.5)
        coord = 2*(1:N)-0.5
        coord1 = 2*(1:N)-1
        coord2 = 2*(1:N)
        m1 = boxplot(vals1$split, at=coord1, names=F, las=2, ylim=ylim, outline=F, col=col1, xlim=xlim)
        m2 = boxplot(vals2$split, at=coord2, names=F, las=2, outline=F, col=col2, add=T)
        axis(1, at=coord, labels=labels, lty=0, las=2)
        if (add.labels) {
            text(x=coord1, y=m1$stats[2,], pos=1, labels=m1$stats[3,], cex=0.75)
            text(x=coord2, y=m2$stats[2,], pos=1, labels=m2$stats[3,], cex=0.75)
        }
        title(main=paste(ext, k, aai.type), cex=0.75)
        fig.end()
    }

    for (k in min.k:max.k) {
        aai.types = if (k == 2) "aai" else c("min_aai", "max_aai")
        for (add.labels in c(F,T)) {
        for (aai.type in aai.types) {

            plot.box(title="ref", idir=ref.dir, k=k, aai.type=aai.type, ext="detailed_count", add.labels=add.labels)
            plot.bars(title="ref", idir=ref.dir, k=k, aai.type=aai.type, ext="observed", add.labels=add.labels)
            plot.bars(title="ref", idir=ref.dir, k=k, aai.type=aai.type, ext="median_count", add.labels=add.labels)
            plot.bars.p(title="ref", idir=ref.dir, k=k, aai.type=aai.type, ext.o="observed", ext.t="total",
                        add.labels=add.labels)

            plot.box(title="anchor", idir=anchor.dir, k=k, aai.type=aai.type, ext="detailed_count", add.labels=add.labels)
            plot.bars(title="anchor", idir=anchor.dir, k=k, aai.type=aai.type, ext="observed", add.labels=add.labels)
            plot.bars(title="anchor", idir=anchor.dir, k=k, aai.type=aai.type, ext="median_count", add.labels=add.labels)
            plot.bars.p(title="anchor", idir=anchor.dir, k=k, aai.type=aai.type, ext.o="observed", ext.t="total",
                        add.labels=add.labels)

            plot.box.compare(title="compare", idir1=anchor.dir, idir2=ref.dir, k=k,
                             aai.type=aai.type, ext="detailed_count", add.labels=add.labels)
            plot.bars.compare(title="compare", idir1=anchor.dir, idir2=ref.dir, k=k,
                              aai.type=aai.type, ext="total", add.labels=add.labels)
            plot.bars.compare(title="compare", idir1=anchor.dir, idir2=ref.dir, k=k,
                              aai.type=aai.type, ext="median_count", add.labels=add.labels)
            plot.bars.compare(title="compare", idir1=anchor.dir, idir2=ref.dir, k=k,
                              aai.type=aai.type, ext="observed", add.labels=add.labels)
            plot.bars.p.compare(title="compare", idir1=anchor.dir, idir2=ref.dir, k=k,
                                aai.type=aai.type, ext.o="observed", ext.t="total",
                                add.labels=add.labels)
        } }
    }
}
