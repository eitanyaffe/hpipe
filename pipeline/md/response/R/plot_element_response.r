plot.elements.internal=function(labels, fdir, lid, rep.table, element.table, anchor.table, q50, norm, min.score, disturb.ids, network, compare)
{
    # remove cluster size
    q50 = q50[,-2]
    element.ids = element.table$id

    M = dim(q50)[2]-1
    coords = 1:M
    ylim = c(1.1*min.score, 1.1*log10(max(q50[,-1])))
    system(paste("rm -rf", fdir))

    sum.color = "darkgreen"
    sum.color1 = "blue"
    sum.color2 = "red"

    dindex = which(is.element(colnames(q50)[-1], disturb.ids))
    drange = c(min(dindex), max(dindex))

    get.anchor.response=function(anchors) {
        rx = log10(q50[match(anchors, q50$cluster), -1])
        rx[rx<min.score] = min.score
        rx
    }
    get.anchor.response.sum=function(anchors) {
        ar = get.anchor.response(anchors)
        rx = log10(colSums(ifelse(ar>min.score,10^ar,0)))
        rx[rx<min.score] = min.score
        rx
    }

    get.response=function(element.id) {
        result = list()
        result$element.id = element.id
        element = element.table$cluster[match(element.id, element.table$id)]
        element.q50 = log10(unlist(q50[element == q50$cluster, -1]))
        element.q50[element.q50<min.score] = min.score
        result$element.response = element.q50

        inetwork = network[network$cluster == element,]
        K = dim(inetwork)[1]
        if (K == 0)
            return (result)
        result$K = K
        anchors = inetwork$anchor
        anchors = anchors[order(match(anchors,anchor.table$set))]
        anchor.response = get.anchor.response(anchors)
        result$anchor.response = anchor.response
        result$anchors = anchors
        result$anchor.desc = paste(rep.table$anchor.id[match(anchors, rep.table$anchor)], ": ",
            rep.table$name[match(anchors, rep.table$anchor)], sep="")

        if (compare) {
            result$anchor.desc = paste(result$anchor.desc,
                ifelse(inetwork$connected1 & inetwork$connected2, "stable", ifelse(inetwork$connected1, "lost", "gained")))
            result$any.connected1 = any(inetwork$connected1)
            result$any.connected2 = any(inetwork$connected2)
            if (result$any.connected1) {
                result$anchors1 = inetwork$anchor[inetwork$connected1]
                result$rsum1 = get.anchor.response.sum(result$anchors1)
                result$delta1 = result$element.response - result$rsum1
                if (element.id == 51)
                    result$rsum1 = get.anchor.response.sum(make.anchor(c("a50", "a57", "a62"), anchor.table))
            }
            if (result$any.connected2) {
                result$anchors2 = inetwork$anchor[inetwork$connected2]
                result$rsum2 = get.anchor.response.sum(result$anchors2)
                if (element.id == 51)
                    result$rsum2 = get.anchor.response.sum(make.anchor(c("a50", "a55", "a57", "a59"), anchor.table))
                result$delta2 = result$element.response - result$rsum2
            }
        } else {
            result$rsum = get.anchor.response.sum(anchors)
            result$delta = result$element.response - result$rsum
        }

        return (result)
    }

    plot.element.legend=function(df, multi, detailed, anchor.colors, color.single, compare, delta) {
        if (!multi)
            axis(side=1, labels=labels, at=1:M, las=2, cex.axis=0.75)
        else
            axis(side=1, labels=F)
        axis(side=2, las=2, labels=!multi)
        leg.colors = c("black")
        leg.names = paste("element", df$element.id)
        if (detailed) {
            leg.colors = c(leg.colors, anchor.colors)
            leg.names = c(leg.names, df$anchor.desc)
        }
        if (compare) {
            leg.colors = c(leg.colors, sum.color1, sum.color2)
            leg.names = c(leg.names, "cummulative1", "cummulative2")
        } else {
            leg.colors = c(leg.colors, sum.color)
            leg.names = c(leg.names, "cummulative")
        }
        if (!multi) {
            par(mai=c(0,0,0,0))
            par(xaxs="i")
            plot.new()
            plot.window(0:1, 0:1)
            legend(x=0, y=0.95, fill=leg.colors, legend=leg.names, border=NA, bty="n")
        }
    }
    plot.element=function(element.id, multi, detailed) {
        df = get.response(element.id)
        element = element.table$cluster[match(element.id, element.table$id)]
        if (df$K > 0)
            colors = colorpanel(df$K, "orange", "purple")
        main = paste(element.id, sep="")
        plot.init(xlim=c(1,M), ylim=ylim, main=main, x.axis=F, y.axis=F)
        abline(h=0)
        at = drange
        abline(v=at, lwd=1, lty=2)
        if (compare && df$K>0) {
            if (df$any.connected1)
                lines(x=coords, y=df$rsum1, col="lightblue", lwd=4)
            if (df$any.connected1) {
                if (!setequal(df$anchors1, df$anchors2))
                    lines(x=coords, y=df$rsum2, col="pink", lwd=4)
            }
        } else {
            if (df$K > 0)
                lines(x=coords, y=df$rsum, col=sum.color, lwd=4)
        }
        lines(x=coords, y=df$element.response, lwd=2)
        if (!multi && df$K > 0 && detailed)
            for (i in 1:df$K)
                lines(x=coords, y=df$anchor.response[i,], col=colors[i], lwd=1)
        plot.element.legend(df=df, multi=multi, detailed=detailed, anchor.colors=colors, compare=compare, delta=F)
    }
    plot.element.delta=function(element.id, multi, detailed, ylim.delta) {
        df = get.response(element.id)
        element = element.table$cluster[match(element.id, element.table$id)]
        if (df$K > 0)
            colors = colorpanel(df$K, "orange", "purple")
        main = paste(element.id, sep="")
        plot.init(xlim=c(1,M), ylim.delta, main=main, x.axis=F, y.axis=F)
        abline(h=0)
        at = drange
        abline(v=at, lwd=1, lty=2)
        if (compare && df$K>0) {
            if (df$any.connected1)
                lines(x=coords, y=df$delta1, col=sum.color1, lwd=2)
            if (df$any.connected1) {
                if (!setequal(df$anchors1, df$anchors2))
                    lines(x=coords, y=df$delta2, col=sum.color2, lwd=2)
            }
        } else {
            if (df$K > 0)
                lines(x=coords, y=df$delta, col="darkgreen", lwd=2)
        }
        plot.element.legend(df=df, multi=multi, detailed=detailed, anchor.colors=colors, compare=compare, delta=T)
    }
    plot.trend=function(trend, title ,color) {
        xtrend = log2(10^trend)
        bp = boxplot(xtrend, outline=F, plot=F)
        ylim = range(bp$stats)
        plot.init(xlim=c(1,M), ylim=ylim, main=title, x.axis=F, y.axis=F, ylab="relative copy, log2", add.grid=F)
        rect(drange[1]-0.5, ylim[1]-1, drange[2]+0.5, ylim[2]+1, col="lightgreen", border=NA)
        grid()
        abline(h=0)
        axis(side=1, labels=labels, at=1:M, las=2, cex.axis=0.75)
        axis(side=2, las=2)
        boxplot(xtrend, col="gray", outline=F, add=T, names=F, axes=F, lty=1)
        box()
    }


    ylim.delta = NULL
    if (compare) {
        trend1 = NULL
        trend2 = NULL
    } else {
        trend = NULL
    }
    for (element.id in element.ids) {
        df = get.response(element.id)
        if (compare) {
            ylim.delta = range(c(ylim.delta, df$delta1, df$delta2))
            trend1 = rbind(trend1, df$delta1 - df$delta1[1])
            trend2 = rbind(trend2, df$delta2 - df$delta2[1])
        } else {
            ylim.delta = range(c(ylim.delta, df$delta))
            trend = rbind(trend, df$delta - df$delta[1])
        }
    }

    # trends
    trend.height = 4
    trend.width = 4
    if (compare) {
        fig.start(fdir=fdir, ofn=paste(fdir, "/trend1.pdf", sep=""), type="pdf", height=trend.height, width=trend.width)
        plot.trend(trend1, "trend1")
        fig.end()
        fig.start(fdir=fdir, ofn=paste(fdir, "/trend2.pdf", sep=""), type="pdf", height=trend.height, width=trend.width)
        plot.trend(trend2, "trend2")
        fig.end()
    } else {
        fig.start(fdir=fdir, ofn=paste(fdir, "/trend.pdf", sep=""), type="pdf", height=trend.height, width=trend.width)
        plot.trend(trend, "trend")
        fig.end()
    }

    # fig per element, detailed
    ffdir = paste(fdir, "/elements_detailed", sep="")
    for (element.id in element.ids) {
        fig.start(fdir=ffdir, ofn=paste(ffdir, "/", element.id, ".pdf", sep=""), type="pdf", height=3, width=8)
        layout(matrix(1:2, 1, 2), widths=c(1,1.3))
        par(mai=c(0.6,0.6,0.4,0))
        plot.element(element.id=element.id, multi=F, detailed=T)
        fig.end()
    }

    # fig per element
    ffdir = paste(fdir, "/elements", sep="")
    for (element.id in element.ids) {
        fig.start(fdir=ffdir, ofn=paste(ffdir, "/", element.id, ".pdf", sep=""), type="pdf", height=3, width=8)
        layout(matrix(1:2, 1, 2), widths=c(1,1.3))
        par(mai=c(0.6,0.6,0.4,0))
        plot.element(element.id=element.id, multi=F, detailed=F)
        fig.end()
    }

    # fig per element, delta
    ffdir = paste(fdir, "/elements_delta", sep="")
    for (element.id in element.ids) {
        fig.start(fdir=ffdir, ofn=paste(ffdir, "/", element.id, ".pdf", sep=""), type="pdf", height=3, width=8)
        layout(matrix(1:2, 1, 2), widths=c(1,1.3))
        par(mai=c(0.6,0.6,0.4,0))
        plot.element.delta(element.id=element.id, multi=F, detailed=F, ylim.delta=ylim.delta)
        fig.end()
    }

    N = length(element.ids)
    Ny = ceiling(sqrt(N))
    Nx = ceiling(N/Ny)
    NN = c(1:N,rep(N+1,Nx*Ny-N))

    # all
    fig.start(fdir=fdir, ofn=paste(fdir, "/all_elements.pdf", sep=""), type="pdf", height=0.5+Ny*1.1, width=0.5+Nx*1.1)
    layout(matrix(NN, Nx, Ny, byrow=T))
    par(mai=c(0.1, 0.1, 0.15, 0.05))
    for (element.id in element.ids)
        plot.element(element.id=element.id, multi=T, detailed=F)
    fig.end()

    # all, delta
    fig.start(fdir=fdir, ofn=paste(fdir, "/all_elements_delta.pdf", sep=""), type="pdf", height=0.5+Ny*1.1, width=0.5+Nx*1.1)
    layout(matrix(NN, Nx, Ny, byrow=T))
    par(mai=c(0.1, 0.1, 0.15, 0.05))
    for (element.id in element.ids)
        plot.element.delta(element.id=element.id, multi=T, detailed=F, ylim.delta=ylim.delta)
    fig.end()
}

plot.elements=function(lid, ifn.reps, ifn.anchors, ifn.elements, ifn.map, ifn.bottom0, ifn.median, ifn.top100,
    ifn.detection, labels, disturb.ids, fdir)
{
    rep.table = load.table(ifn.reps)
    element.table = load.table(ifn.elements)
    anchor.table = load.table(ifn.anchors)
    q50 = load.table(ifn.median)
    min.score = log10(load.table(ifn.detection)[1,1])
    network = load.table(ifn.map)
    network = network[network$type == "connected",]

    plot.elements.internal(labels=labels, fdir=fdir, lid=lid, rep.table=rep.table, disturb.ids=disturb.ids, compare=F,
                           element.table=element.table, anchor.table=anchor.table, q50=q50, min.score=min.score, network=network)
}

plot.elements.compare=function(lid, ifn.reps, ifn.anchors, ifn.elements, ifn.map, ifn.bottom0, ifn.median, ifn.top100,
    ifn.detection, labels, disturb.ids, fdir)
{
    system(paste("rm -rf", fdir))
    rep.table = load.table(ifn.reps)
    element.table = load.table(ifn.elements)
    anchor.table = load.table(ifn.anchors)
    q50 = load.table(ifn.median)
    min.score = log10(load.table(ifn.detection)[1,1])
    network = load.table(ifn.map)

    # limit to changes
    clusters = unique(network$cluster[network$type != "unknown"& network$type != "stable"])
    element.table = element.table[is.element(element.table$cluster, clusters),]

    plot.elements.internal(labels=labels, fdir=fdir, lid=lid, rep.table=rep.table, disturb.ids=disturb.ids, compare=T,
                           element.table=element.table, anchor.table=anchor.table, q50=q50, min.score=min.score, network=network)
}
