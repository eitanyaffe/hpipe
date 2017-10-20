plot.node.degree=function(ifn.network, ifn.majors, fdir)
{
    majors = load.table(ifn.majors)
    anchors = majors$anchor

    df = load.table(ifn.network)
    df$element = df$cluster
    elements = sort(unique(df$element))

    n.anchors = length(anchors)
    n.elements = length(elements)
    n.pairs = dim(df)[1]

    make.matrix=function(df) {
        df$i = match(df$anchor, anchors)
        df$j = match(df$element, elements)
        df$value = 1
        smatrix2matrix(df, dim=c(n.anchors, n.elements))
    }
    element.degree=function(m) {
        result = as.data.frame(table(colSums(m)))
        names(result) = c("degree", "count")
        result$degree = as.numeric(levels(result[,1]))[as.numeric(result$degree)]
        result
    }
    host.degree=function(m) {
        result = as.data.frame(table(rowSums(m)))
        names(result) = c("degree", "count")
        result$degree = as.numeric(levels(result[,1]))[as.numeric(result$degree)]
        result
    }

    # observed
    mm = make.matrix(df)
    element.table = element.degree(mm)
    host.table = host.degree(mm)

    #fig.start(fdir=fdir, ofn=paste(fdir, "/host_per_element.png", sep=""), height=400, width=180)
    #barplot(element.table$count, names.arg=element.table$degree, ylab="#elements", xlab="#hosts", col="red")
    #fig.end()

    #anchor.marg = rowSums(mm)
    #element.marg = colSums(mm)

    # random
    rnumber=function(n,from,to) { ceiling(runif(n, from-1, to)) }

    random.element = data.frame(degree=0:(max(element.table$degree)), count=0)
    random.host = data.frame(degree=0:(max(host.table$degree)), count=0)
    n.random = 1000
    for (i in 1:n.random) {
        # random keeping a minimum degree of 2 per each element
        rmm = matrix(0, n.anchors, n.elements)
        for (j in 1:n.elements)
            rmm[sample.int(size=2, n=n.anchors),j] = 1
        while(sum(rmm) < n.pairs) {
            nn = n.pairs - sum(rmm)
            rmm[sample.int(size=nn, n=length(rmm))] = 1
        }

        relement.table = element.degree(rmm)
        rhost.table = host.degree(rmm)

        mx = match(random.element$degree, relement.table$degree)
        random.element$count = random.element$count + ifelse(!is.na(mx), relement.table$count[mx], 0)

        mx = match(random.host$degree, rhost.table$degree)
        random.host$count = random.host$count + ifelse(!is.na(mx), rhost.table$count[mx], 0)
    }
    random.element$count = random.element$count/n.random
    random.host$count = random.host$count/n.random

    pplot=function(random.table, table, ofn, min.degree=0, max.degree=0, title) {
        xtable = merge(random.table, table, by="degree", all=T)
        xtable$obs = ifelse(!is.na(xtable$count.y), xtable$count.y, 0)
        xtable$exp = ifelse(!is.na(xtable$count.x), xtable$count.x, 0)
        xtable = xtable[,c("degree", "obs", "exp")]
        xtable = xtable[xtable$degree>=min.degree,]
        if (max.degree != 0) {
            xtable.top = xtable[xtable$degree>=max.degree,]
            xtable = rbind(xtable[xtable$degree<max.degree,],
                data.frame(degree=paste(">=", max.degree, sep=""), obs=sum(xtable.top$obs), exp=sum(xtable.top$exp)))
        }
        print(xtable)
        mbar = t(as.matrix(xtable[,c("obs", "exp")]))
        fig.start(fdir=fdir, ofn=ofn, height=400, width=300)
        mp = barplot(mbar, beside=T, names.arg=xtable$degree, ylab=paste(title, "count"), xlab="degree",
            col=c("red", "gray"), main=paste(title, "degree"), border=NA)
        fig.end()
    }
    pplot(random.table=random.host, table=host.table, ofn=paste(fdir, "/host_degree.png", sep=""), max.degree=5, title="host")
    pplot(random.table=random.element, table=element.table, ofn=paste(fdir, "/element_degree.png", sep=""), max.degree=0, min.degree=2, title="element")
}

# collapse to host network
plot.host2host.degree=function(ifn.network, ifn.majors, fdir)
{
    majors = load.table(ifn.majors)
    df = load.table(ifn.network)

    anchors = majors$anchor
    n.anchors = length(anchors)

    s = split(df$anchor, df$cluster)
    mm = matrix(0, n.anchors, n.anchors)
    for (i in 1:length(s)) {
        gg = expand.grid(s[[i]], s[[i]])
        gg = gg[gg[,1] != gg[,2],]
        for (j in 1:dim(gg)[1])
            mm[match(gg[j,1], anchors), match(gg[j,2], anchors)] = 1
    }
    N = sum(mm)/2

    mdegree=function(m) {
        result = as.data.frame(table(rowSums(m)))
        names(result) = c("degree", "count")
        result$degree = as.numeric(levels(result[,1]))[as.numeric(result$degree)]
        result
    }
    table = mdegree(mm)

    random.table = data.frame(degree=0:(max(table$degree)), count=0)
    n.random = 100
    for (i in 1:n.random) {
        # random keeping a minimum degree of 2 per each element
        rmm = matrix(0, n.anchors, n.anchors)
        while (sum(rmm)/2 < N) {
            nx = N - sum(rmm)/2
            rmm[sample.int(size=nx, n=length(rmm))] = 1
            rmm = rmm + t(rmm)
            rmm[rmm>1] = 1
            diag(rmm) = 0
        }
        rtable = mdegree(rmm)
        mx = match(random.table$degree, rtable$degree)
        random.table$count = random.table$count + ifelse(!is.na(mx), rtable$count[mx], 0)
    }
    random.table$count = random.table$count/n.random

    pplot=function(random.table, table, ofn, min.degree=0, max.degree=0, title) {
        xtable = merge(random.table, table, by="degree", all=T)
        xtable$obs = ifelse(!is.na(xtable$count.y), xtable$count.y, 0)
        xtable$exp = ifelse(!is.na(xtable$count.x), xtable$count.x, 0)
        xtable = xtable[,c("degree", "obs", "exp")]
        xtable = xtable[xtable$degree>=min.degree,]
        if (max.degree != 0) {
            xtable.top = xtable[xtable$degree>=max.degree,]
            xtable = rbind(xtable[xtable$degree<max.degree,],
                data.frame(degree=paste(">=", max.degree, sep=""), obs=sum(xtable.top$obs), exp=sum(xtable.top$exp)))
        }
        print(xtable)
        mbar = t(as.matrix(xtable[,c("obs", "exp")]))
        fig.start(fdir=fdir, ofn=ofn, height=400, width=300)
        barplot(mbar, beside=T, names.arg=xtable$degree, ylab=paste(title, "count"), xlab="degree",
                col=c("red", "gray"), main=paste(title, "degree"), border=NA)
        fig.end()
    }
    pplot(random.table=random.table, table=table, ofn=paste(fdir, "/host2host_degree.png", sep=""), max.degree=5, title="host")
}

plot.identity.vs.sharing=function(ifn.network, ifn.majors, ifn.identity, fdir)
{
    table = load.table(ifn.identity)
    table$anchor1 = table$set1
    table$anchor2 = table$set2
    table = table[,c("anchor1" ,"anchor2", "identity")]
    table$connected = F
    table = table[table$anchor1 != "NONE" & table$anchor2 != "NONE",]

    majors = load.table(ifn.majors)
    anchors = majors$anchor
    n.anchors = length(anchors)

    df = load.table(ifn.network)
    df$element = df$cluster

    compute.matrix=function(df) {
        result = table
        s = split(df$anchor, df$cluster)
        for (i in 1:length(s)) {
            gg = expand.grid(s[[i]], s[[i]])
            gg = gg[gg[,1] != gg[,2],]
            for (j in 1:dim(gg)[1]) {
                ix = table$anchor1 == gg[j,1] & table$anchor2 == gg[j,2]
                if (any(ix))
                    result$connected[ix] = T
            }
        }
        result = result[result$anchor1 < result$anchor2,]
        result
    }

    breaks = c(0, 70, 80, 90, 100)
    tt = compute.matrix(df)
    result = sapply(split(tt$connected, cut(tt$identity, breaks=breaks, include.lowest=T)), sum)

    # random
    result.rnd = NULL
    n.random = 1000
    for (i in 1:n.random) {
        rs = sapply(split(sample(tt$connected), cut(tt$identity, breaks=breaks, include.lowest=T)), sum)
        if (i == 1)
            result.rnd = rs
        else
            result.rnd = result.rnd + rs
    }
    result.rnd = result.rnd / n.random

    fig.start(fdir=fdir, ofn=paste(fdir, "/host_per_element_vs_random.png", sep=""), height=400, width=300)
    barplot(rbind(result, result.rnd), beside=T, names.arg=names(s), las=2, ylab="#anchor pairs", xlab="%identity", col=c("red", "gray"), border=NA,
            main=paste("anchor pairs =", sum(result)))
    fig.end()
}
