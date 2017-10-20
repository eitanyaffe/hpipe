
cluster.sets=function(mat.ifn, prefix, min.identity, ofn)
{
    smat = read.delim(mat.ifn)
    smat = smat[smat$set1 != "NONE" & smat$set2 != "NONE",]

    # remove character R from reference ids
#    if (!is.numeric(smat$set1)) {
#        substr(smat$set1,start=0,stop=1) = "0"
#        substr(smat$set2,start=0,stop=1) = "0"
#        smat$set1 = as.numeric(smat$set1)
#        smat$set2 = as.numeric(smat$set2)
#    }

#    sets = sort(as.numeric(unique(smat$set1)))

    sets = unique(smat$set1)

    N = length(sets)
    smat$set1.i = match(smat$set1, sets)
    smat$set2.i = match(smat$set2, sets)
    m = smatrix2matrix(smat, dim=c(N,N), i.field="set1.i", j.field="set2.i", value.field="identity", default.value=min.identity)
    m[m<min.identity] = min.identity
    d = max(m) - as.dist(m)
    hc = hclust(d, method="average")

    #dd = as.matrix(d)
    #o = hc$order
    #image(dd[o,o])

    df = data.frame(id=paste(prefix, 1:N, sep=""), set=sets[hc$order])
    cat(sprintf("saving result: %s\n", ofn))
    write.table(df, ofn, sep="\t", quote=F, col.names=T, row.names=F)
}

distrib.top.hit=function(
    script, self, idir, odir, qsub.dir, batch.max.jobs, total.max.jobs.fn, dtype, jobname)
{
    source("R/distrib.r")

    cat(sprintf("searching for source dirs in: %s\n", idir))
    src.dirs = list.files(idir)
    N = length(src.dirs)
    if (N == 0)
        stop("no files found")
    cat(sprintf("processing %s source directories found\n", N))

    cat(sprintf("cleaning tmp dir: %s\n", qsub.dir))
    system(paste("rm -rf", qsub.dir))
    system(paste("mkdir -p", qsub.dir))

    for (i  in 1:N) {
        commands = NULL
        idir.i = paste(idir, "/", src.dirs[i], sep="")
        odir.i = paste(odir, "/", src.dirs[i], sep="")
        files.i = list.files(idir.i)
        if (length(files.i) == 0)
            stop("no files found")
        system(paste("mkdir -p", odir.i))
        for (j  in 1:length(files.i)) {
            ifn = paste(idir.i, "/", files.i[j], sep="")
            ofn1 = paste(odir.i, "/", files.i[j], "_S", 1, sep="")
            ofn2 = paste(odir.i, "/", files.i[j], "_S", 2, sep="")

            command = sprintf("perl %s %s %s %s %s", script, ifn, self, ofn1, ofn2)
          # cat(sprintf("command: %s\n", command))
            commands = c(commands, command)
        }
        distrib(commands, working.dir=getwd(), qsub.dir=qsub.dir, jobname=jobname, batch.max.jobs=batch.max.jobs,
                total.max.jobs.fn=total.max.jobs.fn, sleep.time=10, rprofile=NULL,
                path=NULL, host.keys=NULL, retries=3, req.mem=NA, dtype=dtype)
    }
}

compute.matrix.single=function(mat.dir, ifn, ifn.genes, field, ofn, identity.threshold, coverage.threshold, min.gene.count, include.shared)
{
    ifn1 = ifn
    ifn2 = ifn
    ifn.genes1 = ifn.genes
    ifn.genes2 = ifn.genes
    field1 = field
    field2 = field
    compute.matrix(mat.dir=mat.dir,
                   ifn1=ifn1, ifn.genes1=ifn.genes1, field1=field1,
                   ifn2=ifn1, ifn.genes2=ifn.genes2, field2=field2,
                   ofn=ofn,
                   identity.threshold=identity.threshold, coverage.threshold=coverage.threshold, min.gene.count=min.gene.count, include.shared=include.shared)
}

compute.breakdown=function(genes, shared, map, identity.threshold, coverage.threshold, side, include.shared, fix.names=T)
{
    total = length(genes)
    map = map[!is.element(map$gene, shared),]
    if (dim(map)[1] != length(unique(map$gene)))
        stop("internal")

    mapped = length(unique(map$gene))
    unmapped = total - mapped - length(shared)

    exact  = map$identity >= 99 & map$coverage >= 0.99
    strong = map$identity >= 70 & map$coverage >= 0.7 & !exact
    weak = !(exact | strong)

    # compute mean on selected genes
    mean.identity.idx = map$identity > identity.threshold & map$coverage > coverage.threshold/100
    mean.identity.count = sum(mean.identity.idx) + length(shared)

    # include shared genes in mean identity
    if (include.shared) {
        n.mean.identity = if (mean.identity.count > 0) mean(c(rep(100,length(shared)), map$identity[mean.identity.idx])) else 0
    } else {
        n.mean.identity = if (mean.identity.count > 0) mean(map$identity[mean.identity.idx]) else 0
    }

    df = data.frame(
        total=total,
        unmapped=unmapped,
        mapped=mapped,
        exact=sum(exact),
        strong=sum(strong),
        weak=sum(weak),
        mean.identity=round(n.mean.identity, 2),
        mean.identity.count=mean.identity.count)

    if (fix.names)
        names(df) = paste(names(df), side, sep="")

    df
}

compute.matrix=function(mat.dir, ifn1, ifn.genes1, field1, ifn2, ifn.genes2, field2, ofn, identity.threshold, coverage.threshold, min.gene.count, include.shared)
{
    table1 = load.table(ifn1)
    table1$set = table1[,field1]
    table1 = unique(table1[,c("set", "gene")])
    s1 = split(table1$gene, table1$set)

    # all genes
    cat(sprintf("reading all genes table1: %s\n", ifn.genes1))
    all.table1 = read.delim(ifn.genes1)
    unassigned.genes1 = setdiff(unique(all.table1$gene), unique(table1$gene))

    table2 = load.table(ifn2)
    table2$set = table2[,field2]
    table2 = unique(table2[,c("set", "gene")])
    s2 = split(table2$gene, table2$set)

    # all genes
    cat(sprintf("reading all genes table2: %s\n", ifn.genes2))
    all.table2 = read.delim(ifn.genes2)
    unassigned.genes2 = setdiff(unique(all.table2$gene), unique(table2$gene))

    N1 = length(s1)
    N2 = length(s2)

    # treat self differently
    self = (ifn1 == ifn2)

    cat(sprintf("going over %d set pairs ", N1*N2))
    count = 0
    result = NULL
    for (i1 in 0:N1) {
        for (i2 in 0:N2) {
            if (i1 == 0 && i2 == 0)
                next
            if (count %% 10 == 0)
                cat(".")
            count = count + 1
            if (self && i1 == i2)
                next

            set1 = if(i1 > 0) names(s1)[i1] else "NONE"
            set2 = if(i2 > 0) names(s2)[i2] else "NONE"

            genes1 = if(i1 > 0) s1[[i1]] else unassigned.genes1
            genes2 = if(i2 > 0) s2[[i2]] else unassigned.genes2

            if (self) {
                shared = unique(intersect(genes1, genes2))
            } else {
                shared = NULL
            }
            n.shared = length(shared)

            ifn1 = paste(mat.dir, "/", set1, "/", set2, "_S1", sep="")
            ifn2 = paste(mat.dir, "/", set1, "/", set2, "_S2", sep="")
            if (file.exists(ifn1) && file.exists(ifn2)) {
                map1 = read.delim(ifn1)
                map2 = read.delim(ifn2)
                df1 = compute.breakdown(genes=genes1, shared=shared, map=map1,
                    identity.threshold=identity.threshold, coverage.threshold=coverage.threshold, side=1, include.shared=include.shared)
                df2 = compute.breakdown(genes=genes2, shared=shared, map=map2,
                    identity.threshold=identity.threshold, coverage.threshold=coverage.threshold, side=2, include.shared=include.shared)
            } else {
                df1 = data.frame(
                    total1=length(genes1), unmapped1=length(genes1)-n.shared,
                    mapped1=0, exact1=0, strong1=0, weak1=0, mean.identity1=0, mean.identity.count1=0)
                df2 = data.frame(
                    total2=length(genes2), unmapped2=length(genes2)-n.shared,
                    mapped2=0, exact2=0, strong2=0, weak2=0, mean.identity2=0, mean.identity.count2=0)
            }

            df = data.frame(set1=set1, df1, set2=set2, df2, shared=n.shared)

            if (df$mean.identity.count1>min.gene.count & df$mean.identity.count2>min.gene.count)
                df$identity = mean(c(df$mean.identity1, df$mean.identity2))
            else
                df$identity = 0

            result = rbind(result, df)
        }
    }
    cat(" done\n")

    cat(sprintf("saving result table: %s\n", ofn))
    write.table(result, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}

compute.order=function(sm, sets1, sets2, field, side=1)
{
    sets.main = if (side==1) sets1 else sets2
    sets.other = if (side==1) sets2 else sets1
    sm$main = if (side==1) sm$set1 else sm$set2
    sm$other = if (side==1) sm$set2 else sm$set1

    # factor
    alpha = 2

    result = NULL
    N = length(sets.other)
    for (set in sets.main) {
        tsm = sm[sm$main == set,]
        weights = (tsm[match(sets.other, tsm$other),field])^alpha
        weights = weights/sum(weights)
        coord = sum(1:N * weights)
        result = rbind(result, data.frame(set=set, coord=coord))
    }
    result[order(result$coord),"set"]
}

plot.sets.order=function(smat, sets1, sets2, str1, str2, self, side, field, no.bg=F)
{
    if (self) {
        sets1.plot = sets1
        sets2.plot = sets2
    } else {
        sets1.passive = compute.order(sm=smat, sets1=sets1, sets2=sets2, field=field, side=1)
        sets2.passive = compute.order(sm=smat, sets1=sets1, sets2=sets2, field=field, side=2)
        if (side == 1) {
            sets1.plot = sets1.passive
            sets2.plot = sets2
            } else {
                sets1.plot = sets1
                sets2.plot = sets2.passive
            }
    }

    sets1.str = str1[match(sets1.plot, sets1)]
    sets2.str = str2[match(sets2.plot, sets2)]

    if (no.bg)
        return (list(sets1=c(sets1.plot), label1=c(sets1.str),
                     sets2=c(sets2.plot), label2=c(sets2.str)))
    else
        return(list(sets1=c(sets1.plot, "NONE"), label1=c(sets1.str, "BG"),
                    sets2=c(sets2.plot, "NONE"), label2=c(sets2.str, "BG")))
}

plot.matrix.single=function(mat.ifn, min.identity, order.ifn, title, fdir)
{
    order.ifn1 = order.ifn
    order.ifn2 = order.ifn
    title1 = title
    title2 = title
    self = T

    plot.matrix(mat.ifn=mat.ifn, self=self, min.identity=min.identity,
                order.ifn1=order.ifn1, order.ifn2=order.ifn2,
                title1=title1, title2=title2, fdir=fdir)
}

plot.matrix=function(mat.ifn, min.identity, order.ifn1, order.ifn2, title1, title2, fdir, self=F, no.bg=T)
{
    library(gplots)

    smat = load.table(mat.ifn)
    smat$exact.ratio1 = 100 * smat$exact1 / smat$total1
    smat$exact.ratio2 = 100 * smat$exact2 / smat$total2
    if (self) {
        smat$shared.ratio1 = 100 * (smat$shared) / smat$total1
        smat$shared.ratio = 100 * (smat$shared) / pmin(smat$total1, smat$total2)
        smat$identity[smat$set1 == smat$set2] = 100
        smat$identity[smat$identity < min.identity] = min.identity
    }

    set.ids1 = read.delim(order.ifn1)$id
    set.ids2 = read.delim(order.ifn2)$id
    sets1 = read.delim(order.ifn1)$set
    sets2 = read.delim(order.ifn2)$set

    N1 = length(sets1) + ifelse(no.bg,0,1)
    N2 = length(sets2) + ifelse(no.bg,0,1)
    dim = c(N1, N2)

    cat(sprintf("plotting results in: %s\n", fdir))
    system(paste("mkdir -p", fdir))

    ncols = 256
    xlim = c(-0.1,N1+.1)
    ylim = c(-0.1,N2+.1)

    xlab = title1
    ylab = title2

    plot.mat=function(smat, breaks, colors, panel, field, default.value=0) {
        if (length(breaks) != length(colors))
            stop("#breaks != #colors")

        sorder = plot.sets.order(smat,
            sets1=sets1, sets2=sets2,
            str1=set.ids1, str2=set.ids2,
            field="exact.ratio2", self=self, side=2, no.bg=no.bg)

        smat$set1.i = match(smat$set1, sorder$sets1)
        smat$set2.i = match(smat$set2, sorder$sets2)
        smat = smat[!is.na(smat$set1.i) & !is.na(smat$set2.i),]

        m = smatrix2matrix(smat, dim=dim, i.field="set1.i", j.field="set2.i", value.field=field, default.value=default.value)
        if (self)
            m = (m + t(m)) / 2
        sm = matrix2smatrix(m)
        if (self)
            sm = sm[sm$i != sm$j,]

        # colors
        sm$col = panel[vals.to.cols(sm$value, breaks, ncols)]

        if (self) {
            ofn = paste(fdir, "/self_", field, ".png", sep="")
        } else {
            by.str = title1
            ofn = paste(fdir, "/", by.str, "_matrix_", field, ".png", sep="")
        }

        cat(sprintf("plotting figure: %s\n", ofn))
        png(ofn, 12*N2 + 100, 12*N1 + 100)
        par(mai=c(0.5,0.5,0.5,0.5))
        plot.new()
        plot.window(xlim=ylim, ylim=xlim)
        title(main=field, xlab=ylab, ylab=xlab)
        rect(sm$j-1, sm$i-1, sm$j, sm$i, col=sm$col, border=NA)
        mtext(text=sorder$label2, side=1, at=1:N2-0.5, las=2, line=0)
        mtext(text=sorder$label1, side=2, at=1:N1-0.5, las=2, line=0)
        # rect(0, 0, N1, N2, col=NA, border="darkgrey")
        dev.off()
    }

    colors.exact = c("white", "blue", "red", "orange")
    breaks.exact = c(0,10,20,100)
    panel.exact = make.color.panel(colors.exact)
    wlegend2(fdir=fdir, panel=panel.exact, breaks=breaks.exact, title="exact", width=300, height=150)

    breaks.shared = c(0,1,10,100)
    colors.shared = c("white","blue","red","orange")
    panel.shared = make.color.panel(colors.shared)
    wlegend2(fdir=fdir, panel=panel.shared, breaks=breaks.shared, title="shared", width=300, height=150)

    breaks.identity = c(60, 70, 80)
    colors.identity = c("blue", "red", "orange")
    breaks.identity = c(30, 50, 60, 80)
    colors.identity = c("darkblue", "blue", "red", "orange")
    panel.identity = make.color.panel(colors.identity)
    wlegend2(fdir=fdir, panel=panel.identity, breaks=breaks.identity, title="identity", width=300, height=150)

    plot.mat(smat=smat, breaks=breaks.identity, colors=colors.identity, panel=panel.identity, field="identity")
    if (self) {
        plot.mat(smat=smat, breaks=breaks.exact, colors=colors.exact, panel=panel.exact, field="shared.ratio1")
        plot.mat(smat=smat, breaks=breaks.shared, colors=colors.shared, panel=panel.shared, field="shared")

        xx = smat[smat$identity < 100 & smat$shared > 0,]
        if (dim(xx)[1] > 0) {
            fig.start(paste(fdir, "/shared_vs_distance.png", sep=""), width=600, height=600)
            x = 100 - xx$identity
            y = xx$shared
            plot.init(xlim=c(0, 45), ylim=c(0,100), xlab="% dissimilarity", ylab="# shared genes")
#            plot.init(xlim=range(x), ylim=range(y), xlab="% dissimilarity", ylab="# shared genes")
            points(x, y, pch="+")
            fig.end()
        }
    } else {
        plot.mat(smat=smat, breaks=breaks.exact, colors=colors.exact, panel=panel.exact, field="exact.ratio1")
    }

}

plot.matrix.full.single=function(mat.ifn, order.ifn, title, fdir)
{
    plot.matrix.full(mat.ifn=mat.ifn, self=T,
                     order.ifn1=order.ifn, order.ifn2=order.ifn,
                     title1=title, title2=title, fdir=fdir)
}

plot.matrix.full=function(mat.ifn, order.ifn1, order.ifn2, title1, title2, fdir, self=F)
{
    options(stringsAsFactors=F)
    library(gplots)

    smat = read.delim(mat.ifn)
    smat$exact.ratio1 = 100 * smat$exact1 / smat$total1
    smat$exact.ratio2 = 100 * smat$exact2 / smat$total2

    sets1 = read.delim(order.ifn1)$set
    sets2 = read.delim(order.ifn2)$set
    set.ids1 = read.delim(order.ifn1)$id
    set.ids2 = read.delim(order.ifn2)$id

    N1 = length(sets1)+1
    N2 = length(sets2)+1

    cat(sprintf("plotting results in: %s\n", fdir))
    system(paste("mkdir -p", fdir))

    ncols = 256
    xlim = c(0,N1)
    ylim = c(0,N2)

    xlab = title1
    ylab = title2

    sorder = plot.sets.order(smat,
        sets1=sets1, sets2=sets2,
        str1=set.ids1, str2=set.ids2,
        field="exact.ratio2", self=self, side=2, no.bg=F)

    smat$set1.i = match(smat$set1, sorder$sets1)
    smat$set2.i = match(smat$set2, sorder$sets2)

    get.coords=function(x, side) {
        gf = function(field) { x[,paste(field, side, sep="")] }
        shared = x$shared / gf("total")
        exact = shared + gf("exact") / gf("total")
        strong = exact + gf("strong") / gf("total")
        weak = strong + gf("weak") / gf("total")
        unmapped = weak + gf("unmapped") / gf("total")
        data.frame(shared=shared, exact=exact, strong=strong, weak=weak, unmapped=unmapped)
    }
    coords1 = get.coords(smat, 1)
    coords2 = get.coords(smat, 2)

    if (self) {
        ofn = paste(fdir, "/self_full.png", sep="")
    } else {
        by.str = title1
        ofn = paste(fdir, "/", by.str, "_matrix_full.png", sep="")
    }

    cat(sprintf("plotting figure: %s\n", ofn))
    png(ofn, 50*N2 + 100, 50*N1 + 100)
    plot.new()
    par(mai=c(0.5,0.5,0.5,0.5))
    plot.window(xlim=ylim, ylim=xlim)
    title(xlab=ylab, ylab=xlab)

    pf = function(field, col, border=NA) {
        rect(smat$set2.i-1, smat$set1.i-1, smat$set2.i-1 + coords2[,field], smat$set1.i-1 + coords1[,field], col=col, border=border)
    }
    pf(field="weak", col="lightgray")
    pf(field="strong", col="blue")
    pf(field="exact", col="red")
    pf(field="shared", col="orange")
    pf(field="unmapped", col=NA, border="black")

    mtext(text=sorder$label2, side=1, at=1:N2-0.5, las=2, line=0)
    mtext(text=sorder$label1, side=2, at=1:N1-0.5, las=2, line=0)

    dev.off()
}

compute.best.match=function(ifn, ofn1, ofn2)
{
    options(stringsAsFactors=F)
    table = load.table(ifn)
    table$ratio1 = table$exact1 / table$total1
    table$ratio2 = table$exact2 / table$total2

    sets1 = sort(unique(table$set1))
    sets2 = sort(unique(table$set2))

    get.best=function(main.i, other.i) {
        table$main.set = table[,paste("set", main.i, sep="")]
        table$other.set = table[,paste("set", other.i, sep="")]
        table$ratio = table[,paste("ratio", main.i, sep="")]

        main.sets = setdiff(sort(unique(table$main.set)), "NONE")
        result = NULL
        for (set in main.sets) {
            x = table[table$main.set == set, c("other.set", "ratio")]
            bg = x[x$other.set == "NONE","ratio"]
            if (length(bg) == 0) {
                bg = 0
            }
            x = x[x$other.set != "NONE",]
            ii = which.max(x$ratio)
            result = rbind(result, data.frame(set=set, best.match=x$other.set[ii], exact.ratio=x$ratio[ii], bg=bg))
        }
        result
    }

    r1 = get.best(main.i=1, other.i=2)
    r2 = get.best(main.i=2, other.i=1)

    save.table(r1, ofn=ofn1)
    save.table(r2, ofn=ofn2)
}

plot.match=function(ifn1, ifn2, fdir)
{
    library(igraph)
    t1 = load.table(ifn1)
    t2 = load.table(ifn2)
    t1$type = T
    t2$type = F
    vertices = c(t1$set, t2$set)
    id = c(rep("R", dim(t1)[1]), rep("P", dim(t2)[1]))
    edges = t1[,1:2]
    g = graph_from_data_frame(edges, directed=T, vertices=vertices)
    cc = as.data.frame(layout_as_bipartite(g, type=c(t1$type, t2$type)))

    vs = data.frame(name=vertices, id=id, x=cc[,1], y=cc[,2])
    vs$col = ifelse(vs$id == "R", "blue", "red")
    vs$pos = ifelse(vs$id == "R", 1, 1)

    es = data.frame(
        x1=vs$x[match(edges[,1], vs$name)],
        x2=vs$x[match(edges[,2], vs$name)],
        y1=vs$y[match(edges[,1], vs$name)],
        y2=vs$y[match(edges[,2], vs$name)])
    es$col = "lightgrey"

    N1 = dim(t1)[1]
    N2 = dim(t2)[1]
    N = max(N1, N2)
    xlim = c(-2, N + 2)
    ylim = c(-0.2, 1.2)

    fig.start(ofn=paste(fdir, "/match.png", sep=""), fdir=fdir, height=200 + N*30, width=400)

    plot.new()
    plot.window(xlim=ylim, ylim=xlim)
    points(vs$y, vs$x, col=vs$col, pch=19)
    segments(es$y1, es$x1, es$y2, es$x2, col=es$col)
    text(vs$y, vs$x, vs$name, pos=vs$pos)

    fig.end()
}

project.genes=function(
    mat.dir, tmp.dir, collapse.script,
    ifn1, field1, all.ifn1, match.ifn1, ifn.genes1,
    ifn2, field2, all.ifn2, match.ifn2, ifn.genes2,
    ofn1, ofn2)
{
    all1 = load.table(all.ifn1)
    all2 = load.table(all.ifn2)

    table1 = load.table(ifn1)
    table1$set = table1[,field1]
    table1 = unique(table1[,c("set", "gene")])
    s1 = split(table1$gene, table1$set)

    # all genes
    all.table1 = load.table(ifn.genes1)
    unassigned.genes1 = setdiff(unique(all.table1$gene), unique(table1$gene))

    table2 = load.table(ifn2)
    table2$set = table2[,field2]
    table2 = unique(table2[,c("set", "gene")])
    s2 = split(table2$gene, table2$set)

    # all genes
    all.table2 = load.table(ifn.genes2)
    unassigned.genes2 = setdiff(unique(all.table2$gene), unique(table2$gene))

    mtable1 = load.table(match.ifn1)
    mtable2 = load.table(match.ifn2)

    all1$best_set = ifelse(is.element(all1$best_hit, table2$gene), table2$set[match(all1$best_hit, table2$gene)], "none")
    all2$best_set = ifelse(is.element(all2$best_hit, table1$gene), table1$set[match(all2$best_hit, table1$gene)], "none")

    ll = list(side1=list(s=s1, ofn=ofn1), side2=list(s=s2, ofn=ofn2))
    for (side in 1:2) {
        s = ll[[side]]$s
        ofn = ll[[side]]$ofn
        N = length(s)
        result = NULL
        for (i in 1:N) {
            set = names(s)[i]
            if (side == 1) {
                set1 = set
                set2 = mtable1$best.match[match(set1, mtable1$set)]
                all = all1
                tfn = paste(mat.dir, "/", set1, "/", set2, "_S1", sep="")
                match = set2
            } else {
                set2 = set
                sets = mtable1$set[mtable1$best.match == set2]
                all = all2
                match = paste(sets, sep="", collapse=",")
                ifns = paste(mat.dir, "/", sets, "/", set2, "_S2", sep="")
                tfn = paste(tmp.dir, "/", paste(sets, sep="", collapse="_"), "_to_", set2, sep="")
                command = sprintf("perl %s coverage max %s %s", collapse.script, tfn, paste(ifns, sep="", collapse=" "))
                if (system(command) != 0)
                    stop(sprintf("command failed: %s\n", command))
            }
            map = read.delim(tfn)
            genes = s[[i]]

            if (!setequal(genes, all$gene[all$set == set]))
                stop("assert")

            in.map = !is.na(match(genes, map$gene))

            df = data.frame(set=set, gene=genes,
                all.target=all$best_hit[match(genes, all$gene)],
                all.identity=all$best_identity[match(genes, all$gene)],
                all.set=all$best_set[match(genes, all$gene)],
                match.target=ifelse(in.map, map$target[match(genes, map$gene)], "none"),
                match.identity=ifelse(in.map, map$identity[match(genes, map$gene)], 0))

            result = rbind(result, df)
        }
        save.table(result, ofn)
    }
}

plot.project=function(ifn.classify.ref.genes, ifn.classify.genes, ifn.classify.contigs, ifn.fwd, ifn.bck, order.ifn.fwd, order.ifn.bck, fdir)
{
    ref.gene.class = load.table(ifn.classify.ref.genes)

    gene.class = load.table(ifn.classify.genes)
    contig.class = load.table(ifn.classify.contigs)
    gene.class$ctype =contig.class$type[match(gene.class$contig, contig.class$contig)]

    table.fwd = load.table(ifn.fwd)
    table.bck = load.table(ifn.bck)

    table.fwd$gtype = ref.gene.class$type[match(table.fwd$gene, ref.gene.class$gene)]

    table.bck$gtype = gene.class$type[match(table.bck$gene, gene.class$gene)]
    table.bck$ctype = gene.class$ctype[match(table.bck$gene, gene.class$gene)]

    sets.fwd = load.table(order.ifn.fwd)
    sets.bck = load.table(order.ifn.bck)

    classify.fwd=function(table) {
        ifelse(table$match.identity >= 90, "match",
               ifelse(table$gtype != "found", paste("g", table$gtype, sep="_"),
                      ifelse(table$all.identity >= 95, ifelse(table$all.set == "none", "unclassified", "misclassified"), "other")))
    }
    classify.bck=function(table) {
        ifelse(table$match.identity >= 90, "match",
               ifelse(table$gtype != "match", paste("g", table$gtype, sep="_"),
                      ifelse(table$ctype != "normal",  paste("c", table$ctype, sep="_"),
                             ifelse(table$all.identity >= 95, "misclassified" ,"other"))))
    }

    table.fwd$type = classify.fwd(table.fwd)
    table.bck$type = classify.bck(table.bck)

    N.fwd = dim(sets.fwd)[1]
    fig.dir(fdir)

    plot.t=function(table, sets, title) {
        types = sort(unique(table$type))
        types = setdiff(types, "match")
        cols = rainbow(length(types))

        N = dim(sets)[1]

        sizes = table(factor(table$set, levels=sets$set))

        fig.start(ofn=paste(fdir, "/project_", title, "_genes_match_percent.png", sep=""), width=200 + N*25, height=600)
        m = table(table$type == "match", factor(table$set, levels=sets$set))
        m = 100 * (m[2,] / sizes)
        barplot(m, ylim=c(0,100), col="grey", border=NA, ylab="%", las=2, main=paste(title, "match %"), names.arg=sets$id)
        fig.end()

        table = table[table$type != "match",]
        fig.start(ofn=paste(fdir, "/project_", title, "_genes.png", sep=""), width=200 + N*25, height=600)
        m = table(factor(table$type, levels=types), factor(table$set, levels=sets$set))
        barplot(m, beside=F, col=cols, border=NA, las=2, names.arg=sets$id)
        fig.end()

        m = 100 * t((t(m)/as.vector(sizes)))
        fig.start(ofn=paste(fdir, "/project_", title, "_perc.png", sep=""), width=200 + N*25, height=600)
        barplot(m, beside=F, col=cols, border=NA, las=2, ylab="%", names.arg=sets$id)
        fig.end()

        wlegend(fdir=fdir, names=types, cols=cols, title=title)
    }

    plot.t(table=table.fwd, sets=sets.fwd, title="fwd")
    plot.t(table=table.bck, sets=sets.bck, title="bck")
}

plot.gene.project=function(fwd.source.ifn, bck.source.ifn, order.ifn1, order.ifn2, fdir)
{
    fwd.source = load.table(fwd.source.ifn)
    bck.source = load.table(bck.source.ifn)

    fwd.sets = load.table(order.ifn1)[,1]
    bck.sets = load.table(order.ifn2)[,1]

    types = c("exact", "weak", "missing")
    tcols = c("lightgreen", "grey", "red")

    fig.dir(fdir)

    ll = list(
        fwd=list(table=fwd.source, sets.order=fwd.sets),
        bck=list(table=bck.source, sets.order=bck.sets))

    result = NULL
    for (i in 1:2) {
        table = ll[[i]]$table
        title = names(ll)[i]
        sets.order = ll[[i]]$sets.order

        sets = sort(unique(table$set))
        N = length(sets)

        tt = table(table$set, table$type)
        df = data.frame(set=rownames(tt), exact=tt[,"found"], weak=tt[,"weak"], missing=tt[,"missing"])
        df = df[match(sets.order, df$set),]
        m = t(as.matrix(df[,types])) + 1

        fig.start(ofn=paste(fdir, "/breakdown_", title, ".png", sep=""), width=300 + N*30, height=600)
        par(mai=c(3,1,1,1))
        mp = barplot(m, beside=T, las=2, col=tcols, log="y", plot=F)
        barplot(m, beside=T, las=2, col=tcols, log="y", xlim=c(0, max(mp) * 1.2), border=NA, main=title)
        legend("topright", fill=tcols, legend=types, bg="white")
        fig.end()

        x = t(colSums(df[,-1]))
        x = x/sum(x)
        result = rbind(result, data.frame(type=title, x))
    }

    fig.start(ofn=paste(fdir, "/summary.png", sep=""), width=400, height=600)
    par(mai=c(3,1,1,1))
    barplot(100*(as.matrix(result[,-1])), beside=T, col=c("blue", "orange"), log="y", ylab="%", ylim=c(10^-2, 100))
    legend("topright", fill=c("blue", "orange"), legend=result$type, bg="white")
    fig.end()

}
