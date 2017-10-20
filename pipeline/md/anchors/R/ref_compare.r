smatrix2matrix=function(smatrix, dim, i.field="i", j.field="j", value.field="value", default.value=0)
{
  indices = smatrix[,i.field] + (smatrix[,j.field]-1) * dim[1]
  v = rep(default.value, dim[1]*dim[2])
  v[indices] = smatrix[,value.field]
  matrix(v, dim[1], dim[2])
}

matrix2smatrix=function(matrix)
{
  dim1 = dim(matrix)[1]
  dim2 = dim(matrix)[2]
  v = as.vector(matrix)
  indices = 1:(dim1*dim2)
  i = (indices-1) %% dim1 + 1
  j = floor((indices-1) / dim1) + 1
  data.frame(i=i, j=j, value=v, stringsAsFactors=F)
}

make.color.panel=function(colors, ncols=256)
{
  panel = NULL
  for (i in 2:length(colors))
    panel = c(panel, colorpanel(ncols, colors[i-1], colors[i]))
  panel
}

# for example vals.to.cols(1:10, c(1, 3, 10), ncols=10) returns:
# [1] 1  6 11 12 14 15 16 17 19 20
vals.to.cols=function(vals, breaks, ncols=256)
{
  min = breaks[1]
  max = breaks[length(breaks)]
  vals = ifelse(vals < min, min, ifelse(vals>max, max, vals))
  n = length(breaks)-1
  cols = rep(-1, length(vals))
  for (i in 1:n)
  {
    ind = (breaks[i] <= vals) & (vals <= breaks[i+1])
    if (!any(ind))
      next
    # normalize to [0,1]
    cols[ind] = (vals[ind] - breaks[i]) / (breaks[i+1] - breaks[i])
    # normalize to [i*ncols,i*(ncols+1)]
    cols[ind] = (i-1)*ncols + cols[ind]*(ncols-1) + 1
    # round
    cols[ind] = round(cols[ind])
  }
  return (cols)
}

compute.matrix=function(mat.dir, ifn1, field1, ifn2, field2, ofn)
{
    options(stringsAsFactors=F)
    cat(sprintf("reading table1: %s\n", ifn1))
    table1 = read.delim(ifn1)
    table1$set = table1[,field1]
    table1 = unique(table1[,c("set", "cgene")])
    s1 = split(table1$cgene, table1$set)

    cat(sprintf("reading table2: %s\n", ifn2))
    table2 = read.delim(ifn2)
    table2$set = table2[,field2]
    table2 = unique(table2[,c("set", "cgene")])
    s2 = split(table2$cgene, table2$set)

    N1 = length(s1)
    N2 = length(s2)

    result = NULL
    for (i1 in 1:N1) {
        for (i2 in 1:N2) {
            set1 = names(s1)[i1]
            set2 = names(s2)[i2]

            ifn = paste(mat.dir, "/", set1, "_", set2, sep="")
            if (!file.exists(ifn))
                next
            map = read.delim(ifn)

            cgenes1 = s1[[i1]]
            cgenes2 = s2[[i2]]

            n1 = length(cgenes1)
            n2 = length(cgenes2)

            # unmapped cgenes
            n.unmapped1 = n1 - length(unique(map$cgene1))
            n.unmapped2 = n2 - length(unique(map$cgene2))

            exact = map$identity > 99 & map$coverage >= 0.99
            strong = map$identity >= 70 & map$coverage >= 0.7 & !exact
            weak = map$identity < 70 & map$coverage < 0.7
            n.exact = sum(exact)
            n.strong = sum(strong)
            n.weak = sum(weak)

            # compute average identity score
            identity = round(mean(map$identity),2)
            identity.count = length(map$identity)

            df = data.frame(set1=set1, set2=set2, total1=n1, total2=n2,
                identity=identity, identity.count=identity.count, exact=n.exact, strong=n.strong, weak=n.weak,
                unmapped1=n.unmapped1, unmapped2=n.unmapped2)
            result = rbind(result, df)
        }
    }

    cat(sprintf("saving result table: %s\n", ofn))
    write.table(result, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}

plot.matrix=function(mat.ifn, order.ifn1, order.ifn2, title1, title2, fdir)
{
    options(stringsAsFactors=F)
    library(gplots)

    smat = read.delim(mat.ifn)
    smat$exact.ratio = 100 * smat$exact / pmin(smat$total1, smat$total2)

    sets1 = read.delim(order.ifn1)[,1]
    sets2 = read.delim(order.ifn2)[,1]

    N1 = length(sets1)
    N2 = length(sets2)
    dim = c(N1,N2)

    cat(sprintf("plotting results in: %s\n", fdir))
    system(paste("mkdir -p", fdir))

    ncols = 256
    xlim = c(0,N1)
    ylim = c(0,N2)

    compute.order=function(sm, field, side=1) {
        sets.main = if (side==1) sets1 else sets2
        sets.other = if (side==1) sets2 else sets1
        sm$main = if (side==1) sm$set1 else sm$set2
        sm$other = if (side==1) sm$set2 else sm$set1
        result = NULL
        N = length(sets.other)
        for (set in sets.main) {
            tsm = sm[sm$main == set,]
            weights = tsm[match(sets.other, tsm$other),field]
            weights = weights/sum(weights)
            coord = sum(1:N * weights)
            result = rbind(result, data.frame(set=set, coord=coord))
        }
        result[order(result$coord),"set"]
    }

    sets1.passive = compute.order(sm=smat, field="exact.ratio", side=1)
    sets2.passive = compute.order(sm=smat, field="exact.ratio", side=2)

    plot.mat=function(breaks, colors, field, default.value=0, side) {
        if (length(breaks) != length(colors))
            stop("#breaks != #colors")

        sm = smat[,c("set1", "set2", field)]

        if (side == 1) {
            sets1.plot = sets1.passive
            sets2.plot = sets2
            by.str = title2
        } else {
            sets1.plot = sets1
            sets2.plot = sets2.passive
            by.str = title1
        }

        sets1.str = paste("R", match(sets1.plot, sort(sets1)), sep="")
        sets2.str = sets2.plot

        sm$set1.i = match(smat$set1, sets1.plot)
        sm$set2.i = match(smat$set2, sets2.plot)
        m = smatrix2matrix(sm, dim=dim, i.field="set1.i", j.field="set2.i", value.field=field, default.value=default.value)
        sm = matrix2smatrix(m)

        # colors
        panel = make.color.panel(colors, ncols)
        sm$col = panel[vals.to.cols(sm$value, breaks, ncols)]

        ofn = paste(fdir, "/matrix_", field, "_by_", by.str, ".png", sep="")
        cat(sprintf("plotting figure: %s\n", ofn))
        png(ofn, 12*N1 + 100, 12*N2 + 100)
        par(mai=c(0.5,0.5,0.5,0.5))
        plot.new()
        plot.window(xlim=xlim, ylim=ylim)
        title(main=field)
        rect(sm$i-1, sm$j-1, sm$i, sm$j, col=sm$col, border=NA)
        mtext(text=sets1.str, side=1, at=1:N1-0.5, las=2, line=0)
        mtext(text=sets2.str, side=2, at=1:N2-0.5, las=2, line=0)
        dev.off()
    }

    for (side in 1:2) {
        plot.mat(breaks=c(60,80,95,100), colors=c("white","blue","red","orange"), field="identity", side=side)
        plot.mat(breaks=c(0,50,90,100), colors=c("white", "blue", "red", "orange"), field="exact.ratio", side=side)
    }
}
