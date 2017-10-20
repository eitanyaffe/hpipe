plot.compare.cc.matrix=function(adir, aid, ids, titles, fdir)
{
    ll = list()
    keys = NULL
    for (i in 1:length(ids)) {
        id = ids[i]
        ifn = paste(adir, "/datasets/", id, "/contig_mat/table", sep="")
        x = load.table(ifn)
        x = x[x$contig1 != x$contig2 ,]
        key = paste(x$contig1, x$contig2, sep="_")
        ll[[id]] =  data.frame(key=key, value=x$contacts)
        keys = unique(c(keys, key))
    }

    table = data.frame(key=keys)
    for (i in 1:length(ids)) {
        id = ids[i]
        x = ll[[id]]
        table[,id] = 0
        ix = match(x$key, table$key)
        table[ix,id] = x$value
    }
    N = length(ids)
    fig.start(fdir=fdir, ofn=paste(fdir, "/all.png", sep=""), height=600, width=600)
    par(mai=c(0.5,0.5,0.5,0.5))
    layout(matrix(1:N^2,N,N))
    for (i1 in 1:N)
    for (i2 in 1:N) {
        id1 = ids[i1]
        id2 = ids[i2]
        if (i1 == i2) {
            plot.new()
            plot.window(xlim=0:1, ylim=0:1)
            text(0.5, 0.5, titles[i1])
        }
        if (i1 < i2) {
            ix = table[,id1] > 0 & table[,id2] > 0
            cc = round(cor(table[ix,id1], table[ix,id2], method="spearman"),2)
            plot(table[ix,id1], table[ix,id2], pch=".", main=cc, log="xy", xlab="", ylab="")
        }
        if (i1 > i2) {
            plot.new()
        }
    }
    fig.end()
}
