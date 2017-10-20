plot.stats=function(input.ifn, filter.ifn, pair.ifn, parse.ifn, id, fdir)
{
    input = as.matrix(load.table(input.ifn))*2
    parse = load.table(parse.ifn)
    filter = load.table(filter.ifn)
    pair = load.table(pair.ifn)

    united = as.matrix(data.frame(filter, parse[,-1]))
    ucols.def = unlist(list(
        ok="darkgreen",
        short_match_length="darkred",
        high_edit_dist="red",
        low_score="orange",
        unmapped="gray",
        multi_segment="darkgrey",
        duplicate="darkgrey"))
    ucols = ucols.def[match(colnames(united), names(ucols.def))]
    ylim = c(0, 1.1 * max(c(sum(input), sum(parse), sum(filter))))
    fig.start(fdir=fdir, ofn=paste(fdir, "/side_map_stats.png", sep=""))
    layout(matrix(1:2, 1, 2, ), widths=c(1,2))
    barplot(t(united), las=2, ylim=ylim, border=NA, col=ucols, xlim=c(0, 1.25))
    abline(h=input, lwd=2, lty=2)
#    abline(h=parse$ok, lty=2, lwd=2)
    plot.legend(fill=ucols, text=names(ucols), main=paste("side", id))
    fig.end()

    ylim = c(0, 1.1 * max(c(filter$ok/2, sum(pair))))
    pcols.def = unlist(list(
        ok="darkgreen",
        no_pair="orange"))
    pcols = pcols.def[match(colnames(pair), names(pcols.def))]
    fig.start(fdir=fdir, ofn=paste(fdir, "/pair_map_stats.png", sep=""))
    layout(matrix(1:2, 1, 2, ), widths=c(1,2))
    barplot(t(pair), las=2, ylim=ylim, border=NA, col=pcols, xlim=c(0, 1.25))
#    abline(h=filter$ok/2, lty=2, lwd=2)
    plot.legend(fill=pcols, text=names(pcols), main=paste("pair", id))
    fig.end()
}

plot.stats.summary=function(map.tag, adir, aid, ids, titles, filter.id, fdir)
{
    table = NULL
    table.pair = NULL
    for (i in 1:length(ids)) {
        id = ids[i]
        title = titles[i]

        input.ifn = paste(adir, "/datasets/", id, "/map_", map.tag, "/stats_input", sep="")
        parse.ifn = paste(adir, "/datasets/", id, "/map_", map.tag, "/parse_stat.table", sep="")
        filter.ifn = paste(adir, "/datasets/", id, "/map_", map.tag, "/filter_stat_", filter.id, ".table", sep="")
        pair.ifn = paste(adir, "/datasets/", id, "/map_", map.tag, "/pair_stat_", filter.id, ".table", sep="")

        input = as.matrix(load.table(input.ifn))*2
        parse = load.table(parse.ifn)
        filter = load.table(filter.ifn)
        pair = load.table(pair.ifn)

        united = as.matrix(data.frame(filter, parse[,-1]))
        table = rbind(table, data.frame(id=title, input=input[1], united))

        table.pair = rbind(table.pair, data.frame(id=title, as.matrix(pair)))
    }
    N = dim(table)[1]

    ucols.def = unlist(list(
        ok="darkgreen",
        short_match_length="darkred",
        high_edit_dist="red",
        low_score="orange",
        unmapped="grey",
        multi_segment="black",
        duplicate="black"))

    pcols.def = unlist(list(
        ok="darkgreen",
        no_pair="orange"))

    #####################
    # one side
    #####################

    ucols = ucols.def[match(colnames(table[,-(1:2)]), names(ucols.def))]

    wbarplot(m=t(as.matrix(table[,-(1:2)]))/10^6, names=titles, main=paste("# sides", aid),
             fdir=fdir, ofn=paste(fdir, "/", aid, "_side.png", sep=""),
             beside=F, normalize.columns=F, cols=ucols, ylab="M sides")
    wbarplot(m=t(as.matrix(table[,-(1:2)]))/10^6, names=titles, main=paste("% sides", aid),
             fdir=fdir, ofn=paste(fdir, "/", aid, "_side_percent.png", sep=""),
             beside=F, normalize.columns=T, cols=ucols, ylab="% sides")
    wlegend(fdir=fdir, names=names(ucols), cols=ucols, title="sides")

    #####################
    # pairs
    #####################

    pcols = pcols.def[match(colnames(table.pair[,-1]), names(pcols.def))]

    wbarplot(m=t(as.matrix(table.pair[,-1]))/10^6, names=titles, main=paste("# reads", aid),
             fdir=fdir, ofn=paste(fdir, "/", aid, "_read_pair.png", sep=""),
             beside=F, normalize.columns=F, cols=pcols, ylab="M reads")
    wbarplot(m=t(as.matrix(table.pair[,-1]))/10^6, names=titles, main=paste("% reads", aid),
             fdir=fdir, ofn=paste(fdir, "/", aid, "_read_pair_percent.png", sep=""),
             beside=F, normalize.columns=T, cols=pcols, ylab="% reads")
    wlegend(fdir=fdir, names=names(pcols), cols=pcols, title="reads")

    # add total and save to table
    table = rbind(table, c(data.frame(id="total"), colSums(table[,-1])))
    table.pair = rbind(table.pair, c(data.frame(id="total"), colSums(table.pair[,-1])))

    save.table(table, paste(fdir, "/read_sides.txt", sep=""))
    save.table(table.pair, paste(fdir, "/pairing.txt", sep=""))
}
