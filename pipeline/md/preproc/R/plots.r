plot.complex=function(ldir, fdir, ids, titles, cols)
{
    N = length(ids)
#    cols = c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")[1:N]
    cols = cols[1:N]

    tt = NULL
    for (id in ids) {
        ifn = paste(ldir, "/", id, "/complexity.table", sep="")
        x = load.table(ifn)
        ii = which(diff(x$multi) != 1)[1]
        if (!is.na(ii)) {
            if (ii > 80) ii = 80
            head = x[1:ii,]
            tail = x[(ii+1):dim(x)[1],]
            tail.m = head$multi[dim(head)[1]]
            tail.v = sum(tail$count)
            names = c(head$multi, paste(">", tail.m, sep=""))
            values = c(head$count, tail.v)
        } else {
            names = x$multi
            values = x$count
        }

        dups = sum(x$count[x$multi>1])
        all = sum(x$count)

        tt = rbind(tt, data.frame(id=id, dup=100 * dups/all))
        main = sprintf("%s, dup molecules: %d (%%%.2f)", id, dups, 100 * dups/all)
        fig.start(fdir=fdir, width=1000, height=400, ofn=paste(fdir, "/", id, ".png", sep=""))
        barplot(values+1, names.arg=names, main=main, log="y", las=2, border=NA, ylim=c(1, max(values+1)))
        fig.end()

    }

    fig.start(fdir=fdir, width=600, height=400, ofn=paste(fdir, "/summary.png", sep=""))
    par(mai=c(2,1.5,1,0.5))
    barplot(tt$dup, names.arg=tt$id, ylab="%", border=NA, main="percentage of duplicated molecules", col=cols, las=2)
    fig.end()

}

plot.stats=function(ldir, fdir, ids, titles, cols)
{
    N = length(ids)
#    cols = c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")[1:N]
    cols = cols[1:N]

    fig.start(fdir=fdir, width=600, height=200 + N*40, ofn=paste(fdir, "/id_legend.png", sep=""))
    par(mai=c(2,1.5,1,0.5))
    plot.new()
    legend("center", fill=cols, legend=ids)
    fig.end()

    par(mai=c(2,1.5,1,0.5))

    get.reads=function(x, side) { x[x[,2] == side,3] }
    get.bps=function(x, side) { x[x[,2] == side,4] }

    reads.count.m = NULL
    reads.loss.m = NULL
    bps.count.m = NULL
    bps.loss.m = NULL

    width = 350

    for (id in ids) {
        reads = NULL
        bps = NULL

        x = load.table(paste(ldir, "/", id, "/.count_input", sep=""))
        reads = rbind(reads, data.frame(id=id, type="input_P", count=x$input))
        reads = rbind(reads, data.frame(id=id, type="dups_P", count=x$output))

        x = load.table(paste(ldir, "/", id, "/.count_dup", sep=""), header=F)
        bps = rbind(bps, data.frame(id=id, type="dups_R1", count=get.bps(x,"R1")))
        bps = rbind(bps, data.frame(id=id, type="dups_R2", count=get.bps(x,"R2")))

        x = load.table(paste(ldir, "/", id, "/.count_trim", sep=""), header=F)
        bps = rbind(bps, data.frame(id=id, type="trim_R1", count=get.bps(x,"R1")))
        bps = rbind(bps, data.frame(id=id, type="trim_R2", count=get.bps(x,"R2")))
        reads = rbind(reads, data.frame(id=id, type="trim_P", count=get.reads(x,"R1")))

        x = load.table(paste(ldir, "/", id, "/.count_adaptor", sep=""), header=F)
        bps = rbind(bps, data.frame(id=id, type="adaptor_R1", count=get.bps(x,"R1")))
        bps = rbind(bps, data.frame(id=id, type="adaptor_R2", count=get.bps(x,"R2")))
        reads = rbind(reads, data.frame(id=id, type="adaptor_P", count=get.reads(x,"R1")))

        x = load.table(paste(ldir, "/", id, "/.count_human_simple", sep=""), header=F)
        bps = rbind(bps, data.frame(id=id, type="human_R1", count=get.bps(x,"R1")))
        bps = rbind(bps, data.frame(id=id, type="human_R2", count=get.bps(x,"R2")))
        reads = rbind(reads, data.frame(id=id, type="human_R1", count=get.reads(x,"R1")))
        reads = rbind(reads, data.frame(id=id, type="human_R2", count=get.reads(x,"R2")))

        gloss = function(n, total) { 100 * (total-n) / total }
        reads$loss = 0
        reads$loss[reads$type == "dups_P"] = gloss(reads$count[reads$type == "dups_P"], reads$count[reads$type == "input_P"])
        reads$loss[reads$type == "trim_P"] = gloss(reads$count[reads$type == "trim_P"], reads$count[reads$type == "dups_P"])
        reads$loss[reads$type == "adaptor_P"] = gloss(reads$count[reads$type == "adaptor_P"], reads$count[reads$type == "trim_P"])
        reads$loss[reads$type == "human_R1"] = gloss(reads$count[reads$type == "human_R1"], reads$count[reads$type == "adaptor_P"])
        reads$loss[reads$type == "human_R2"] = gloss(reads$count[reads$type == "human_R2"], reads$count[reads$type == "adaptor_P"])

        bps$loss = 0
        bps$loss[bps$type == "trim_R1"] = gloss(bps$count[bps$type == "trim_R1"], bps$count[bps$type == "dups_R1"])
        bps$loss[bps$type == "trim_R2"] = gloss(bps$count[bps$type == "trim_R2"], bps$count[bps$type == "dups_R2"])
        bps$loss[bps$type == "adaptor_R1"] = gloss(bps$count[bps$type == "adaptor_R1"], bps$count[bps$type == "trim_R1"])
        bps$loss[bps$type == "adaptor_R2"] = gloss(bps$count[bps$type == "adaptor_R2"], bps$count[bps$type == "trim_R2"])
        bps$loss[bps$type == "human_R1"] = gloss(bps$count[bps$type == "human_R1"], bps$count[bps$type == "adaptor_R1"])
        bps$loss[bps$type == "human_R2"] = gloss(bps$count[bps$type == "human_R2"], bps$count[bps$type == "adaptor_R2"])

        reads.count.m = rbind(reads.count.m, reads$count)
        colnames(reads.count.m) = reads$type

        reads.loss.m = rbind(reads.loss.m, reads$loss[-1])
        colnames(reads.loss.m) = reads$type[-1]

        bps.count.m = rbind(bps.count.m, bps$count)
        colnames(bps.count.m) = bps$type

        bps.loss.m = rbind(bps.loss.m, bps$loss)
        colnames(bps.loss.m) = bps$type

        ffdir = paste(fdir, "/datasets/", id, sep="")

        fig.start(fdir=ffdir, width=width, height=400, ofn=paste(ffdir, "/read_count.png", sep=""))
        barplot(reads$count, names.arg=reads$type, border=NA, main="read count (M)", las=2, col="darkblue", ylab="M reads")
        fig.end()
        fig.start(fdir=ffdir, width=width, height=400, ofn=paste(ffdir, "/read_loss.png", sep=""))
        barplot(reads$loss[-1], names.arg=reads$type[-1], border=NA, main="read loss", las=2, col="darkblue", ylab="%")
        fig.end()
        fig.start(fdir=ffdir, width=width, height=400, ofn=paste(ffdir, "/bp_count.png", sep=""))
        barplot(bps$count, names.arg=bps$type, border=NA, main="bp count (M)", las=2, col="darkblue", ylab="M bps")
        fig.end()
        fig.start(fdir=ffdir, width=width, height=400, ofn=paste(ffdir, "/bp_loss.png", sep=""))
        barplot(bps$loss, names.arg=bps$type, border=NA, main="bp loss", las=2, col="darkblue", ylab="%")
        fig.end()

    }

    # plot comparison
    fig.start(fdir=fdir, width=width, height=400, ofn=paste(fdir, "/read_count.png", sep=""))
    par(mai=c(2,1.5,1,0.5))
    barplot(reads.count.m/10^6, beside=T, border=NA, main="read count (M)", las=2, col=cols)
    fig.end()

    fig.start(fdir=fdir, width=width, height=400, ofn=paste(fdir, "/read_loss.png", sep=""))
    par(mai=c(2,1.5,1,0.5))
    barplot(reads.loss.m, beside=T, border=NA, main="read loss %", las=2, col=cols)
    fig.end()

    fig.start(fdir=fdir, width=width, height=400, ofn=paste(fdir, "/bp_count.png", sep=""))
    par(mai=c(2,1.5,1,0.5))
    barplot(bps.count.m/10^9, beside=T, border=NA, main="nt count (G)", las=2, col=cols)
    fig.end()

    fig.start(fdir=fdir, width=width, height=400, ofn=paste(fdir, "/bp_loss.png", sep=""))
    par(mai=c(2,1.5,1,0.5))
    barplot(bps.loss.m, beside=T, border=NA, main="bps loss %", las=2, col=cols)
    fig.end()
}
