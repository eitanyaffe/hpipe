compute.rarify=function(folds, base.dir, fold0.table, ofn)
{
    result = NULL
    for (fold in sort(folds)) {
        dir = paste(base.dir, "/fold_", fold, sep="")
        cat(sprintf("processing: %s\n", dir))

        # count ontigs
        x = load.table(paste(dir, "/contig_table", sep=""))
        long = x$length > 1000

        # count reads
        if (fold == 0) {
            command = sprintf("cat `cat %s` | wc -l", fold0.table)
        } else
            command = paste("cat ", dir, "/rarify/reads.fasta | wc -l", sep="")
        cat(sprintf("command: %s\n", command))
        N = as.numeric(system(command, intern=T))

        df = data.frame(fold=fold, reads=N,
            contigs=dim(x)[1], median=median(x$length), total=sum(x$length),
            contigs.long=sum(long), median.long=median(x$length[long]), total.long=sum(x$length[long]))
        result = rbind(result, df)
    }
    save.table(result, ofn=ofn)
}

plot.rarify=function(ifn, fdir)
{
    table = load.table(ifn)
    fig.dir(fdir)

    fields = c("contigs", "median", "total", "contigs.long", "median.long", "total.long")
    names = c("#contigs", "median(contig)", "total bp", "#contigs >1k", "median(contig >1k)", "total bp >1k")
    for (i in 1:length(fields)) {
        field = fields[i]
        ylim = range(c(0, table[,field]))
        fig.start(ofn=paste(fdir, "/", field, ".png", sep=""))
        par(mai=c(1, 0.5, 1.2, 0.5))
        plot(table$reads/10^6, table[,field], type="b", log="x", xlab="M reads", ylab=names[i], main=names[i], ylim=ylim)
        grid()
        axis(side=3, at=table$reads/10^6, labels=table$fold, cex.axis=0.8)
        fig.end()
    }
}
