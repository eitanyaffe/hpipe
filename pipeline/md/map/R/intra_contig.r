plot.intra=function(ifn, idir.coverage, idir.contacts, threshold, odir)
{
    contigs = load.table(ifn)

#    x = read.delim("/relman02/users/eitany/bcc_output/Relman/assembly/pre/fold_0/datasets/pre_hic/anchors/pre_hic/ca_matrix/contigs")
#    x = field.count(x, "contig")
 #   contigs$count = x$count[match(contigs$contig, x$contig)]
 #   contigs = contigs[contigs$count > 1 & !is.na(contigs$count),]

#    contigs = contigs[contigs$length > threshold,]
    cat(sprintf("number of contigs: %d\n", dim(contigs)[1]))

    fig.dir(odir)
    for (contig in contigs$contig) {
        cov = load.table(paste(idir.coverage, "/", contig, sep=""))
        contacts = load.table(paste(idir.contacts, "/", contig, sep=""))

        rot = data.frame(x=(contacts$coord1 + contacts$coord2)/2, y=(contacts$coord2 - contacts$coord1)/2)

        fig.start(paste(odir, "/", contig, ".png", sep=""), width=800, height=400)

        par(mai=c(0,0.8,0,0))
        layout(matrix(1:2, 2, 1), height=c(10,3))
        lim = range(cov$coord)

        # contacts
        plot.new()
        plot.window(xlim=lim, ylim=c(0, max(lim)/2))
        points(rot$x, rot$y, pch=".")
        grid()
        box()

        # coverage
        par(mai=c(0.5,0.8,0,0))
        plot(cov$coord, cov$count, type="l", ylim=c(0,max(cov$count)), las=1, ylab="", xlab="")
        lines(cov$coord, smooth(cov$count), col=2)
        grid()

        fig.end()
    }
}
