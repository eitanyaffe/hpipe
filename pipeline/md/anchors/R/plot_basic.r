plot.total.size=function(ifn.contigs, ifn.ca, fdir)
{
    ca = load.table(ifn.ca)
    contigs = load.table(ifn.contigs)
    ca$length = contigs$length[match(ca$contig,contigs$contig)]
    s = sapply(split(ca$length, ca$anchor), sum)/1000
    s = sort(s,decreasing=T)
    fig.start(fdir=fdir, ofn=paste(fdir, "/union_size.pdf", sep=""), type="pdf", width=10, height=4)
    par(mai=c(1,1,0.5,0.2))
    barplot(s, names.arg=1:length(s), las=2, ylab="kb")
    grid()
    fig.end()
}

