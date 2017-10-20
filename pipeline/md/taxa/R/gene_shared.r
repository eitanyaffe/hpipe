plot.shared=function(ifn, ga.ifn, fdir)
{
    table = read.delim(ifn)
    system(paste("mkdir -p", fdir))

    ga = read.delim(ga.ifn)
    anchors = sort(unique(ga$anchor))

    s = split(ga$gene, ga$anchor)
    N = length(s)
    for (i in 1:N)
        for (j in 1:N) {
            if (i == j)
                next
            anchor1 = names(s)[i]
            anchor2 = names(s)[j]
            genes = intersect(s[[as.character(anchor1)]], s[[as.character(anchor2)]])
            count = length(genes)
            if (count <= 2)
                next
            ttable = table[is.element(table$gene, genes),]
            main =  paste(anchor1, " ", anchor2, sep="")
            x = table(ttable$prot_desc)
            x = x[order(x, decreasing=T)]
            x = x[x>0]
            if (length(x) < 2)
                next
            if (length(x) > 40)
                x = x[1:40]
            ofn = paste(fdir, "/", anchor1, "_", anchor2, ".png", sep="")
            cat(sprintf("plotting: %s\n", ofn))
            png(ofn, units="in", width=12, height=1 + 0.3*length(x), res=96)
#            png(ofn, width=800, height=200 + 12*length(x))
            par(mai=c(0.5, 8, 0.1, 0.1))
            barplot(x, horiz=T, las=2, main=main, log="x", cex.names=1.2)
            dev.off()
        }
}
