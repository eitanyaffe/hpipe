compare=function(
    ifn1="/relman02/users/eitany/bcc_output/Relman/assembly/pre_minia_contig/fold_0/datasets/pre_hic/anchors/pre_hic/ca_matrix/contigs",
    ifn2="/relman02/users/eitany/bcc_output/Relman/assembly/pre_minia_contig/fold_0/datasets/post_lib_hic_simple/anchors/pre_hic/ca_matrix/contigs",
)
{
    title1="pre"
    title2="post"
    fdir="figs"
    ifn1 = "/relman02/users/eitany/bcc_output/Relman/assembly/pre_minia_contig/fold_0/datasets/pre_hic/anchors/pre_hic/ca_matrix/contigs"
    ifn2 = "/relman02/users/eitany/bcc_output/Relman/assembly/pre_minia_contig/fold_0/datasets/post_lib_hic_simple/anchors/pre_hic/ca_matrix/contigs"
    t1 = read.delim(ifn1)
    t2 = read.delim(ifn2)


    m = merge(t1, t2, by=c("anchor", "contig"), all=T)
     m = m[m$contig_total_count.x > 30 | m$contig_total_count.y > 30,]
#    m[,title1] = m$contig_total_count.x
#    m[,title2] = m$contig_total_count.y

    m[,title1] = m$enrichment.x
    m[,title2] = m$enrichment.y

    m[is.na(m[,title1]),title1] = 1
    m[is.na(m[,title2]),title2] = 1

    s = split(m, m$anchor)
    anchors = names(s)

    lim = range(c(m[,title1], m[,title2]))
    system(paste("mkdir -p", fdir))
    for (anchor in anchors) {
        ofn = paste(fdir, "/", anchor, ".png", sep="")
        x = s[[anchor]]
        png(ofn, 400, 400)
        plot(x[,title1], x[,title2], type="p", pch="+", main=anchor, xlab=title1, ylab=title2, xlim=lim, ylim=lim)
        abline(a=0,b=1, col="grey")
        dev.off()
    }

}
