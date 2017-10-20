plot.class=function(ifn.elements, ifn.class, fdir)
{
    clusters = load.table(ifn.elements)$cluster
    df = load.table(ifn.class)
    df = df[is.element(df$cluster, clusters),]
    df$type =
        ifelse(df$phage & !df$plasmid, "phage",
               ifelse(!df$phage & df$plasmid, "plasmid",
                             ifelse(df$phage | df$plasmid, "mixed",
                                    "none")))

    df$f = factor(df$type, levels=c("phage", "plasmid", "mixed", "none"))
    tt = table(df$f)
    an = sum(tt) - tt["none"]
    main = paste("annotated=", an, " \n%", round(100*an/sum(tt),1), sep="")
    fig.start(fdir=fdir, ofn=paste(fdir, "/element_class.pdf", sep=""), type="pdf", width=1+length(tt)*0.3, height=4)
    barplot(tt, las=2, col="darkgreen", border=NA, ylab="#elements", ylim=c(0,1.1*max(tt)), main=main, cex.main=0.75)
    fig.end()
}
