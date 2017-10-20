plot.GO.internal=function(df, type, min.count, fdir)
{
    df = df[df$count >= min.count & df$type == type,]
    df$enrichment = log2(df$enrichment)
    N = dim(df)[1]

    pos = rep(4, N)
    # pos = round(runif(N,0,3))+1

    s1 = c("methyltransferase activity")
    s2 = c("integral component of membrane", "DNA-templated transcription, initiation", "nucleic acid phosphodiester bond hydrolysis", "C-5 methylation of cytosine", "cytoplasm")
    s3 = c("kinase activity", "methylation")
    s4 = c("DNA (cytosine-5-)-methyltransferase activity")
    pos[is.element(df$desc, s1)] = 1
    pos[is.element(df$desc, s2)] = 2
    pos[is.element(df$desc, s3)] = 3
    pos[which.max(df$count)] = 2
    # pos[df$count == min.count] = 2

    df$str = df$desc


    # for (i in 1:nchar(x)) print(paste(substr(x,i,i), substr(x,i,i) == substr(y,i,i)))
    xcolor = 0.75
    xlim = range(df$count)
    xlim[1] = min.count
    xlim[2] = xlim[2]
    ylim = range(df$enrichment)
    fig.start(fdir=fdir, ofn=paste(fdir, "/GO_", type, "_diagram.pdf", sep=""), type="pdf", width=10, height=6)
    par(xaxs="i")
    par(yaxs="i")
    plot.init(xlim=xlim, ylim=ylim, xlab="gene count (log scale)", ylab="enrichment (log2)", log="x", axis.las=1, y.axis=F, x.axis=F, add.grid=F, main=type)
    # at.y = axTicks(2)
    at.x = unique(sort(c(1:4, 1:5*2, 1:10*10)))
    at.y = -10:10
    axis(1, at=at.x)
    axis(2, at=at.y, las=1)

    rect(xleft=0.01, xright=xlim[2]*100, ybottom=0, ytop=ylim[2]+1, col=rgb(1,xcolor,xcolor), border=NA)
    rect(xleft=0.01, xright=xlim[2]*100, ybottom=ylim[1]-1, ytop=0, col=rgb(xcolor,xcolor,1), border=NA)
    abline(h=0, col="lightgray")
    abline(h=at.y, v=at.x, col="darkgray", lty=3)
    box()
    points(df$count, df$enrichment, col="lightgray", pch=19, cex=1.1)
    points(df$count, df$enrichment, col=1, cex=1.1)
    text(df$count, df$enrichment, df$str, pos=pos, cex=0.75)
    fig.end()
}

plot.GO=function(ifn, min.count, fdir)
{
    df = load.table(ifn)
    plot.GO.internal(df=df, type="func", fdir=fdir, min.count=min.count)
    plot.GO.internal(df=df, type="component", fdir=fdir, min.count=min.count)
    plot.GO.internal(df=df, type="process", fdir=fdir, min.count=min.count)
}
