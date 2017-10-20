
is.defined=function(panel, field)
{
    is.element(field, names(panel))
}

get.field=function(panel, field, default)
{
    if (is.element(field, names(panel)))
        return (panel[[field]])
    default
}

plot.bar.f=function(panel)
{
    props = panel$props
    coords = 1:props$N - 0.5
    width = 0.8
    xleft = coords-width/2
    xright = coords+width/2

    if(!is.defined(panel, "log")) {
        values = props$df[,panel$field]
    } else {
        values = log10(1+props$df[,panel$field])
    }
        vrange = range(c(0,values))

    ylim = get.field(panel, "ylim", vrange)

    par(mai=props$mai)
    plot.init(xlim=props$xlim, ylim=ylim, ylab=panel$ylab, x.axis=F, y.axis=F, axis.las=2)
    rect(xleft=xleft, xright=xright, ybottom=0, ytop=values, col=props$bar.color, border=NA)

    if (is.element("field.add", names(panel))) {
        if(!is.defined(panel, "log")) {
            values.add = props$df[,panel$field.add]
        } else {
            values.add = log10(1+props$df[,panel$field.add])
        }
        rect(xleft=xleft, xright=xright, ybottom=0, ytop=values.add, col=props$bar.color.add, border=NA)
    }

    if(is.defined(panel, "log")) {
        at = 0:ceiling(max(ylim))
        labels = 10^at
    } else {
        at = axTicks(side=2)
        labels = at
    }

    axis(side=2, at=at, labels=labels, las=2)
}

plot.bar.range.f=function(panel)
{
    props = panel$props
    coords = 1:props$N - 0.5
    width = 0.8
    xleft = coords-width/2
    xright = coords+width/2

    values.bottom = props$df[,panel$field.bottom]
    values.top = props$df[,panel$field.top]
    vrange = range(c(0, values.bottom, values.top))
    ylim = get.field(panel, "ylim", vrange)

    par(mai=props$mai)
    plot.init(xlim=props$xlim, ylim=ylim, ylab=panel$ylab, x.axis=F, axis.las=2)
    rect(xleft=xleft, xright=xright, ybottom=values.bottom, ytop=values.top, col=props$bar.color, border=NA)
}

plot.ids.f=function(props)
{
    coords = 1:props$N - 0.5
    mai = props$mai
    par(mai=c(0, mai[2], 0, mai[4]))
    plot.new()
    plot.window(xlim=props$xlim, ylim=c(0, 1))
    text(x=coords, y=1, labels=props$id, srt=90, adj=1)
}

plot.elements=function(ifn.elements, fdir)
{
    df = load.table(ifn.elements)
    ids = df$id
    N = length(ids)
    xlim = c(0, N)
    props = list(df=df, xlim=xlim, N=N, bar.color="blue", bar.color.add="red", mai=c(0.1, 1.4, 0.1, 0.4), ids=ids)
    panels = list(
        host.count=list(props=props, ylab="#hosts", field="host.count", plot.f=plot.bar.f, height=1),
        identity.median=list(props=props, ylab="%diameter", field="identity.diameter", plot.f=plot.bar.f, height=1, ylim=c(40,80)),
        gene.count=list(props=props, ylab="#genes", field="gene.count", plot.f=plot.bar.f, height=1, log=T),
        novel.score=list(props=props, ylab="%unknown", field="novel.score", plot.f=plot.bar.f, height=1))

    M = length(panels)
    heights = NULL
    for (i in 1:M)
        heights = c(heights, panels[[i]]$height)
    heights = c(heights, 0.25)
    fig.start(fdir=fdir, ofn=paste(fdir, "/elements.pdf", sep=""), type="pdf", height=4, width=12)

    layout(matrix(1:(M+1), M+1, 1), heights=heights)
    for (i in 1:M) {
        panel = panels[[i]]
        cat(sprintf("plotting panel: %s\n", names(panels)[i]))
        panel$plot.f(panel)
    }

    plot.ids.f(props)
    fig.end()
}
