#########################################################################################################
# Plot matrix : utility function
#########################################################################################################

smatrix2matrix=function(smatrix, dim, i.field="i", j.field="j", value.field="value", default.value=0)
{
  indices = smatrix[,i.field] + (smatrix[,j.field]-1) * dim[1]
  v = rep(default.value, dim[1]*dim[2])
  v[indices] = smatrix[,value.field]
  matrix(v, dim[1], dim[2])
}

matrix2smatrix=function(matrix)
{
  dim1 = dim(matrix)[1]
  dim2 = dim(matrix)[2]
  v = as.vector(matrix)
  indices = 1:(dim1*dim2)
  i = (indices-1) %% dim1 + 1
  j = floor((indices-1) / dim1) + 1
  data.frame(i=i, j=j, value=v, stringsAsFactors=F)
}

make.color.panel=function(colors, ncols=256)
{
  panel = NULL
  for (i in 2:length(colors))
    panel = c(panel, colorpanel(ncols, colors[i-1], colors[i]))
  panel
}

# for example vals.to.cols(1:10, c(1, 3, 10), ncols=10) returns:
# [1] 1  6 11 12 14 15 16 17 19 20
vals.to.cols=function(vals, breaks, ncols=256)
{
  min = breaks[1]
  max = breaks[length(breaks)]
  n = length(breaks)-1
  cols = rep(-1, length(vals))
  for (i in 1:n)
  {
    ind = (breaks[i] <= vals) & (vals <= breaks[i+1])
    if (!any(ind))
      next
    # normalize to [0,1]
    cols[ind] = (vals[ind] - breaks[i]) / (breaks[i+1] - breaks[i])
    # normalize to [i*ncols,i*(ncols+1)]
    cols[ind] = (i-1)*ncols + cols[ind]*(ncols-1) + 1
    # round
    cols[ind] = round(cols[ind])
  }
  return (cols)
}

plot.colored.matrix=function(x, y, values,
                            size=0.5, breaks=NULL, ncols=256, border=NA, auto.axis=F, labels=NA, at=NA,
                            colors=c("blue", "cyan", "green", "yellow", "orange", "red"), rotated=F, plot.middle.axes=F, main=NULL, text=F, offset=0)
{
  if (!any(is.na(at))) {
    axis(1, at=at, labels=labels, las=2)
    axis(2, at=at, labels=labels, las=2)
  } else if (auto.axis) {
    axis(1)
    axis(2)
  }
  title(main=main)

  # temporarily put values instead of NAs
  na.ind = is.na(values)
  values[na.ind] = min(c(na.omit(values), 0))
  if (is.null(breaks))  {
    breaks = seq(min(values), max(values), length.out=length(colors))
    cat(sprintf("value range: (%f,%f)\n", min(values), max(values)))
  }
  min.value = breaks[1]
  max.value = breaks[length(breaks)]
  if (length(breaks) != length(colors))
    stop("#breaks and #colors must be equal")
  values[na.ind] = min.value

  # report on trimmed values
  if (min(values) < min.value || max(values) > max.value)
  {
      cat(sprintf("Warning: trimming range (%.1f,%.1f) into range (%.1f, %.1f)\n",
                  min(values), max(values), min.value, max.value))
      indices = which(values > max.value | values < min.value)
      if (length(indices) < 10) {
          for (ind in indices)
              cat(sprintf("  (%.1f, %.1f): %.1f\n", x[ind], y[ind], values[ind]))
      } else {
          cat(sprintf("range: (%.1f,%.1f)\n", min(values), max(values)))
      }
  } else {
      cat(sprintf("range: (%.1f,%.1f)\n", min(values), max(values)))
  }

  # trim values
  values[values < min.value] = min.value
  values[values > max.value] = max.value

  # make color panel
  panel = make.color.panel(colors, ncols)

  # get color indices
  col = vals.to.cols(values, breaks, ncols)

  if (!all(is.finite(col)))
    stop("missing values")

  # set colors
  rect.col = panel[col]

  # set NA colors
  rect.col[na.ind] = "lightgrey"

  # for middle axes
  mid = mean(x)

  if (!rotated) {
    if (plot.middle.axes) {
      abline(h=mid, lty=3)
      abline(v=mid, lty=3)
    }
    rect(x-size, y-size, x+size, y+size, col=rect.col, border=border)

    if (text)
      text(x, y, round(values,2))
  } else {
    s2 = 2
    xleft = x-size
    xright = x+size
    ybottom = y-size
    ytop = y+size
    xx = NULL
    yy = NULL
    polcol = NULL

    for (i in 1:length(xleft)) {
      if (xleft[i] < ybottom[i])
        next
      xx = c(xx, NA, (xleft[i]+ybottom[i])/s2, (xleft[i]+ytop[i])/s2, (xright[i]+ytop[i])/s2, (xright[i]+ybottom[i])/s2)
      yy = c(yy, NA, (xleft[i]-ybottom[i])/s2, (xleft[i]-ytop[i])/s2, (xright[i]-ytop[i])/s2, (xright[i]-ybottom[i])/s2)
      polcol = c(polcol, rect.col[i])
    }
    if (plot.middle.axes) {
      abline(a=mid, b=-1, lty=3)
      abline(a=-mid, b=1, lty=3)
    }
    polygon(x=xx, y=yy, col=polcol, border=NA)
  }

}

plot.colored.matrix.legend=function(values=1000*1:10, breaks=c(2,7,10)*1000, ncols=256, equal.spaced=F,
                                    colors=c("green", "red", "blue"),
                                    width=0.1, height=2, rotate=F)
{
  if (is.null(breaks))
    breaks = seq(min(na.omit(values)), max(na.omit(values)), length.out=length(colors))

  min.value = breaks[1]
  max.value = breaks[length(breaks)]
  if (length(breaks) != length(colors))
    stop("#breaks and #colors must be equal")

  # make color panel
  panel = make.color.panel(colors, ncols)

  nlegend = 100
  legend.values = seq(min.value, max.value, (max.value - min.value) / (nlegend-1))
  legend.cols = vals.to.cols(legend.values, breaks, ncols)

  par(mai=c(0.1, 0.1, 0.1, 0.1))
  plot.new()
  fin = par("fin")

  if (rotate) {
    plot.width = fin[2]
    plot.height = fin[1]
  } else {
    plot.width = fin[1]
    plot.height = fin[2]
  }

  # x axis
  xlim=c(0,1)
  offset.x = ((plot.width - width)/2) / plot.width
  if (offset.x < 0)
    stop("legend too wide")
  left = offset.x
  right = 1 - offset.x

  # y axis
  ylim = c(min.value, max.value)
  delta.y = ylim[2]-ylim[1]
  ratio = delta.y / height
  offset.y = ratio * (((plot.height - height)/2))
  if (offset.y < 0)
    stop("legend too high")
  ylim[1] = ylim[1] - offset.y
  ylim[2] = ylim[2] + offset.y

  if (rotate)
    plot.window(xlim=ylim, ylim=xlim)
  else
    plot.window(xlim=xlim, ylim=ylim)

  rect.col = panel[legend.cols]

  y = legend.values
  y.bottom = y[1:(length(y)-1)]
  y.top = y[2:length(y)]

  if (rotate)
    rect(y.bottom, left, y.top, right, col=rect.col, border=NA)
  else
    rect(left, y.bottom, right, y.top, col=rect.col, border=NA)

  if (equal.spaced)
    axis.values = seq(min.value, max.value,length.out=4)
  else
    axis.values = breaks
  tick.size = 0.07 / plot.width

  if (rotate) {
    segments(axis.values, right, axis.values, right+tick.size)
    text(axis.values, right+tick.size, axis.values, pos=3, cex=0.75)
  } else {
    segments(right, axis.values, right+tick.size, axis.values)
    text(right+tick.size, axis.values, axis.values, pos=4, cex=0.75)
  }
}

# wrapper function
plot.smatrix.simple=function(smat, size=0.5, visitor.f=NULL, title="",
                             ncols=256, border=NA, auto.axis=F, labels=NA, at=NA, equal.spaced=F,
                             breaks=NULL, colors=c("blue", "cyan", "green", "yellow", "orange", "red"), rotated=F, plot.middle.axes=F)
{
  par(xaxs="i")
  par(yaxs="i")

  plot.new()
  title(title)

  x = smat[,1]
  y = smat[,2]
  values = smat[,3]
  xlim=range(c(x+size, x-size))
  ylim=range(c(y+size, y-size))

  plot.window(xlim=xlim, ylim=ylim)
  plot.colored.matrix(x=x, y=y, size=size, auto.axis=auto.axis, labels=labels, at=at,
                      values=values, breaks=breaks, ncols=ncols, colors=colors, border=border,
                      rotated=rotated, plot.middle.axes=plot.middle.axes)
  if (!is.null(visitor.f))
    visitor.f()
}

# wrapper function
plot.smatrix=function(smat, size=0.5, visitor.f=NULL,
                      title="", fn="", width=5, height=4,
                      mai.matrix=c(0.7, 0.5, 0.4, 0.1),
                      mai.legend=c(0.5, 0.1, 0.5, 0.6),
                      legend.width=1.2,
                      breaks=NULL, ncols=256, border=NA, auto.axis=F, labels=NA, at=NA, equal.spaced=F,
                      colors=c("blue", "cyan", "green", "yellow", "orange", "red"), discard.diag=F)
{
  par(xaxs="i")
  par(yaxs="i")

  start.plot(fn, width=width, height=height, device="postscript")
  layout(matrix(c(1,1,2,2),2), width=c(width-legend.width,legend.width))
  op = par(mai=mai.matrix)
  plot.new()
  title(title)

  x = smat[,1]
  y = smat[,2]
  values = smat[,3]
  xlim=range(c(x+size, x-size))
  ylim=range(c(y+size, y-size))

  plot.window(xlim=xlim, ylim=ylim)

  plot.colored.matrix(x=x, y=y, size=size, auto.axis=auto.axis, labels=labels, at=at,
                      values=values, breaks=breaks, ncols=ncols, colors=colors, border=border)
  if (!is.null(visitor.f))
    visitor.f()
  par(mai=mai.legend)
  plot.colored.matrix.legend(values=values, breaks=breaks, ncols=ncols, colors=colors, equal.spaced=equal.spaced)
  end.plot(fn)
  par(op)
}
