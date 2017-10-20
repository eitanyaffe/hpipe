##################################################################
# Source file managment
##################################################################

# remove all files
reset.source=function()
{
  g.source.files <<- list()
}

# register a file
register.source=function(file)
{
  if (!is.element(file, g.source.files))
    g.source.files <<- c(g.source.files, list(file))
  rl()
}

go=function()
{
  graphics.off()
}

# reload source files
rl=function()
{
  #X11.options(type="Xlib")
  X11.options(type="nbcairo")
  options(stringsAsFactors=F)
  options(warn=1)
  library(gplots)
  library(igraph)
  library(Rsge)
#  library(maxLik)

  for (f in g.source.files)
    source(f)
  
  #options(error=NULL)
  options(error=recover)

  g.current.plot.fn <<- ""
}

##################################################################
# Global cache
##################################################################

ca.clean=function(context=NULL)
{
  if (is.null(context))
    g.ca <<- list()
  else
    g.ca[[context]] <<- list()

  g.ca.context <<- context  
}

ca.clean.tag=function(tags)
{
  for (tag in tags)
    g.ca[[ca.get.context()]][[tag]] <<- NULL
}

ca.get=function(tag, value, verbose=F)
{
  context = g.ca.context
  if (is.null(context))
    context = "global"
  
  if (!exists("g.ca"))
    g.ca <<- list()
  if (is.null(g.ca[[context]]))
    g.ca[[context]] <<- list()
  
  if (is.null(g.ca[[context]][[tag]]))
  {
    if (verbose) cat(sprintf("ca.get: tag=%s\n", tag))
    tmp = value
    g.ca[[context]][[tag]] <<- tmp
  }
  
  return (g.ca[[context]][[tag]])
}

# shortcut to read tables
ca.get.table=function(table.fn, verbose=F)
{
  return (ca.get(table.fn, read.delim(table.fn), verbose=verbose))
}

ca.read=function(tag)
{
  context = g.ca.context
  if (is.null(context))
    context = "global"
  return (g.ca[[context]][[tag]])
}

ca.context=function(context=NULL, clean=F)
{
  if (clean)
    ca.clean(context)
  g.ca.context <<- context
}

ca.get.context=function()
{
  if (is.null(g.ca.context))
    return ("global")
  else
    return (g.ca.context)
}

ca.names=function(context = ca.get.context())
{
  if (is.null(context))
    context = "global"
  return (names(g.ca[[context]]))
}

ca.set.root=function(dir="/net/mraid14/export/data/db/tgdb/mm9/trackdb/", gmax.data.size=NA)
{
  if (!exists(".groot_set", where=.GlobalEnv)) {
    cat ("setting root...\n")
    gsetroot(dir)
    options(gparam.type="string")

    if (!is.na(gmax.data.size))
      options(gmax.data.size = gmax.data.size)
    
    assign(".groot_set", T, envir=.GlobalEnv)
  }
}

##################################################################
# Utility functions
##################################################################

assert=function (cond, msg="assert failed") 
{
  if (!cond) {
    .Internal(stop(as.logical(T), .makeMessage(msg)) )   
  }
}

get.fn.dir=function(fn)
{  
  grep = gregexpr("/", fn)[[1]]
  substr(fn, 1, grep[length(grep)])
}

start.plot=function(fn="", width=5, height=5, device="postscript", png.units="in", png.res=72)
{  
  if (fn != "") {
    cat(sprintf("plot start: %s\n", fn))
    if (device == "postscript") {
      postscript(file=fn, width=width, height=height, horizontal = FALSE, onefile = FALSE,
                 paper = "special", pointsize=10)
    } else if (device == "png") {
      png(file=fn, width=width, height=height, units=png.units, res=png.res)
    } else if (device == "bmp") {
      bmp(file=fn, width=width, height=height, units=png.units, res=png.res)
    } else stop("unknown device")

    g.current.plot.fn <<- fn
    g.current.plot.device <<- device
  }
}

end.plot=function(fn)
{
  if (fn != "")
  {  
    assert(g.current.plot.fn == fn, "cannot end plot in middle of other plot")
    dev.off()
    if (g.current.plot.device == "postscript") {
      # cat(sprintf("Converting: %s ... ", fn))
      # system(paste("convert", fn, sub("eps", "png", fn)))
      # cat("Done\n")
      #system(paste("eps2png", fn))
    }
    g.current.plot.fn <<- ""
    cat(sprintf("plot done: %s\n", fn))
  }
}

n2str=function(n)
{
  if (n >= 10^6)
    n.str = paste(n/10^6, "M", sep="")
  else
    n.str = paste(n/10^3, "K", sep="")

  n.str
}

between=function(a,b,x)
{
  x >= a & x <= b
}

interv.intersect=function(a1,a2,b1,b2)
{
  between(a1,a2,b1) | between(b1,b2,a1)
}
