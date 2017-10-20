source("R/colored_matrix.r")
source("R/single.r")
library(gplots)

plot.feature.matrices = function(
  idir="/home/lubling/storage14/datasets/sexton_mouse_hic/ESA/model/",
  dataset="ESA", fdir="/home/lubling/storage14/datasets/sexton_mouse_hic/ESA/figures/model/",
  model.fn="/home/lubling/storage14/o3c_pipeline/models/map_len_gc.mdl",
  wd="/home/lubling/storage14/o3c_pipeline/", width=4, height=4,
  breaks=c(-4, -2, -1, 0, 0.5, 1), break.cols=c("black", "blue", "cyan", "white", "orange", "red"))
{
  odir.dataset = paste(fdir, "/", dataset, sep="")
  odir.summary = paste(fdir, "/summary", sep="")

  system(paste("mkdir -p", odir.dataset))
  system(paste("mkdir -p", odir.summary))

  size = 0.5
  units = "in"

  #border = "white"
  border = NA

  mf = read.delim(model.fn)
  for (seed in c(T, F)) {
    for (f in mf$raw_field) {


      ifn = paste(idir, "/main_", f, if (seed ) "_bin_seed.f" else "_bin.f", sep="")
      if (!file.exists(ifn)) {
          cat(sprintf("skipping file: %s\n", ifn))
          next
      }
      vals = read.delim(ifn)
      vals[,3] = vals[,3] / median(vals[,3])

      ranges = read.delim(paste(idir, "/main_", f, ".bin_ranges", sep=""))

      m = matrix(0, nrow(ranges), nrow(ranges))
      rownames(m) = 1:nrow(ranges)
      colnames(m) = 1:nrow(ranges)

      m[ cbind(vals[,1], vals[,2]) ] = log10(vals[,3])
      m[ cbind(vals[,2], vals[,1]) ] = log10(vals[,3])
      smat = matrix2smatrix(m)

      main = if (seed) "seed" else "model"
      main = paste(dataset, f, main)
      rr = range(smat[,3])
      main = sprintf("%s min=%g, max=%g", f, rr[1], rr[2])
      lplot = function(fn) {
        cat(sprintf("generating figure: %s\n", fn))
        png(fn, width=width, height=height, units=units, res=100)

#      plot.colored.matrix(x=smat$i, y=smat$j, size=0.5, values=smat$value, breaks=breaks, xlabels=ranges$start, at.x=seq_along(ranges$start), ylabels=ranges$start, at.y=seq_along(ranges$start), colors=break.cols, colnote.from.num=NA, ylab=f, xlab=f, main=if (seed) "seed" else "model",  border="white")

        par(xaxs="i")
        par(yaxs="i")
        plot.new()
        title(main=main)
        x = smat[,1]
        y = smat[,2]
        values = smat[,3]
        xlim=range(c(x+size, x-size))
        ylim=range(c(y+size, y-size))
        plot.window(xlim=xlim, ylim=ylim)
        plot.colored.matrix(x=smat$i, y=smat$j, size=0.5, values=smat$value, breaks=breaks, labels=ranges$start, at=seq_along(ranges$start),
                            colors=break.cols, border=border)
        dev.off()
      }

      fn = paste(odir.dataset, "/", f, if (seed) "_seed" else "_model", ".png", sep="")
      lplot(fn)

      # plot only final matrices in summary dir
      if (!seed) {
        fn = paste(odir.summary, "/", f, "_", dataset, ".png", sep="")
        lplot(fn)
      }

      # legend
      mai.legend = c(0.5, 0.1, 0.5, 0.6)
      fn = paste(fdir, "/", if (seed) "seed" else "model", "_legend.png", sep="")
      cat(sprintf("generating figure: %s\n", fn))
      par(mai=mai.legend)
      png(fn, width=4, height=4, units=units, res=100)
        plot.colored.matrix.legend(breaks=breaks, colors=break.cols, height=2, rotate=T)
      dev.off()
    }
  }

  return(0)
}
