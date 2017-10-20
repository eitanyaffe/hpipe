source("R/distrib.r")

# directory with files which have matching R1/R2 in names
list.file.pairs=function(path, pattern)
{
    fns = list.files(path, pattern)
    if (length(fns) == 0)
        stop(paste("no files found in", path))
    fns = fns[grepl("R[12]", fns)]
    fns = unique(gsub("R[12]", "@@@@", fns))
    bname = gsub("@@@@", "", fns)
    fns = paste(path, "/", fns, sep="")
    data.frame(base=bname, fn1=gsub("@@@@", "R1", fns), fn2=gsub("@@@@", "R2", fns))
}


distrib.sickle=function(
	  idir, odir, qsub.dir, sickle, type,
	  batch.max.jobs, total.max.jobs.fn, dtype, jobname)
{
  system(paste("rm -rf", qsub.dir))
  system(paste("mkdir -p", qsub.dir, odir))

  files = list.file.pairs(path=idir, pattern="*fastq*")

  N = dim(files)[1]
  cat(sprintf("processing %s paired files found in directory: %s\n", N, idir))

  commands = NULL
  for (i  in 1:N) {
      ifn1 = files$fn1[i]
      ifn2 = files$fn2[i]

      ofn1 = paste(odir, "/R1", files$base[i], sep="")
      ofn2 = paste(odir, "/R2", files$base[i], sep="")
      ofns = paste(odir, "/RS", files$base[i], sep="")

      log.file = paste(odir, "/.log", files$base[i], sep="")

      command = sprintf("%s pe -f  %s -r %s -t %s -o %s -p %s -s %s > %s 2>&1",
          sickle, ifn1, ifn2, type, ofn1, ofn2, ofns, log.file)
      # cat(sprintf("command: %s\n", command))
      commands = c(commands, command)
  }
  distrib(commands, working.dir=getwd(), qsub.dir=qsub.dir, jobname=jobname, batch.max.jobs=batch.max.jobs,
	  total.max.jobs.fn=total.max.jobs.fn, sleep.time=10, rprofile=NULL,
	  path=NULL, host.keys=NULL, retries=3, req.mem=NA, dtype=dtype)
}

