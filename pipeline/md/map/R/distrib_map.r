source("R/distrib.r")

# directory with files which have matching R1/R2 in names
list.file.pairs=function(path, pattern)
{
    fns = list.files(path, pattern)
    if (length(fns) == 0)
        stop(paste("no files found in", path))
    if (!all(grepl("R[12]", fns)))
        stop("not all files have R[12] in name")
    fns = unique(gsub("R[12]", "@@@@", fns))
    bname = gsub("@@@@", "", fns)
    fns = paste(path, "/", fns, sep="")
    data.frame(base=bname, fn1=gsub("@@@@", "R1", fns), fn2=gsub("@@@@", "R2", fns))
}


distrib.filter=function(
	  script, idir, odir, sdir, sfile, qsub.dir, min.score, min.length, min.distance,
	  batch.max.jobs, total.max.jobs.fn, dtype, jobname)
{
  system(paste("rm -rf", qsub.dir))
  system(paste("mkdir -p", qsub.dir, odir))

  files = list.files(idir, pattern="*q")

  N = length(files)
  if (N == 0)
      stop("no files found")
  cat(sprintf("processing %s files found in directory: %s\n", N, idir))

  commands = NULL
  for (i  in 1:N) {
      ifn = paste(idir, "/", files[i], sep="")
      ofn = paste(odir, "/", files[i], sep="")
      sfn = paste(sdir, "/", files[i], sep="")

      command = sprintf("perl %s %s %d %d %d %s %s", script, ifn, min.score, min.length, min.distance, ofn, sfn)
      commands = c(commands, command)
  }
  distrib(commands, working.dir=getwd(), qsub.dir=qsub.dir, jobname=jobname, batch.max.jobs=batch.max.jobs,
	  total.max.jobs.fn=total.max.jobs.fn, sleep.time=10, rprofile=NULL,
	  path=NULL, host.keys=NULL, retries=3, req.mem=NA, dtype=dtype)

  # unite stat files
  result = NULL
  for (i  in 1:N) {
      sfn = paste(sdir, "/", files[i], sep="")
      x = read.delim(sfn)

      if (i == 1) {
          result = x
      } else {
          result = result + x
      }
  }
  save.table(result, sfile)
}

distrib.pair=function(
	  script, idir, odir, sdir, sfile, qsub.dir, batch.max.jobs, total.max.jobs.fn, dtype, jobname)
{
  system(paste("rm -rf", qsub.dir))
  system(paste("mkdir -p", qsub.dir, odir))

  files = list.file.pairs(idir, pattern="*q")
  N = dim(files)[1]
  if (N == 0)
      stop("no files found")

  cat(sprintf("processing %s files found in directory: %s\n", N, idir))

  commands = NULL
  for (i  in 1:N) {
      ifn1 = files$fn1[i]
      ifn2 = files$fn2[i]
      ofn = paste(odir, "/P", files$base[i], ".pair", sep="")
      sfn = paste(sdir, "/P", files$base[i], ".pair", sep="")

      command = sprintf("perl %s %s %s %s %s", script, ifn1, ifn2, ofn, sfn)
      commands = c(commands, command)
  }
  distrib(commands, working.dir=getwd(), qsub.dir=qsub.dir, jobname=jobname, batch.max.jobs=batch.max.jobs,
	  total.max.jobs.fn=total.max.jobs.fn, sleep.time=10, rprofile=NULL,
	  path=NULL, host.keys=NULL, retries=3, req.mem=NA, dtype=dtype)


  # unite stat files
  result = NULL
  for (i  in 1:N) {
      sfn = paste(sdir, "/P", files$base[i], ".pair", sep="")
      x = read.delim(sfn)

      if (i == 1) {
          result = x
      } else {
          result = result + x
      }
  }
  save.table(result, sfile)
}
