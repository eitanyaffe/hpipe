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

distrib.remove.human=function(
    deconseq, dbs, identity, coverage, idir, odir, wdir,
    qsub.dir, batch.max.jobs, total.max.jobs.fn, dtype, jobname)
{
  system(paste("rm -rf", qsub.dir))
  system(paste("mkdir -p", qsub.dir, odir))

  files = list.files(path=idir, pattern="*fastq*")

  N = length(files)
  cat(sprintf("processing %s files found in directory: %s\n", N, idir))
  if (N == 0)
      stop()
  commands = NULL
  for (i  in 1:N) {
      ifn = paste(idir, "/", files[i], sep="")
      log.file = paste(odir, "/.log_", files[i], sep="")
      command = sprintf(
          "perl %s -c %d -i %d -f %s -id %s -dbs %s -out_dir %s > %s 2>&1",
          deconseq, coverage, identity, ifn, files[i], dbs, odir, log.file)

#      cat(sprintf("command: %s\n", command))
      commands = c(commands, command)
  }
  distrib(commands, working.dir=wdir, qsub.dir=qsub.dir, jobname=jobname, batch.max.jobs=batch.max.jobs,
	  total.max.jobs.fn=total.max.jobs.fn, sleep.time=10, rprofile=NULL,
	  path=NULL, host.keys=NULL, retries=3, req.mem=NA, dtype=dtype)
}

