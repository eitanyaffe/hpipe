source("R/distrib.r")

distrib.parse.bwa=function(
    script, idir, odir, qsub.dir, batch.max.jobs, total.max.jobs.fn, dtype, jobname)
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

      command = sprintf("perl %s %s %s", script, ifn, ofn)
#      cat(sprintf("command: %s\n", command))
      commands = c(commands, command)
  }
  distrib(commands, working.dir=getwd(), qsub.dir=qsub.dir, jobname=jobname, batch.max.jobs=batch.max.jobs,
	  total.max.jobs.fn=total.max.jobs.fn, sleep.time=10, rprofile=NULL,
	  path=NULL, host.keys=NULL, retries=3, req.mem=NA, dtype=dtype)
}
