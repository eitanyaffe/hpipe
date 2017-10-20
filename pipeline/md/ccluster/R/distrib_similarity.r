source("R/distrib.r")

distrib.similarity=function(
	  mummer, mummer.parse, split.dir, odir, 
	  qsub.dir, batch.max.jobs, total.max.jobs.fn, dtype, show.coords, parse.show.coords.script, jobname)
{
  system(paste("rm -rf", qsub.dir))
  system(paste("mkdir -p", qsub.dir, odir))

  files = list.files(split.dir)

  N = length(files)
  cat(sprintf("Using %s contig files found in directory: %s\n", N, split.dir))

  commands = NULL
  for (i  in 1:N) {
  for (j  in 1:N) {
    if (i == j) next
    file.i =  paste(split.dir, "/", i, sep="")
    file.j =  paste(split.dir, "/", j, sep="")

    if (!file.exists(file.i) || !file.exists(file.j))
      stop("input file missing")

    ofile =  paste(odir, "/", i, "_", j, ".table", sep="")
    command = paste(mummer, "-maxmatch -b -c -F", file.i, file.j, "|", mummer.parse, ofile)
    commands = c(commands, command)
  } }
  distrib(commands, working.dir=getwd(), qsub.dir=qsub.dir, jobname=jobname, batch.max.jobs=batch.max.jobs, 
	  total.max.jobs.fn=total.max.jobs.fn, sleep.time=10, rprofile=NULL, 
	  path=NULL, host.keys=NULL, retries=3, req.mem=NA, dtype=dtype)
}

