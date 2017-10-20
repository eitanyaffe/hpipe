source("R/distrib.r")

###########################################################
# fragment input
###########################################################

distrib.fragment=function(
    script,
    table.ifn,
    ref.dir,
    anchor.dir,
    read.length,
    step,
    reads.per.file,
    fragment.dir,
    qsub.dir,
    nthreads,
    batch.max.jobs,
    total.max.jobs.fn,
    dtype,
    jobname)
{
  system(paste("rm -rf", qsub.dir))
  system(paste("mkdir -p", qsub.dir))

  table = load.table(table.ifn)
  anchors = sort(unique(table$anchor))
  ids = sort(unique(table$accession))

  commands = NULL

  type = "file"
  ext = "NA"
  format = "fasta"

  for (id in ids) {
      input = paste(ref.dir, "/", id, "/genome.fasta", sep="")
      odir = paste(fragment.dir, "/refs/", id, sep="")
      command = paste("perl", script, type, input, ext, format, read.length, step, reads.per.file, 0, odir)
      # cat(sprintf("command: %s\n", command))
      commands = c(commands, command)
  }

  for (anchor in anchors) {
      input = paste(anchor.dir, "/", anchor, ".fasta", sep="")
      odir = paste(fragment.dir, "/anchors/", anchor, sep="")
      command = paste("perl", script, type, input, ext, format, read.length, step, reads.per.file, 0, odir)
      # cat(sprintf("command: %s\n", command))
      commands = c(commands, command)
  }
  distrib(commands, working.dir=getwd(), qsub.dir=qsub.dir, jobname=jobname, batch.max.jobs=batch.max.jobs,
	  total.max.jobs.fn=total.max.jobs.fn, sleep.time=10, rprofile=NULL,
	  path=NULL, host.keys=NULL, retries=3, req.mem=NA, dtype=dtype)
}

###########################################################
# bwa index
###########################################################

distrib.bwa.index=function(
    table.ifn,
    ref.dir,
    anchor.dir,
    bwa.binary,
    qsub.dir,
    nthreads,
    batch.max.jobs,
    total.max.jobs.fn,
    dtype,
    jobname)
{
  system(paste("rm -rf", qsub.dir))
  system(paste("mkdir -p", qsub.dir))

  table = load.table(table.ifn)
  anchors = sort(unique(table$anchor))
  ids = sort(unique(table$accession))

  commands = NULL

  # add refs
  for (id in ids) {
      fn = paste(ref.dir, "/", id, "/genome.fasta", sep="")
      command = paste(bwa.binary, "index", fn)
      # cat(sprintf("command: %s\n", command))
      commands = c(commands, command)
  }

  for (anchor in anchors) {
      fn = paste(anchor.dir, "/", anchor, ".fasta", sep="")
      command = paste(bwa.binary, "index", fn)
      # cat(sprintf("command: %s\n", command))
      commands = c(commands, command)
  }

  distrib(commands, working.dir=getwd(), qsub.dir=qsub.dir, jobname=jobname, batch.max.jobs=batch.max.jobs,
	  total.max.jobs.fn=total.max.jobs.fn, sleep.time=10, rprofile=NULL,
	  path=NULL, host.keys=NULL, retries=3, req.mem=NA, dtype=dtype)
}

###########################################################
# mapping
###########################################################

distrib.ref2anchor=function(
    script,
    table.ifn,
    ref.dir,
    anchor.dir,
    bwa.binary,
    parse.bwa,
    profile.binary,
    read.length,
    fragment.dir,
    odir,
    qsub.dir,
    nthreads,
    batch.max.jobs,
    total.max.jobs.fn,
    dtype,
    jobname)
{
    cat(sprintf("clearing output directory: %s\n", odir))
    system(paste("rm -rf", qsub.dir, odir))
    system(paste("mkdir -p", qsub.dir, odir))

    table = load.table(table.ifn)

    commands = NULL
    for (i  in 1:dim(table)[1]) {
        anchor = table$anchor[i]
        accession = table$accession[i]

        src.contigs = paste(ref.dir, "/", accession, "/contig_table", sep="")
        tgt.contigs = paste(anchor.dir, "/", anchor, ".contig_table", sep="")

        idir = paste(fragment.dir, "/refs/", accession, sep="")
        target_fasta = paste(anchor.dir, "/", anchor, ".fasta", sep="")
        tdir = paste(odir, "/", anchor, "_", accession, sep="")
        bwa.index.prefix = target_fasta

        parse.src = "T"

        command = paste("perl",
            script, idir, target_fasta, bwa.binary, parse.bwa,
          profile.binary, src.contigs, tgt.contigs, read.length, parse.src,
            nthreads, bwa.index.prefix, tdir)
        # cat(sprintf("command: %s\n", command))
        commands = c(commands, command)
    }
    distrib(commands, working.dir=getwd(), qsub.dir=qsub.dir, jobname=jobname, batch.max.jobs=batch.max.jobs,
            total.max.jobs.fn=total.max.jobs.fn, sleep.time=10, rprofile=NULL,
            path=NULL, host.keys=NULL, retries=3, req.mem=NA, dtype=dtype)
}

distrib.anchor2ref=function(
    script,
    table.ifn,
    ref.dir,
    anchor.dir,
    bwa.binary,
    parse.bwa,
    profile.binary,
    read.length,
    fragment.dir,
    odir,
    qsub.dir,
    nthreads,
    batch.max.jobs,
    total.max.jobs.fn,
    dtype,
    jobname)
{
    system(paste("rm -rf", qsub.dir, odir))
    system(paste("mkdir -p", qsub.dir, odir))

    table = load.table(table.ifn)

    commands = NULL
    for (i  in 1:dim(table)[1]) {
        anchor = table$anchor[i]
        accession = table$accession[i]

        src.contigs = paste(anchor.dir, "/", anchor, ".contig_table", sep="")
        tgt.contigs = paste(ref.dir, "/", accession, "/contig_table", sep="")

        idir = paste(fragment.dir, "/anchors/", anchor, sep="")
        target_fasta = paste(ref.dir, "/", accession, "/genome.fasta", sep="")
        tdir = paste(odir, "/", anchor, "_", accession, sep="")
        bwa.index.prefix = target_fasta

        parse.src = "T"

        command = paste("perl",
            script, idir, target_fasta, bwa.binary, parse.bwa,
            profile.binary, src.contigs, tgt.contigs, read.length, parse.src,
            nthreads, bwa.index.prefix, tdir)
        # cat(sprintf("command: %s\n", command))
        commands = c(commands, command)
    }
    distrib(commands, working.dir=getwd(), qsub.dir=qsub.dir, jobname=jobname, batch.max.jobs=batch.max.jobs,
            total.max.jobs.fn=total.max.jobs.fn, sleep.time=10, rprofile=NULL,
            path=NULL, host.keys=NULL, retries=3, req.mem=NA, dtype=dtype)
}

distrib.read2ref=function(
    script,
    table.ifn,
    ref.dir,
    bwa.binary,
    parse.bwa,
    profile.binary,
    read.length,
    idir,
    odir,
    qsub.dir,
    nthreads,
    batch.max.jobs,
    total.max.jobs.fn,
    dtype,
    jobname)
{
    system(paste("rm -rf", qsub.dir, odir))
    system(paste("mkdir -p", qsub.dir, odir))

    table = load.table(table.ifn)

    commands = NULL
    ids = sort(unique(table$accession))

    commands = NULL
    for (accession in ids) {
        src.contigs = "NA"
        tgt.contigs = paste(ref.dir, "/", accession, "/contig_table", sep="")

        target_fasta = paste(ref.dir, "/", accession, "/genome.fasta", sep="")
        tdir = paste(odir, "/", accession, sep="")
        bwa.index.prefix = target_fasta

        parse.src = "F"

        command = paste("perl",
            script, idir, target_fasta, bwa.binary, parse.bwa,
            profile.binary, src.contigs, tgt.contigs, read.length, parse.src,
            nthreads, bwa.index.prefix, tdir)
        # cat(sprintf("command: %s\n", command))
        commands = c(commands, command)
    }
    distrib(commands, working.dir=getwd(), qsub.dir=qsub.dir, jobname=jobname, batch.max.jobs=batch.max.jobs,
            total.max.jobs.fn=total.max.jobs.fn, sleep.time=10, rprofile=NULL,
            path=NULL, host.keys=NULL, retries=3, req.mem=NA, dtype=dtype)
}

###########################################################
# collect results
###########################################################

distrib.ref2anchor.table=function(
    binary,
    ref.dir,
    anchor.dir,
    table.ifn,
    wdir,
    length,
    qsub.dir,
    nthreads,
    batch.max.jobs,
    total.max.jobs.fn,
    dtype,
    jobname)
{
  system(paste("rm -rf", qsub.dir))
  system(paste("mkdir -p", qsub.dir))

  table = load.table(table.ifn)

  commands = NULL
  for (i  in 1:dim(table)[1]) {
      anchor = table$anchor[i]
      accession = table$accession[i]

      idir = paste(wdir, "/", anchor, "_", accession, "/parsed", sep="")
      src.contigs = paste(ref.dir, "/", accession, "/contig_table", sep="")
      tgt.contigs = paste(anchor.dir, "/", anchor, ".contig_table", sep="")
      src.ofn = paste(wdir, "/", anchor, "_", accession, "/src_table", sep="")
      tgt.ofn = paste(wdir, "/", anchor, "_", accession, "/tgt_table", sep="")
      src.summary.ofn = paste(wdir, "/", anchor, "_", accession, "/src_summary", sep="")
      tgt.summary.ofn = paste(wdir, "/", anchor, "_", accession, "/tgt_summary", sep="")

      command = sprintf("%s -parse_src T -idir %s -src_contig_table %s -tgt_contig_table %s -read_length %d -src_ofn %s -tgt_ofn %s -src_summary_ofn %s -tgt_summary_ofn %s",
          binary, idir, src.contigs, tgt.contigs, length, src.ofn, tgt.ofn, src.summary.ofn, tgt.summary.ofn);
      # cat(sprintf("command: %s\n", command))
      commands = c(commands, command)
  }
  distrib(commands, working.dir=getwd(), qsub.dir=qsub.dir, jobname=jobname, batch.max.jobs=batch.max.jobs,
	  total.max.jobs.fn=total.max.jobs.fn, sleep.time=10, rprofile=NULL,
	  path=NULL, host.keys=NULL, retries=3, req.mem=NA, dtype=dtype)
}

distrib.anchor2ref.table=function(
    binary,
    ref.dir,
    anchor.dir,
    table.ifn,
    wdir,
    length,
    qsub.dir,
    nthreads,
    batch.max.jobs,
    total.max.jobs.fn,
    dtype,
    jobname)
{
  system(paste("rm -rf", qsub.dir))
  system(paste("mkdir -p", qsub.dir))

  table = load.table(table.ifn)

  commands = NULL
  for (i  in 1:dim(table)[1]) {
      anchor = table$anchor[i]
      accession = table$accession[i]

      idir = paste(wdir, "/", anchor, "_", accession, "/parsed", sep="")
      src.contigs = paste(anchor.dir, "/", anchor, ".contig_table", sep="")
      tgt.contigs = paste(ref.dir, "/", accession, "/contig_table", sep="")
      src.ofn = paste(wdir, "/", anchor, "_", accession, "/src_table", sep="")
      tgt.ofn = paste(wdir, "/", anchor, "_", accession, "/tgt_table", sep="")
      src.summary.ofn = paste(wdir, "/", anchor, "_", accession, "/src_summary", sep="")
      tgt.summary.ofn = paste(wdir, "/", anchor, "_", accession, "/tgt_summary", sep="")

      command = sprintf("%s -parse_src T -idir %s -src_contig_table %s -tgt_contig_table %s -read_length %d -src_ofn %s -tgt_ofn %s -src_summary_ofn %s -tgt_summary_ofn %s",
          binary, idir, src.contigs, tgt.contigs, length, src.ofn, tgt.ofn, src.summary.ofn, tgt.summary.ofn);
      # cat(sprintf("command: %s\n", command))
      commands = c(commands, command)
  }
  distrib(commands, working.dir=getwd(), qsub.dir=qsub.dir, jobname=jobname, batch.max.jobs=batch.max.jobs,
	  total.max.jobs.fn=total.max.jobs.fn, sleep.time=10, rprofile=NULL,
	  path=NULL, host.keys=NULL, retries=3, req.mem=NA, dtype=dtype)
}

distrib.read2ref.table=function(
    binary,
    ref.dir,
    anchor.dir,
    table.ifn,
    wdir,
    length,
    qsub.dir,
    nthreads,
    batch.max.jobs,
    total.max.jobs.fn,
    dtype,
    jobname)
{
  system(paste("rm -rf", qsub.dir))
  system(paste("mkdir -p", qsub.dir))

  table = load.table(table.ifn)

  commands = NULL
  for (i  in 1:dim(table)[1]) {
      anchor = table$anchor[i]
      accession = table$accession[i]

      idir = paste(wdir, "/", accession, "/parsed", sep="")
      tgt.contigs = paste(ref.dir, "/", accession, "/contig_table", sep="")
      tgt.ofn = paste(wdir, "/", accession, "/tgt_table", sep="")
      tgt.summary.ofn = paste(wdir, "/", accession, "/tgt_summary", sep="")

      command = sprintf("%s -parse_src F -idir %s -tgt_contig_table %s -read_length %d %s -tgt_ofn %s%s -tgt_summary_ofn %s",
          binary, idir, tgt.contigs, length, tgt.ofn, tgt.summary.ofn);
      #cat(sprintf("command: %s\n", command))
      commands = c(commands, command)
  }
  distrib(commands, working.dir=getwd(), qsub.dir=qsub.dir, jobname=jobname, batch.max.jobs=batch.max.jobs,
	  total.max.jobs.fn=total.max.jobs.fn, sleep.time=10, rprofile=NULL,
	  path=NULL, host.keys=NULL, retries=3, req.mem=NA, dtype=dtype)
}
