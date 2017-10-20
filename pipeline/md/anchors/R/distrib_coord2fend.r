#!/bin/env Rscript

#! /net/mraid04/export/users/eitany/root/bin/Rscript --vanilla

options(warn=1)

# get script name
all.args = commandArgs(F)
fn.arg = "--file="
script.name = sub(fn.arg, "", all.args[grep(fn.arg, all.args)])

args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <title> <input dir> <fends table> <read length> <trim length> <facing threshold> <max segment length> <split dir> <output dir> <qsub directory> <working directory> <dtype>\n", script.name))
  q(status=1)
}

title = args[1]
input.dir = args[2]
fends.fn = args[3]
read.length = as.numeric(args[4])
trim.length = as.numeric(args[5])
facing.threshold = as.numeric(args[6])
max.seg.len = as.numeric(args[7])
split.dir = args[8]
output.dir = args[9]
qsub.dir = args[10]
wd = args[11]
dtype = args[12]
total.max.jobs.fn = args[13]
chr.size.fn = args[14]
correct.trimmed = args[15]

cat(sprintf("cleaning qsub dir: %s\n", qsub.dir))
system(paste("rm -rf", qsub.dir))
system(paste("mkdir -p", qsub.dir))

cat(sprintf("cleaning split dir: %s\n", split.dir))
system(paste("rm -rf", split.dir))
system(paste("mkdir -p", split.dir))

ifiles = list.files(input.dir)
if (length(ifiles) == 0) {
  stop(sprintf("no files in input directory: %s\n", input.dir))
} else {
  cat(sprintf("found %d input files in input directory: %s\n", length(ifiles), input.dir))
}

commands = NULL
for (i in seq_along(ifiles)) {
  oprefix.job = paste(split.dir, "/", i, sep="")
  file = paste(input.dir, ifiles[i], sep="/")
  command = sprintf("perl %s/pl/coords2fends.pl %s %d %d %d %d %s %s %s %s", wd, fends.fn, read.length, trim.length, facing.threshold, max.seg.len, chr.size.fn, oprefix.job, correct.trimmed, file)
  commands = c(commands, command)
}

time = proc.time()
source(paste(wd, "/R/distrib.r", sep=""))
distrib(commands, working.dir=wd, qsub.dir=qsub.dir, jobname=paste("coord2fends", title, sep="_"), batch.max.jobs=400, total.max.jobs.fn=total.max.jobs.fn, sleep.time=10, dtype=dtype)
cat(sprintf("distributed computation took %.1f seconds\n", (proc.time() - time)[3]))

# give time for nfs to update files
max.sleeps = 0
while (1) {
  command = sprintf("grep done %s/*.stderr | wc -l", qsub.dir)
  cat(sprintf("command: %s\n" ,command))
  done = as.numeric(system(command, intern=T))
  if (done == length(ifiles))
    break
  cat("All jobs not done yet, sleeping for 30 secs...\n")
  Sys.sleep(30)
  max.sleeps = max.sleeps + 1
  if (max.sleeps > 10)
    stop("not all distributed tasks finished successfully, stopping")
}

# unite general stats
all.fns = NULL
for (i in seq_along(ifiles)) {
  all.fns = c(all.fns, paste(split.dir, "/", i, "_all.mat.stats", sep=""))
}

unite.stats = function(ifns, ofn) {
  if (length(ifns) == 0)
    stop("no files to unite")
  for (i in seq_along(ifns)) {
    data = read.delim(ifns[i])

    # when uniting chroms might have different order
    if (i > 1) {
      if (!(all(result$chr == data$chr)))
	stop("different chrom order")
      result[,-1] = result[,-1] + data[,-1]
    } else
      result = data
  }
  write.table(result, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}

cat("merging stat files...\n")
unite.stats(all.fns, paste(output.dir, "/all.stats", sep=""))

# unite type specific stats
types = c("s0", "s1", "s2")
chr.size = read.delim(chr.size.fn, header=F)
for (type in types) {
  read.stats.fns = NULL
  for (i in seq_along(ifiles)) {
    oprefix.job = paste(split.dir, "/", i, "_", type, sep="")

    read.stats.fns = c(read.stats.fns, paste(oprefix.job, ".read.stats", sep=""))
  }
  unite.stats(read.stats.fns, paste(output.dir, "/", type, ".read.stats", sep=""))
}

cat("merging results...\n")
for (type in types) {
  oprefix = paste(output.dir, "/", type, sep="")
  command = sprintf("%s/pl/merge_mat_files.pl %s %s %s %s", wd, fends.fn, oprefix, split.dir, type)
  cat(sprintf("command: %s\n", command))
  if (system(command) != 0)
    stop(sprintf("error in command: %s", command))
}
