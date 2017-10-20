#!/bin/env Rscript

#! /net/mraid04/export/users/eitany/root/bin/Rscript --vanilla

options(warn=1)

cat(sprintf("current dir: %s\n", getwd()))

# get script name
all.args = commandArgs(F)
fn.arg = "--file="
script.name = sub(fn.arg, "", all.args[grep(fn.arg, all.args)])

args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <input prefix> <model file> <tmp dir> <log dir> <model_integrate binary file> <dtype> <mjobs> <fends suffix> <Rcall>\n",
              script.name))
  q(status=1)
}

ifn.prefix = args[1]
model.ifn = args[2]
tmp.dir = args[3]
log.dir = args[4]
binary = args[5]
wd = args[6]
dtype = args[7]
mjobs = args[8]
fends.suffix = args[9]
Rcall = args[10]

mtable = read.delim(model.ifn)
mfields = mtable$field
maxvals = mtable$size

source(paste(wd, "/R/model_predict.r", sep=""))
cat(sprintf("calling: source(\"%s/R/model_predict.r\")\n", wd))
cat(sprintf("calling: compute.total.counts(binary=\"%s\", tmp.dir=\"%s\", log.dir=\"%s\", prefix=\"%s\", ofields=c(%s), max.vals=c(%s), wd=%s, dtype=%s, Rcall=%s)\n",
            binary, tmp.dir, log.dir, ifn.prefix, paste(mfields, sep="", collapse=","), paste("\"", maxvals, "\"", sep="", collapse=","), wd, dtype, Rcall))

compute.total.counts(binary=binary, tmp.dir=tmp.dir, log.dir=log.dir, prefix=ifn.prefix, ofields=mfields, max.vals=maxvals,
		     wd=wd, dtype=dtype,  max.jobs.fn=mjobs, fends.suffix=fends.suffix, Rcall=Rcall)

q(status=0)
