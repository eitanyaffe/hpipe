#########################################################################################################
# Internal functions
#########################################################################################################

get.short.fn=function(fn)
{
  grep = gregexpr("/", fn)[[1]]
  substr(fn, grep[length(grep)]+1, nchar(fn))
}

get.ofield.args=function(ranges, fields, suffix)
{
  from = ranges[paste(fields, suffix, "from", sep="_")]
  to = ranges[paste(fields, suffix, "to", sep="_")]
  paste(length(fields), paste(fields, from, to, collapse=" "))
}

get.ofield.range=function(ranges, fields, suffix)
{
  from = ranges[paste(fields, suffix, "from", sep="_")]
  to = ranges[paste(fields, suffix, "to", sep="_")]
  list(from=from, to=to)
}

model.predict=function(
    ofield.x.from, ofield.x.to, ofield.y.from, ofield.y.to, index, from.fend, to.fend,
    fends.fn="/home/eitany/storage/o3c/results/h1_cl_gc_1000.binned",
    tmp.dir="/home/eitany/storage/o3c/tmp/", log.dir="/home/eitany/storage/o3c/log/",
    binary = "/home/eitany/storage/o3c/bin/model_integrate",
    mfields=NULL, mfields.maxvals=NULL,
    mfields.fns=NULL,
    scope, model, prior,
    ofields.x=c("f1", "f2"), ofields.y=c("f1", "f2"))
{
  if (is.null(ofields.x)) stop("ofields.x not defined")
  if (is.null(ofields.y)) stop("ofields.y not defined")

  ofn=paste(tmp.dir, "/", get.short.fn(fends.fn), ".", index, sep="")
  log=paste(log.dir, "/", get.short.fn(fends.fn), ".", index, sep="")

  # basic parameters
  command = sprintf("%s -fends %s -from_fend %d -to_fend %d -output %s", binary, fends.fn, from.fend, to.fend, ofn)

  # scope and model
  command = sprintf("%s -scope %s -model %s -prior %g", command, scope, model, prior)

  # model features
  length.mfields = if (!is.null(mfields)) length(mfields) else 0
  command = sprintf("%s -model_num %d", command, length.mfields)
  if (length.mfields > 0)
  for (i in 1:length.mfields)
    command = sprintf("%s -f_model_field %s -f_model_fn %s -f_model_size %d",
      command, mfields[i], mfields.fns[i], mfields.maxvals[i])

  # x output features
  command = sprintf("%s -output_x_num %d", command, length(ofields.x))
  for (i in 1:length(ofields.x))
    command = sprintf("%s -f_x_field %s -f_x_from %d -f_x_to %d",
      command, ofields.x[i], ofield.x.from[i], ofield.x.to[i])

  # y output features
  command = sprintf("%s -output_y_num %d", command, length(ofields.y))
  for (i in 1:length(ofields.y))
    command = sprintf("%s -f_y_field %s -f_y_from %d -f_y_to %d",
      command, ofields.y[i], ofield.y.from[i], ofield.y.to[i])

  system(paste("echo command: ", command, " >> ", log, sep=""))

  # log
  command = sprintf("%s >> %s 2>&1", command, log)

  cat(sprintf("%s\n", command))

  system(paste("echo hostname: `hostname` > ", log, sep=""))
  system(paste("echo command: ", command, " >> ", log, sep=""))

  # remove ofn before run
  system(paste("rm -rf", ofn))

  sleep = round(runif(1, 1, 3))
  system(paste("echo sleeping ",sleep," before command ... >> ", log, sep=""))
  system(paste("sleep", sleep))
  system(paste("echo sleep done, running command >> ", log, sep=""))

  ec = system(command)
  if (ec != 0)
  {
    cat(sprintf("command: %s\n", command))
    stop(paste("Command failed, error code:", ec))
  }

  return (0)
}

expand.ranges=function(fields=c("f1", "f2") , from.vals=c(1,1), to.vals=c(10,10), num.splits=4)
{
  dim = length(fields)
  num.splits.dim = num.splits^(1/dim)
  field.nparts = rep(1, length(fields))
  for (i in 1:dim)
    field.nparts[i] = min(to.vals[i] - from.vals[i] + 1, ceiling(num.splits.dim))

  params = list()
  for (i in 1:dim)
  {
    levels = seq(from.vals[i], to.vals[i])
    if (field.nparts[i] > 1) {
      lookup.table = data.frame(t(sapply(split(levels ,cut(levels, field.nparts[i])), range)))
      rownames(lookup.table) = NULL
      colnames(lookup.table) = c("from", "to")
    } else {
      lookup.table = data.frame(from=from.vals[i], to=to.vals[i])
    }
    params[[fields[i]]] = lookup.table
  }
  index.list = rev(lapply(params, function(x) 1:(dim(x)[1])))
  egrid = expand.grid(index.list, stringsAsFactors=F)

  result = NULL
  for (i in 1:dim)
  {
    ind = egrid[,fields[i]]
    lookup.table = params[[fields[i]]]
    t = lookup.table[ind,]
    colnames(t) = paste(fields[i], c("from", "to"), sep="_")
    rownames(t) = NULL
    result = if(i==1) t else cbind(result, t)
  }

  result
}

model.predict.split=function(
    tmp.dir = "/home/eitany/storage/o3c/mp_tmp/",
    log.dir = "/home/eitany/storage/o3c/mp_log/",
    wd="/home/eitany/storage/o3c/",
    binary = "/home/eitany/storage/o3c/bin/model_integrate",
    fends.fn="/home/eitany/storage/o3c/results/h1_cl_gc_1000.binned",
    ofn="/home/eitany/storage/o3c/results/h1_cl_gc_1000.expected_counts",
    num.splits=10, req.mem=2000, scope, model,
    max.njobs=400, max.jobs.fn="/home/eitany/maxjobs", prior=1, fends=NULL,
    mfields=NULL, mfields.maxvals=NULL,
    mfields.fns=NULL,
    dtype=NULL,
    ofields.x=c("frag_len_bin","frag_gc_bin"), from.vals.x=c(1,1), to.vals.x=c(5,5),
    ofields.y=c("frag_len_bin","frag_gc_bin"), from.vals.y=c(1,1), to.vals.y=c(5,5), Rcall)
{
  system(paste("rm -rf", tmp.dir, log.dir))
  system(paste("mkdir -p", tmp.dir))
  system(paste("mkdir -p", log.dir))

  cat(sprintf("binary: %s\n", binary))
  cat(sprintf("log dir: %s\n", log.dir))
  cat(sprintf("tmp dir: %s\n", tmp.dir))
  cat(sprintf("fends file: %s\n", fends.fn))
  cat(sprintf("ofn: %s\n", ofn))
  cat(sprintf("prior: %g\n", prior))
  cat(sprintf("scope: %s\n", scope))
  cat(sprintf("model: %s\n", model))
  if (!is.null(mfields))
    cat(sprintf("model fields: %s (%s)\n", paste(mfields, collapse=","), paste(mfields.maxvals, collapse=",")))
  else
    cat(sprintf("no model fields\n", paste(mfields, collapse=",")))
  cat(sprintf("output fields X: %s (%s)\n", paste(ofields.x, collapse=","),
              paste(from.vals.x, to.vals.x, sep="-", collapse=",")))
  cat(sprintf("output fields Y: %s (%s)\n", paste(ofields.y, collapse=","),
              paste(from.vals.y, to.vals.y, sep="-", collapse=",")))

  ofields = c(paste(ofields.x, "x", sep="_"), paste(ofields.y, "y", sep="_"))
  from.vals = c(from.vals.x, from.vals.y)
  to.vals = c(to.vals.x, to.vals.y)
  nbins = prod(to.vals - from.vals + 1)

  cat(sprintf("Splitting jobs according to output bin space\n"))
  ranges = expand.ranges(fields=ofields, from.vals=from.vals, to.vals=to.vals, num.splits=num.splits)
  ranges$from.fend = rep(0, dim(ranges)[1])
  ranges$to.fend = rep(0, dim(ranges)[1])
  concat = T

  # add running index field
  ranges$index = 1:dim(ranges)[1]

  # build commands
  commands = c()
  ofns = c()

  vec2str=function(x) {
    paste(x, collapse=" ")
  }

  Rcall = paste(Rcall, collapse=" ")

  for (i in 1:nrow(ranges)) {
    r = ranges[i,]

    ofields.x.ranges = get.ofield.range(r, ofields.x, "x")
    ofields.y.ranges = get.ofield.range(r, ofields.y, "y")

    command = sprintf("%s %s R/model_predict.r model.predict ofields.x=\"%s\" ofield.x.from=\"%s\" ofield.x.to=\"%s\" ofields.y=\"%s\" ofield.y.from=\"%s\" ofield.y.to=\"%s\" index=%d from.fend=%d to.fend=%d fends.fn=\"%s\" log.dir=\"%s\" tmp.dir=\"%s\" binary=\"%s\" prior=%g scope=\"%s\" model=\"%s\"",
        Rcall, wd, vec2str(ofields.x), vec2str(ofields.x.ranges$from), vec2str(ofields.x.ranges$to),
        vec2str(ofields.y), vec2str(ofields.y.ranges$from), vec2str(ofields.y.ranges$to),
        r$index, r$from.fend, r$to.fend, fends.fn, log.dir, tmp.dir,
        binary, prior, scope, model)

    if (!is.null(mfields)) {
      command = sprintf("%s mfields=\"%s\" mfields.maxvals=\"%s\" mfields.fns=\"%s\"", command, vec2str(mfields), vec2str(mfields.maxvals), vec2str(mfields.fns))
    }

    commands = c(commands, command)

    ofns = c(ofns, paste(tmp.dir, "/", get.short.fn(fends.fn), ".", r$index, sep=""))
  }

  ntasks = dim(ranges)[1]
  njobs = min(ntasks, max.njobs)
  cat(sprintf("Number of commands: %d\n", length(commands)))
  cat(sprintf("Number of bins: %g\n", nbins))
  cat(sprintf("Number of tasks: %d\n", ntasks))
  cat(sprintf("Average square bins per tasks: %d\n", round(nbins/ntasks)))
  cat(sprintf("Number of running jobs: %d\n", njobs))
  system(paste("rm -rf ", tmp.dir, get.short.fn(fends.fn), "*", sep=""))

  system(paste("rm -rf ", log.dir, get.short.fn(fends.fn), "*", sep=""))
  system(paste("rm -rf ", tmp.dir, get.short.fn(fends.fn), "*", sep=""))
  system(paste("rm -rf ", ofn, sep=""))

  source(paste(wd, "/R/distrib.r", sep=""))
  time = proc.time()
  distrib(commands=commands, qsub.dir=log.dir, jobname="econtact", batch.max.jobs=njobs, total.max.jobs.fn=max.jobs.fn, sleep.time=1, req.mem=req.mem,
	  dtype=dtype)

  cat(sprintf("distributing expected counts: %f minutes\n", (proc.time() - time)[3]/60))

  cat(sprintf("Merging results into %s\n", ofn))
  if (concat) {
    for (i in seq_along(ofns)) {

      rep = 10
      while (!file.exists(ofns[i])) {
        Sys.sleep(60)
        rep = rep - 1
        if (rep == 0)
          stop(paste("file does not exist:", ofns[i]))
      }

      if (i==1)
        ec = system(paste("cat ", ofns[i], " > ", ofn, sep=""))
      else
        ec = system(paste("tail -n +2 ", ofns[i], " >> ", ofn, sep=""))
      if (ec != 0)
        stop("Error merging results")
    }
  } else {
    ofields = c(paste(ofields.x, "1", sep=""), paste(ofields.y, "2", sep=""))
    table = NULL
    for (i in seq_along(ofns)) {
      if (!file.exists(ofns[i]))
        stop(paste("file does not exist:", ofns[i]))

      ttable = read.delim(ofns[i], stringsAsFactors=F)
      if (i == 1) {
        table = ttable
        next
      }
      if (any(ttable[,ofields] != table[,ofields]))
        stop("result does not have the same field order")
      table$value = table$value + ttable$value
    }
    write.table(table, ofn, quote=F, col.names=T, row.names=F, sep="\t")
  }
}

omit.dups=function(table, nfields)
{
  diff = table[,1:nfields+nfields] - table[,1:nfields]

  if (nfields == 1)
    return (table[diff>=0,])

  ind = rep(F, dim(table)[1])
  eq = rep(T, dim(table)[1])
  for (i in 1:nfields)
  {
    if (i<nfields)
      ind = ind | (eq & (diff[,i] > 0))
    else
      ind = ind | (eq & (diff[,i] >= 0 ))
    eq = eq & (diff[,i] == 0)
  }
  table[ind,]
}

#########################################################################################################
# Wrapper functions
#########################################################################################################

# used to setup the model - get total possible counts for each bin
compute.total.counts=function(tmp.dir = "/home/eitany/storage/o3c/mp_tmp/",
                              log.dir = "/home/eitany/storage/o3c/mp_log/",
                              wd = "/home/eitany/storage/o3c/",
                              binary = "/home/eitany/storage/o3c/bin/model_integrate",
                              prefix="/home/eitany/storage/o3c/results/h1_cl_gc_200",
                              num.splits=400, max.jobs.fn="/home/eitany/maxjobs", fends.suffix="binned",
                              scope, req.mem=2000, dtype="sge", model.ifn, Rcall)
{
  fends.fn = paste(prefix, ".", fends.suffix, sep="")
  table = load.table(fends.fn)
  fends = table$fend
  model = "count"

  mtable = load.table(model.ifn)
  mfields = mtable$field
  maxvals = mtable$size
  cat(sprintf("mfn fields: %s\n", paste(mfields, collapse=",")))
  ofields = mfields

  ofn = paste(prefix, ".total_counts", sep="")
  model.predict.split(tmp.dir=tmp.dir, log.dir=log.dir, wd=wd, binary=binary,
                      fends.fn=fends.fn, ofn=ofn, fends=fends, num.splits=num.splits,
                      scope=scope, model=model, max.jobs.fn=max.jobs.fn,
                      ofields.x=ofields, from.vals.x=rep(1, length(ofields)), to.vals.x=maxvals,
                      ofields.y=ofields, from.vals.y=rep(1, length(ofields)), to.vals.y=maxvals,
                      req.mem=req.mem, dtype=dtype, Rcall=Rcall)

  table = read.delim(ofn, stringsAsFactors=F)
  result = omit.dups(table, length(ofields))
  names(result)[dim(result)[2]] = "total"
  write.table(result, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}

# used to check model performance in a user-defined xy plane
compute.expected.counts=function(
    tmp.dir = "/home/eitany/storage/o3c/mp_tmp/",
    log.dir = "/home/eitany/storage/o3c/mp_log/",
    wd = "/home/eitany/storage/o3c/",
    binary = "/home/eitany/storage/o3c/bin/model_integrate",
    model.prefix="/home/eitany/storage/o3c/results/h1_cl_gc_1000_1M",
    fends.ifn, num.splits=100, max.jobs.fn="/home/eitany/maxjobs",
    scope, model,
    model.ifn="/home/eitany/storage/o3c/results/map_len_gc.model",
    ofn="/home/eitany/storage/o3c/results/h1_cl_gc_1000_1M.e_contact",
    ofields.x=c("coord_bin"), ofields.y=c("coord_bin"), req.mem=2000, dtype,
    omit.x.zero=F, omit.y.zero=F, Rcall)
{
  prior.ifn = paste(model.prefix, "prior", sep=".")
  cat(sprintf("reading prior file %s\n", prior.ifn))
  prior = read.delim(prior.ifn, header=F)
  prior = prior[1,1]

  # get all bin ranges
  cat(sprintf("reading fends table %s\n", fends.ifn))
  table = read.delim(fends.ifn, stringsAsFactors=F)
  fends = table$fend
  for (field in c(ofields.x, ofields.y))
    if (!is.element(field, names(table)))
      stop(paste("field", field, "not found"))
  from.vals.x = numeric(length(ofields.x))
  to.vals.x = numeric(length(ofields.x))
  from.vals.y = numeric(length(ofields.y))
  to.vals.y = numeric(length(ofields.y))
  for (i in seq_along(ofields.x))
  {
    vx = unique(table[,ofields.x[i]])
    if (omit.x.zero)
      vx = setdiff(vx, 0)
    rx = range(vx)
    from.vals.x[i] = rx[1]
    to.vals.x[i] = rx[2]
  }
  for (i in seq_along(ofields.y))
  {
    vy = unique(table[,ofields.y[i]])
    if (omit.y.zero)
      vy = setdiff(vy, 0)
    ry = range(vy)
    from.vals.y[i] = ry[1]
    to.vals.y[i] = ry[2]
  }

  # model parameters
  cat(sprintf("reading mtable file %s\n", model.ifn))
  mtable = read.delim(model.ifn)
  if (dim(mtable)[1] > 0)
  {
    mfields = mtable$field
    mfields.maxvals = mtable$size
    mfields.fns = paste(paste(model.prefix, mtable$field, sep="_"), ".f", sep="")
  }
  else {
    mfields = NULL
    mfields.maxvals = NULL
    mfields.fns = NULL
  }
  model.predict.split(tmp.dir=tmp.dir, log.dir=log.dir, wd=wd, binary=binary,
                      fends.fn=fends.ifn, ofn=ofn, fends=fends, num.splits=num.splits,
                      scope=scope, model=model, prior=prior,
                      max.jobs.fn=max.jobs.fn,
                      mfields=mfields, mfields.maxvals=mfields.maxvals,
                      mfields.fns=mfields.fns,
                      ofields.x=ofields.x, from.vals.x=from.vals.x, to.vals.x=to.vals.x,
                      ofields.y=ofields.y, from.vals.y=from.vals.y, to.vals.y=to.vals.y, req.mem=req.mem, dtype=dtype, Rcall=Rcall)
}
