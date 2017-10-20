#####################################################################################
# Background model of an open ended experiment
#
# We optimize the following function:
#
# P(Xij) = Pr * F^L(xi, xj) * F^E(yi, yj)
#
# Xij := RV representing fragends i,j being connected
# Pr  := Prior on pairs beings connected
# F^L(xi, xj) := connection bias of bins xi,xj (binned according to fragment length)
# F^E(yi, yj) := connection bias of bins yi,yj (binned according to fragend length)
#
# In our implementation we alternate between holding F^L and optimizing F^E
# and vice versa.
#
#####################################################################################

#####################################################################################
# Utility functions
#####################################################################################

format.number=function(num)
{
  size = c(10^9, 10^6, 10^3)
  name = c("G", "M", "K")
  for (i in 1:length(size))
    if (num > size[i])
      return (paste(round(num/size[i]), name[i], sep=""))
  return (num)
}

counts2labels=function(F.seed)
{
  N = dim(F.seed)[1]
  result = vector("character", N)
  for (i in 1:N)
  {
    count = F.seed[i,"count"]
    total = F.seed[i,"total"]
    count = format.number(count)
    total = format.number(total)
    result[i] = paste(count, "/", total, sep="")
  }
  result
}


counts2probs=function(F.seed, Pr)
{
  N = dim(F.seed)[1]
  probs = vector("numeric", N)
  for (i in 1:N)
    probs[i] = F.seed[i,"count"] / F.seed[i,"total"]
  probs[!is.finite(probs)] = 0
  probs = probs / Pr
  result = cbind(F.seed[,1:2], probs)
  result
}

get.seeds=function(count.table, field="frag_len", is.const)
{
  cat(sprintf("computing seed function for %s\n", field))
  bfields = paste(field, c(1, 2), sep="")
  keys = count.table[, bfields]
  keys = unique(keys)
  ind = keys[,1] <= keys[,2]
  keys = keys[ind,]
  N = dim(keys)[1]
  count = vector("integer", N)
  total = vector("integer", N)

  for (i in 1:N)
  {
    bin1 = keys[i,1]
    bin2 = keys[i,2]
    ind = (count.table[,bfields[1]] == bin1 & count.table[,bfields[2]] == bin2) |
          (count.table[,bfields[1]] == bin2 & count.table[,bfields[2]] == bin1)
    tmp = count.table[ind,]
    count[i] = sum(as.numeric(tmp$count)) + ifelse(is.const,0,1)
    total[i] = sum(as.numeric(tmp$total)) + ifelse(is.const,0,1)
  }
  result = cbind(keys, count, total)
  rownames(result)=NULL
  result
}

extend.func=function(func)
{
  ind = func[,1] != func[,2]
  tmp = func[ind,c(2,1,3)]
  names(tmp) = names(func)
  rbind(func, tmp)
}

append.func.to.counts=function(count.table, fields, functions)
{
  for (i in length(fields):1)
  {
    field = fields[i]
    func = extend.func(functions[[i]])
    bfields = paste(field, c(1, 2), sep="")

    keys = paste(count.table[,bfields[1]], count.table[,bfields[2]], sep="_")
    lu.keys = paste(func[,1], func[,2], sep="_")
    lu.values = func[,3]
    indices = as.numeric(factor(keys, levels=lu.keys))
    values = lu.values[indices]
    count.table[[field]] = values
  }

  fields1 = paste(fields, 1, sep="")
  fields2 = paste(fields, 2, sep="")
  count.table = count.table[,c(fields1, fields2, "count", "total", fields)]
  count.table
}

get.func.value=function(func, key1, key2)
{
  ind = func[,1] == key1 & func[,2] == key2
  if (sum(ind) != 1)
    stop("more than one value")
  func[ind, "probs"]
}

#####################################################################################
# Algorithm functions
#####################################################################################

# use ML to compute seed functions
compute.ll=function(count.table, fields, functions, Pr)
{
  cat(sprintf("appending function values to table\n"))
  table = append.func.to.counts(count.table, fields, functions)

  # omit zero probability events
  cat(sprintf("omit zero prob\n"))
  table = table[table$total!=0,]

  N = dim(table)[1]

  cat(sprintf("compute sum log and prod\n"))

  # log sum vector
  sum.log = rep(0.0, N)
  for (i in 1:length(fields))
    sum.log = sum.log + log2(table[,fields[i]])

  # product vector
  prod = rep(1.0, N)
  for (i in 1:length(fields))
    prod = prod * table[,fields[i]]

  n = table[,"count"]
  m = table[,"total"] - n
  result = n*(log2(Pr) + sum.log) + m*log2(1 - Pr*prod)

  result = sum(result)
  if (is.nan(result)) {
    browser()
    stop("reduce number of model parameters")
  }

  result
}

# single step in alternating algorithm:
#  maximize one function, hold other const
maximize.ll=function(count.table, fields, functions, Pr, func.i)
{
  table = append.func.to.counts(count.table, fields, functions)

  # omit zero probability events
  table = table[table$total!=0,]

  # indices w/o function we are optimizing
  const.ind = setdiff(1:length(fields), func.i)

  # optimize each key separately
  bfields = paste(fields[func.i], 1:2, sep="")
  func = functions[[func.i]]
  keys = func[,1:2]

  print(paste("optimizing", dim(keys)[1], "values for field", fields[func.i]))
  probs = vector("numeric", dim(keys)[1])
  for (k in 1:dim(keys)[1])
  {
    key1 = keys[k,1]
    key2 = keys[k,2]
    prev.value = get.func.value(functions[[func.i]], key1, key2)

    ind = (table[,bfields[1]] == key1 & table[,bfields[2]] == key2) |
          (table[,bfields[2]] == key1 & table[,bfields[1]] == key2)
    table.key = table[ind,]
    N = dim(table.key)[1]

    if (N == 0)
    {
      probs[k] = prev.value
      next
    }

    # log sum vector
    sum.log = rep(0.0, N)
    for (i in const.ind)
      sum.log = sum.log + log2(table.key[,fields[i]])

    # product vector
    prod = rep(1.0, N)
    for (i in const.ind)
      prod = prod * table.key[,fields[i]]

    n = table.key[,"count"]
    m = table.key[,"total"] - n

    ll = function(alpha) {
      sapply(alpha, function(alpha) sum(n*(log2(Pr) + sum.log + log2(alpha)) + m*log2(1 - Pr*prod*alpha)))
    }
    ll.grad = function(alpha) {
      sapply(alpha, function(alpha) sum(n/(log(2)*alpha) - (m*Pr*prod)/(log(2)*(1-Pr*prod*alpha))))
    }

    lower.bound = 0
    upper.bound = min(1/(Pr*prod))
    start = prev.value
    if (start > upper.bound) start = upper.bound
    if (start < lower.bound) start = lower.bound

    A = matrix(c(1, -1), 2, 1)
    B = matrix(c(-lower.bound, upper.bound), 2, 1)

    # use analytic gradient
    mnr = maxBFGS(ll, start=start, grad=ll.grad, constraints=list(ineqA=A, ineqB=B))
    if (mnr$code != 0)
      stop(mnr$message)

    probs[k] = mnr$estimate
  }

  result = cbind(keys, probs)
  return (result)
}

#####################################################################################
# Wrapper functions
#####################################################################################

learn.model.file=function(model.prefix, model.params, fields)
{
  ifn.count = paste(model.prefix, ".nm", sep="")
  count.table = read.delim(ifn.count, stringsAsFactors=F)

  if (sum(count.table$count) < 10)
    stop("<10 contacts found, cannot learn model")

  # prior on interaction
  Pr = sum(as.numeric(count.table$count)) / sum(as.numeric(count.table$total))

  pfn = paste(model.prefix, ".prior", sep="")
  write.table(list(prior=Pr), pfn,
                   quote=F, col.names=F, row.names=F, sep="\t")
  cat(sprintf("generating prior file %s\n", pfn))

  F.seed = list()
  functions = list()
  seed.labels = list()
  for (field in fields)
  {
      is.const = (model.params[[field]]$type == "const")
      F.seed[[field]] = get.seeds(count.table, field, is.const=is.const)
      functions[[field]] = counts2probs(F.seed[[field]], Pr)
      seed.labels[[field]] = counts2labels(F.seed[[field]])
  }

  # save functions to file
  for (i in 1:length(fields))
  {
    field = fields[i]
    func = functions[[i]]
    func = func[order(func[,1], func[,2]),]
    func[,3] = round(func[,3],8)
    write.table(func, paste(model.prefix, "_", fields[i], "_seed.f", sep=""),
                quote=F, col.names=T, row.names=F, sep="\t")
  }

  # replace seed function for const functions (e.g. mappability bias)
  for (field in fields)
  {
    if (model.params[[field]]$type == "const")
    {
      seed.fn = paste(model.prefix, "_", field, ".f", sep="")
      cat(sprintf("loading const function for field: %s\n", field))
      functions[[field]] = read.delim(seed.fn)
      seed.labels[[field]] = NULL
    }
  }

  ll = compute.ll(count.table=count.table, fields=fields, functions=functions, Pr=Pr)
  print(paste("Initial LL=", ll, sep=""))
  system(paste("echo ", ll, " > ", model.prefix, ".initial_ll", sep=""))

  max.iter = 10
  if (length(fields) > 1)
  {
    for (i in 1:max.iter)
    {
      print(paste(">>> Iteration #", i, sep=""))
      for (j in 1:length(fields))
      {
        field = fields[j]
        if (model.params[[field]]$type == "const")
          next

        result = maximize.ll(count.table=count.table, fields=fields, functions=functions,
                             Pr=Pr, func.i=j)
        functions[[j]] = result
        ll.new = compute.ll(count.table, fields=fields, functions=functions, Pr)
        delta = round(ll.new-ll,2)
        print(paste("LL delta=", delta, sep=""))
        ll = ll.new
      }
      if (delta < 1)
      {
        system(paste("echo ", ll, " > ", model.prefix, ".final_ll", sep=""))
        break
      }
    }
  }

  # save functions to file
  for (i in 1:length(fields))
  {
    field = fields[i]
    if (model.params[[field]]$type == "const")
      next
    func = functions[[i]]
    func = func[order(func[,1], func[,2]),]
    func[,3] = round(func[,3],8)
    write.table(func, paste(model.prefix, "_", fields[i], ".f", sep=""),
                quote=F, col.names=T, row.names=F, sep="\t")
  }
}

learn.model=function(model.prefix="h_results/h1_cl_gc_200",
                     model.fn="h_results/map_len_gc.model")
{
    library(maxLik)

    model.table = load.table(model.fn)
    model.params = list()
    fields = NULL
    for (i in 1:dim(model.table)[1]) {
        field = model.table[i,"field"]
        type = model.table[i,"type"]
        if (type != "const" && type != "optimize")
            stop("in model table type must be 'optimize' or ' const'")
        model.params[[ field ]] = list(size=model.table[i,"size"], type=type)
        fields = c(fields, field)
    }

    learn.model.file(model.prefix=model.prefix, model.params=model.params, fields=fields)
}

rl=function() {
    source("md/anchors/R/model.r")
}
