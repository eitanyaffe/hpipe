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

counts2factors2=function(table, Pr)
{
  N = dim(table)[1]
  probs = vector("numeric", N)
  for (i in 1:N)
    probs[i] = table[i,"count"] / table[i,"total"]
  probs[!is.finite(probs)] = 0
  result = cbind(table[,1:2], probs)
  result$probs = result$probs/Pr
  result
}

get.seeds2=function(count.table, field="frag_len", is.const)
{
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
    result = count.table
    for (i in 1:length(fields)) {
        field = fields[i]
        func = as.data.frame(functions[[i]])

        if (dim(func)[2] == 2) {
            func$key = func[,1]
            keys = result[,field]
        } else {
            func = extend.func(func)
            func$key = paste(func[,1], func[,2], sep="_")
            bfields = paste(field, c(1, 2), sep="")
            keys = paste(result[,bfields[1]], result[,bfields[2]], sep="_")
        }

        result[,paste(field,"prob",sep="_")] = func$probs[match(keys,func$key)]
    }
    result
}

#####################################################################################
# Algorithm functions
#####################################################################################

compute.prior=function(count.table=count.table, fields=fields, functions=functions, min.probability)
{
    table = append.func.to.counts(count.table, fields, functions)
    N = dim(table)[1]
    pfields = paste(fields, "prob", sep="_")
    table$product = rep(1.0, N)
    for (i in 1:length(pfields))
        table$product = table$product * table[,pfields[i]]

    if (any(table$total == 0))
        stop(sprintf("zero total found\n"))

    bottom = max(min.probability/table$product)
    top = min((1-min.probability)/table$product)
    optimal = sum(table$count) / sum(table$prod)
    middle = (top + bottom) / 2

    if (optimal > bottom && optimal < top) {
        cat(sprintf("constant valid range=(%g,%g), middle optimal=%g\n", bottom, top, middle))
        result = optimal
    } else {
        cat(sprintf("constant valid range=(%g,%g), optimal prior=%g\n", bottom, top, optimal))
        result = middle
    }
    result
}

compute.ll=function(count.table, fields, functions, Pr)
{
  table = append.func.to.counts(count.table, fields, functions)

  pfields = paste(fields, "prob", sep="_")

  if (any(table$total == 0))
      stop(sprintf("zero total found\n"))

  N = dim(table)[1]

  # product vector
  prod = rep(1.0, N)
  for (i in 1:length(pfields))
    prod = prod * table[,pfields[i]]
  rr = range(Pr*prod)
  # cat(sprintf("ll prob range=(%g,%g)\n", rr[1], rr[2]))

  n = table[,"count"]
  m = table[,"total"] - n
  result = n*log2(Pr*prod) + m*log2(1 - Pr*prod)

  result = sum(result)
  if (is.nan(result)) {
    browser()
    stop("LL failed, try to reduce number of model parameters")
  }

  result
}

# single step in alternating algorithm:
#  maximize one function, hold other const
maximize.ll=function(count.table, fields, functions, Pr, func.i, min.probability)
{
  table = append.func.to.counts(count.table, fields, functions)
  pfields = paste(fields, "prob", sep="_")

  if (any(table$total == 0))
      stop(sprintf("zero total found\n"))

  # indices w/o function we are optimizing
  const.ind = setdiff(1:length(fields), func.i)

  func = functions[[func.i]]

  bfields = paste(fields[func.i], 1:2, sep="")
  keys = paste(func[,1], func[,2], sep="_")
  table$key = paste(table[,bfields[1]], table[,bfields[2]], sep="_")
  func$key = paste(func[,bfields[1]], func[,bfields[2]], sep="_")

  M = length(keys)
  cat(paste("optimizing", M, "values for field", fields[func.i]), "\n")

  result = func
  prev.values = result$probs[match(keys,result$key)]

  for (k in 1:M) {
      key = keys[k]
      prev.value = prev.values[k]
      table.key = table[key == table$key,]
      N = dim(table.key)[1]

      if (N == 0) {
          result$probs[k] = prev.value
          next
      }

      # product vector
      prod = rep(1.0, N)
      for (i in const.ind)
          prod = prod * table.key[,pfields[i]]

      n = table.key[,"count"]
      m = table.key[,"total"] - n

      # rr = range(Pr*prod*prev.value)
      # cat(sprintf("range(Pr*prod*prev.value)=(%g,%g)\n", rr[1], rr[2]))

      ll = function(alpha) {
          sapply(alpha, function(alpha) sum(n*(log2(Pr*prod*alpha)) + m*log2(1 - Pr*prod*alpha)))
      }
      ll.grad = function(alpha) {
          sapply(alpha, function(alpha) sum(n/(log(2)*alpha) - (m*Pr*prod)/(log(2)*(1-Pr*prod*alpha))))
      }

      lower.bound = max(min.probability/Pr*prod)
      upper.bound = min((1-min.probability)/Pr*prod)
      start = prev.value
      if (start > upper.bound) start = upper.bound
      if (start < lower.bound) start = lower.bound

      A = matrix(c(1, -1), 2, 1)
      B = matrix(c(-lower.bound, upper.bound), 2, 1)

      # use analytic gradient
      mnr = maxBFGS(ll, start=start, grad=ll.grad, constraints=list(ineqA=A, ineqB=B))
#      mnr = maxBFGS(ll, start=start, grad=NULL, constraints=list(ineqA=A, ineqB=B))
      if (mnr$code != 0) {
          browser()
          stop(mnr$message)
      }

      result$probs[k] = mnr$estimate
  }

  return (result[,-match("key",names(result))])
}

#####################################################################################
# Wrapper functions
#####################################################################################

learn.model.file=function(model.prefix, model.params, fields)
{
    count.table = load.table(paste(model.prefix, ".nm", sep=""))

    # prior
    Pr = sum(as.numeric(count.table$count)) / sum(as.numeric(count.table$total))
    pfn = paste(model.prefix, ".prior", sep="")
    cat(sprintf("generating prior file %s\n", pfn))
    write.table(list(prior=Pr), pfn, quote=F, col.names=F, row.names=F, sep="\t")

    cat("computing seed functions\n");
    seed.counts = list()
    seed.functions = list()
    for (field in fields) {
        type = model.params[[field]]$type
        is.const = (type == "const")
        seed.counts[[field]] = get.seeds2(count.table, field, is.const=is.const)
        seed.functions[[field]] = counts2factors2(seed.counts[[field]], Pr)
    }

    # save functions to file
    for (i in 1:length(fields)) {
        field = fields[i]
        func = seed.functions[[i]]
        func = func[order(func[,1], func[,2]),]
        func[,3] = round(func[,3],8)
        save.table(func, paste(model.prefix, "_", fields[i], "_seed.f", sep=""))
    }

    # replace seed function for const functions (e.g. mappability bias)
    for (field in fields) {
        if (model.params[[field]]$type == "const") {
            seed.fn = paste(model.prefix, "_", field, ".f", sep="")
            seed.functions[[field]] = load.table(seed.fn)
        }
    }

    # print out function ranges
    for (field in fields) {
        type = model.params[[field]]$type
        rr = range(seed.functions[[field]]$probs)
        cat(sprintf("seed field=%s, type=%s, range: (%f,%f)\n", field, type, rr[1], rr[2]))
    }

    min.probability = 10^-100
    min.probability = 0

    # Pr = compute.prior(count.table=count.table, fields=fields, functions=seed.functions, min.probability=min.probability)
    # pfn = paste(model.prefix, ".prior", sep="")
    # write.table(list(prior=Pr), pfn, quote=F, col.names=F, row.names=F, sep="\t")

    ll = compute.ll(count.table=count.table, fields=fields, functions=seed.functions, Pr=Pr)
    cat(paste("Initial LL=", ll, sep=""), "\n")
    system(paste("echo ", ll, " > ", model.prefix, ".initial_ll", sep=""))

    functions = seed.functions

    max.iter = 10
    for (iter.index in 1:max.iter) {
        cat(paste(">>> Iteration #", iter.index, sep=""), "\n")
        for (func.i in 1:length(fields)) {
            field = fields[func.i]
            if (model.params[[field]]$type == "const")
                next

            result = maximize.ll(
                count.table=count.table, fields=fields, functions=functions,
                Pr=Pr, func.i=func.i, min.probability=min.probability)

            cat(sprintf("field: %s, sqrt(mean(delta^2))=%g\n", field, sqrt(mean((functions[[func.i]]$probs - result)^2))))
            functions[[func.i]] = result
            ll.new = compute.ll(count.table=count.table, fields=fields, functions=functions, Pr)
            delta = round(ll.new-ll,2)
            cat(paste("LL delta=", delta, sep=""), "\n")
            ll = ll.new
        }
        if (delta < 1) {
            system(paste("echo ", ll, " > ", model.prefix, ".final_ll", sep=""))
            break
        }
    }

    # print out function ranges
    for (field in fields) {
        type = model.params[[field]]$type
        rr = range(seed.functions[[field]]$probs)
        cat(sprintf("final field=%s, type=%s, range: (%f,%f)\n", field, type, rr[1], rr[2]))
    }

    # save functions to file
    for (i in 1:length(fields)) {
        field = fields[i]
        if (model.params[[field]]$type == "const")
            next
        func = functions[[i]]

        func = func[order(func[,1], func[,2]),]
        func[,3] = round(func[,3],8)
        ofn = paste(model.prefix, "_", fields[i], ".f", sep="")
        save.table(func, ofn)
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
        model.params[[ field ]] = list(size=model.table[i,"size"], type=type, dim=2)
        fields = c(fields, field)
    }

    learn.model.file(model.prefix=model.prefix, model.params=model.params, fields=fields)
}

rl=function() {
    source("md/anchors/R/model.r")
}
