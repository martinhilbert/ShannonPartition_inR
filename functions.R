############################################################################################################
# SHANNON INFORMATION PARTITIONS: by Martin Hilbert & Mahima Agarwal @UCDavis.edu (unreleased version 0.0) #
############################################################################################################
#  It was adopted from an earlier Python implementation by Ryan James: https://github.com/Autoplectic/dit 
############################################################################################################
# This file contains helper functions required for creating the distributions 
# and information partitions
###############################################################################

# 1. Construct alphabets from outcomes for distribution

construct.alphabets<- function(x, sort=TRUE)
{
  l = length(x)
  outcome.class = class(x[0])
  if (l == 0)
  {
    stop('Outcomes must not be empty.')
  }
  if (class(x) == 'character')
  {
    x =strsplit(x, '')
  }
  # Confirm that all outcomes have the same length
  if (length(unique(unlist(lapply(x, length)))) != 1)
  {
    stop('All outcomes must have equal length')
  }
  if (outcome.class == 'character')
  {
    x = matrix(unlist(x), ncol = length(x[[1]]), byrow = TRUE)
  }
  outcome.length = ncol(x)
  alphabets = apply(x, 2, unique)
  if (class(alphabets) != 'list')
  {
    alphabets<- as.list(as.data.frame(alphabets, stringsAsFactors = F))
  }
  if (sort)
  {
    alphabets = lapply(alphabets, sort)
  }
  return(list(alphabet=alphabets, outcome_length=outcome.length, outcome_class=outcome.class))
}

###############################################################################

# 2. Calculate shannon entropy for sets of random variables for a distribution

Shannon.Entropy<- function(distribution, rvs=NULL, rv.mode=NULL)
{
  if (distribution$is_joint)
  {
    if(is.null(rvs))
    {
      rvs <- seq(distribution$outcome_length)
      rv.mode <- 1
    }
    d = distribution.marginal(distribution, rvs, rv.mode)
    
  }else
  {
    d = distribution
  }
  pmf = d$pmf
  if (d$is_log)
  {
    base = d$base
    terms = -base^pmf * pmf
  }else
  {
    terms = -pmf * log(pmf, 2)
  }
  H = sum(terms, na.rm = TRUE)
  return(H)
}

###############################################################################

# 3. Parse random variables, and return corresponding indexes in the distribution 

Parse.RVs<- function(rvs, distribution, rv.mode=NULL, unique=TRUE, sort=TRUE)
{
  if (is.null(rv.mode))
  {
    rv.mode <- distribution$rv_mode
  }
  rv.mode<- RV.MODES[[rv.mode]]
  if (length(rvs) == 0)
  {
    return(data.frame(rvs=numeric(0), indexes=numeric(0)))
  }
  
  #Check that all random variables are unique
  if (sum(duplicated(rvs)) > 0)
  {
    error.msg = 'Duplicate values not allowed in rvs'
    stop(error.msg)
  }
  
  if (rv.mode == RV.MODES$NAMES)
  {
    # Convert random variable names to indices
    if (is.null(distribution$rvs))
    {
      stop('There are no random variable names to use.')
    }
    indexes = which(distribution$rvs %in% rvs)
    if (length(indexes) != length(rvs))
    {
      error.msg = 'rvs contain invalid random variable names'
      stop(error.msg)
    }
  }else
  {
    indexes=rvs
  }
  
  # Check validity of indexes
  all_indexes = unique(seq(distribution$outcome_length))
  if (length(which(!(indexes %in% all_indexes))) > 0)
  {
    error.msg = paste('rvs contain invalid random variables:', 
                      paste(indexes[which(!(indexes %in% all_indexes))], collapse=', '))
    stop(error.msg)
  }
  
  # Sort random variables by their index
  out<- data.frame(rvs=rvs, indexes=indexes)
  if (sort)
  {
    out <- out[order(out[,1]),]
  }
  return(out)
}

###############################################################################

# 4. Generate sample space for the distribution from alphabets

generate.sample.space<- function(alphabet, ctor, indexes=NULL)
{
  if (!(is_null(indexes)))
  {
    alphabet = alphabet[indexes]
  }
  sample_space = sapply(cross(alphabet), ctor)
  return(sample_space)
}

###############################################################################
# 5. Create a power set of a set

power.set<- function(x)
{
  N<- length(x)
  if (N == 1)
  {
    return(list(c(), x))
  }
  null<- list(c())
  l<- vector('list', N)
  for (i in (1:N))
  {
    l[[i]] = combn(x,i, simplify=FALSE)
  }
  l = append(null, unlist(l, recursive=FALSE))
  return(l)
}

###############################################################################

# 6. Reverse a graph

reverse<- function(g)
{
  el<- as_edgelist(g)
  reversed.graph<- graph(rbind(el[,2],el[,1]))
  return(reversed.graph)
}

###############################################################################

# 7. Convert names of vertices to random variable names

convert.vertex.name<- function(x, len=FALSE, str=FALSE)
{
  x<- unlist(strsplit(x, ','))
  if (length(x) == 1)
  {
    if (x == '()')
    {
      if (len)
      {
        return(0)
      }
      if (str)
      {
        return(character(0))
      }
      return(numeric(0))
    }
  }
  res<-as.character(x)
  if (!str)
  {
    res<-as.numeric(res)
  }
  if (len)
  {
    return(length(res))
  }
  return(res)
}

###############################################################################

# 8. Get appropriate symbol (H/I)

symbol<- function(rvs, crvs)
{
  if (length(rvs) == 1)
  {
    return('H')
  }
  return('I')
}

###############################################################################

# 9. Return elements of a tuple

get.tuple.nos<- function(x)
{
  if (class(x) == 'data.frame')
  {
    x<- unlist(x)
  }
  x<- trimws(as.character(x))
  parts<- strsplit(x, '[(,)]',perl=TRUE)[[1]][-1]
  return(as.numeric(parts))
}
###############################################################################