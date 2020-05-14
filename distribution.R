############################################################################################################
# SHANNON INFORMATION PARTITIONS: by Martin Hilbert & Mahima Agarwal @UCDavis.edu (unreleased version 0.0) #
############################################################################################################
#  It was adopted from an earlier Python implementation by Ryan James: https://github.com/Autoplectic/dit 
############################################################################################################
# This file contains functions around distributions.
###############################################################################

# 1. Create distribution in the form of a list

# Required input: Matrix or dataframe
# pmf: If pmf is not supplied, last column f outcomes is assumed to describe the pmf
create.distribution<- function(outcomes, pmf=NULL, rvs=NULL, sample.space=NULL, sort=TRUE)
{
  dist = list()
  dist$is_joint = TRUE
  if (is.null(pmf))
  {
    p<- as.numeric(outcomes[,ncol(outcomes)])
    outcomes<- outcomes[,-ncol(outcomes)]
  }else
  {
    p<- as.numeric(pmf)
  }
  alphabets = construct.alphabets(outcomes, sort)
  dist$alphabet = alphabets$alphabet
  dist$outcome_class = alphabets$outcome_class
  dist$outcome_length = alphabets$outcome_length
  if (dist$outcome_class == 'character')
  {
    dist$x = data.frame(matrix(unlist(strsplit(outcomes, '')), ncol = dist$outcome_length, byrow=TRUE))
    dist$outcome_ctor<- function(x){paste(x, collapse='')}
  }else
  {
    dist$x = data.frame(outcomes, stringsAsFactors = FALSE)
    dist$outcome_ctor<- function(x){paste('(', paste(x, collapse=', '), ')', sep='')}
  }
  dist$p_table = data.frame(x=outcomes, p_x=p)
  dist$rv_mode = 1
  dist$pmf = p
  dist$rvs = rvs
  dist$base = 'linear'
  if (is.null(sample.space))
  {
    sample.space = generate.sample.space(dist$alphabet, dist$outcome_ctor)
  }
  dist$sample_space = sample.space
  dist$mask = rep(FALSE, dist$outcome_length)
  dist$is_log=FALSE
  return(dist)
}

###############################################################################
# 2. Create distribution for subsets of random variables

# This function does not handle more than one set of random variables at a time
distribution.coalesce<- function(distribution, rvs, rv.mode=NULL, extract=FALSE)
{
  parse<- function(x){Parse.RVs(x, distribution, rv.mode)[,2]}
  indexes<- lapply(rvs, parse)
  if (length(rvs) == 1)
  {
    if (extract)
    {
      ctor.o = function(x){return(x[[1]])}
    }
  }else
  {
    if (extract)
    {
      stop('Cannot extract with more than 1 rv')
    }
    stop('Cannot handle more than 1 random variable set at a time.')
  }
  ctor.i <- distribution$outcome_ctor
  d = data.frame('x'=character(), 'p_x'=integer(), stringsAsFactors = FALSE)

  # Build the distribution
  for (i in 1:nrow(distribution$x))
  {
    outcome = unlist(distribution$x[i,])
    c.outcome = lapply(indexes, function(x){ctor.i(outcome[x])})
    c.outcome = ctor.o(c.outcome)
    d[i,]<- c(c.outcome, as.numeric(distribution$pmf[i]))
  }
  
  outcomes = unique(d[,1])
  pmf = aggregate(as.numeric(unlist(d[,2])), list(d[,1]), sum)$x
  
  sample_spaces = lapply(indexes, function(x){generate.sample.space(distribution$alphabet, ctor.i, x)})
  
  if (extract)
  {
    sample_space = sample_spaces[[1]]
  }else
  {
    sample_space = sample_spaces
  }
  new.distribution = create.distribution(as.data.frame(outcomes), pmf, sort=TRUE, sample.space=sample_space)
  return(new.distribution)
}

###############################################################################
# 3. Create marginal distribution from random variables

distribution.marginal<- function(distribution, rvs, rv.mode=NULL)
{
  rvs.index<- Parse.RVs(rvs, distribution, rv.mode, unique=TRUE, sort=TRUE)
  d = distribution.coalesce(distribution, list(rvs.index$indexes), rv.mode=1, extract=TRUE)
  # Set the random variable names  
  if (is.null(distribution$rvs))
  {
    names <- NULL
  }else
  {
    all_rvs = distribution$rvs
    rvs = all_rvs[c(rvs.index[,2])]
    d$rvs = rvs
  }
  L = distribution$outcome_length
  mask = rep(TRUE, L)
  mask[rvs.index[,1]] = FALSE
  d$mask = mask
  return(d)
}

###############################################################################
# 4. Modify outcomes of a distribution

# This may require combinations of pmfs when repeated values are found in the
# new outcomes.
# Input: ctor defines the function to use for transforming the outcomes
distribution.modify.outcomes<- function(distribution, ctor)
{
  outcomes<- distribution$x
  new.outcomes<- apply(outcomes, 1, ctor)
  if (ncol(new.outcomes) == nrow(outcomes))
  {
    new.outcomes<- t(new.outcomes)
  }
  new.outcome.names<- apply(new.outcomes, 1, distribution$outcome_ctor)
  new.dist = rep(0, length(new.outcome.names))
  names(new.dist)<- new.outcome.names
  pmf<- distribution$pmf
  for (i in 1:length(new.outcome.names))
  {
    new.dist[new.outcome.names[i]] = pmf[i] + new.dist[new.outcome.names[i]]
  }
  return(create.distribution(as.data.frame(new.outcomes), pmf=unname(new.dist)))
}
###############################################################################

# 5. Set names for random variables of a distribution

distribution.set.rv.names<- function(distribution, rv.names)
{
  if (is.null(rv.names))
  {
    # Explicit clearing of rv names
    rvs<- NULL
  }else
  {
    L = distribution$outcome_length
    if (length(unique(rv.names)) < L)
    {
      stop('Too few unique random variables')
    }else if (length(unique(rv.names)) > L)
    {
      stop('Too many unique random variables')
    }
    if (L>0)
    {
      rvs = rv.names
      distribution$rv_mode = 2
    }
    else
    {
      rvs = NULL
    }
    distribution$rvs = rvs
  }
  return(distribution)
}
###############################################################################