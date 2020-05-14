############################################################################################################
# SHANNON INFORMATION PARTITIONS: by Martin Hilbert & Mahima Agarwal @UCDavis.edu (unreleased version 0.0) #
############################################################################################################
#  It was adopted from an earlier Python implementation by Ryan James: https://github.com/Autoplectic/dit 
############################################################################################################
# This file contains functions required for calculation of Shannon Information partitions
# and functions for working with and printing these partitions
###############################################################################

# 1. Create poset lattice for all random variables

Poset.Lattice<- function(elements)
{
  child = function(a,b){(sum(!(a %in% b)) == 0) && (length(b) - length(a) == 1)}
  all.combinations = combn(power.set(elements), 2, simplify=FALSE)
  edges.index<- which(sapply((all.combinations), function(x){child(x[[1]], x[[2]])}))
  edge.list<- matrix(rep('()'), ncol=2, nrow=length(edges.index))
  for (i in (1:length(edges.index)))
  {
    ind <- edges.index[i]
    a<- paste(paste(all.combinations[[ind]][[1]], collapse = ','))
    b<- paste(paste(all.combinations[[ind]][[2]], collapse = ','))
    if (a == '')
    {
      a = '()'
    }
    if (b == '')
    {
      b = '()'
    }
    edge.list[i,]<- c(b, a)
  }
  
  lattice = graph_from_edgelist(edge.list)
  return(lattice)
}

###############################################################################

# 2. Calculate shannon entropies

Shannon.Partition<- function(distribution)
{
  part<-list()
  part$measure<- Shannon.Entropy
  part$unit<- 'bits'
  part$dist<- distribution
  rvs<- distribution$rvs
  if (is.null(rvs))
  {
    rvs = seq(distribution$outcome_length)
  }
  part$lattice = Poset.Lattice(rvs)
  part$rlattice<- reverse(part$lattice)
  vertices<- V(part$lattice)$name
  Hs = rep(0, length(vertices))
  names(Hs)<- vertices
  Is = rep(0, length(vertices))
  names(Is)<- vertices
  rv.mode<- RV.MODES[[distribution$rv_mode]]
  get.string<- rv.mode == RV.MODES$NAMES
  # Calculate entropies
  for (i in 1:length(vertices))
  {
    Hs[i] = part$measure(part$dist, convert.vertex.name(vertices[i], str = get.string), distribution$rv_mode)
  }
  #Calculate co-information
  for (i in 1:length(vertices))
  {
    dfs.res = na.omit(dfs(part$lattice,vertices[i], unreachable = FALSE)$order)
    H = Hs[dfs.res]
    len.rv = sapply(names(dfs.res), convert.vertex.name, TRUE, get.string)
    Is[i] = sum((-1) ^ (len.rv + 1) * H)
  }
  # Mobius inversion of the above, resulting in the Shannon atoms.  
  topo.sort<- names(topo_sort(part$lattice, 'in'))
  atoms<- rep(0, length(topo.sort) - 1)
  names(atoms)<- topo.sort[-1]
  for (i in (length(topo.sort):2))
  {
    node<- topo.sort[i]
    kids<- names(na.omit(dfs(part$rlattice, node, unreachable = F)$order)[-1])
    atoms[node]<- Is[node] - sum(atoms[kids])
  }
  new.atoms = matrix(NA, nrow = length(atoms), ncol=4)
  rownames(new.atoms) = names(atoms)
  for (i in names(atoms))
  {
    rv<- convert.vertex.name(i, str=get.string)
    crvs<- sort(rvs[which(!(rvs %in% convert.vertex.name(i, str=get.string)))])
    a.rvs<- paste(rv, collapse=':')
    a.crvs<- paste(crvs, collapse = ', ')
    sym<- symbol(rv, crvs)
    new.atoms[i,]<- c(sym, a.rvs, a.crvs, atoms[i])
  }
  rownames(new.atoms)<- NULL
  part$atoms = new.atoms
  return(part)
}

###############################################################################

# 3. Print table of entropies from list from Shannon.Partition function

print.table<- function(partition, digits=3)
{
  atoms<- partition$atoms
  measures<- str_replace(paste(atoms[,1], '[', atoms[,2], '|', atoms[,3], ']', sep=''), '\\|]', ']')
  values<- formatC(round(as.numeric(atoms[,4]), digits), format='f', digits=digits, flag='0')
  colnames<- c('measure', partition$unit)
  table<- matrix(c(measures, values), ncol=2)
  pandoc.table(table, style='grid', col.names= colnames, split.cells=100)
}

###############################################################################  

# 4. Print subsets of partitions from list from the Shannon.Partition function

get.entropy.subset<- function(partition, var, cond.var, digits=3)
{
  atoms<- partition$atoms
  rows.with.var<- sapply(strsplit(atoms[,2], ':'), function(x){var %in% x})
  rows.with.cond.var<- sapply(strsplit(atoms[,3], ','), function(x){cond.var %in% trimws(x)})
  row.indexes<- as.logical(rows.with.cond.var * rows.with.var)
  atoms.table<- atoms[row.indexes,]
  emph.row <- nrow(atoms.table) + 1
  sum.entropy<- c('H', var, cond.var, sum(as.numeric(atoms[row.indexes,4])))
  atoms.table<- rbind(atoms.table, sum.entropy)
  measures<- str_replace(paste(atoms.table[,1], '[', atoms.table[,2], '|', atoms.table[,3], ']', sep=''), '\\|]', ']')
  values<- formatC(round(as.numeric(atoms.table[,4]), digits), format='f', digits=digits, flag='0')
  colnames<- c('measure', partition$unit)
  table<- matrix(c(measures, values), ncol=2)
  pandoc.table(table, style='grid', col.names= colnames, split.cells=100, emphasize.strong.rows=emph.row)
}
###############################################################################