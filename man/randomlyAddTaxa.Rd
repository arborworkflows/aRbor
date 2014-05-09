\name{randomlyAddTaxa}
\alias{randomlyAddTaxa}
\title{Add species into phylogeny near hypothesized relatives}
\usage{
  randomlyAddTaxa(tree, realGroupings, possGroupings,
    noRandomTrees, saveToFile = FALSE)
}
\arguments{
  \item{tree}{An ape-style phylogenetic tree}

  \item{realGroupings}{A data frame with values of class
  character. First column, "species", must match a value in
  tree. Second column, "group", provides the taxonomic
  assignment of that species.}

  \item{possGroupings}{A data frame with values of class
  character. First column, "species", are tips to be added
  into the input tree. Second column, "group", provides the
  taxonomic assignment of that species.}

  \item{noRandomTrees}{The number of desired final trees
  with all species from possGroupings added.}

  \item{saveToFile}{Default is FALSE. To automatically save
  the multiPhylo object to file set to TRUE. The object
  will not be stored directly in memory, but it is
  calculated in memory. The function could probably be sped
  up by having this save directly to file, appending each
  tree instead of first saving them all to memory.}
}
\value{
  A multiPhylo object with number of trees as determined by
  noRandomTrees
}
\description{
  Given two data frames, one of species in an accepted
  phylogeny and their corresponding taxonomic assignment,
  and another of species to be bound in, along with their
  taxonomic assignment, will bind the missing species in
  next to a taxonomic relative.
}
\details{
  Given two data frames with values class CHARACTER (it is
  critical they are not of class factor), both with column
  names c("species", "group"), will take a species, e.g.,
  spA from the possGroupings data frame, find a species in
  the realGroupings frame, e.g., spB with the the same
  taxonomic group, and bind spA into the tree by creating a
  new node at half the distance between spB and its most
  recent common ancestor. spA is assigned a branch length
  from the new node of half the original distance of the
  spB to its original most recent common ancestor. Thus, if
  the input trees are ultrametric, the output trees should
  also remain ultrametric. Currently, no effort is made to
  ensure that the input data frames truly are of class
  character, nor is any effort made to confirm that all
  taxonomic groups from possGroupings are to be found in
  realGroupings. These checks should be built in at some
  point.
}
\references{
  Eliot Miller unpublished
}

