#' Add species into phylogeny near hypothesized relatives
#'
#' Given a data frame of species and taxonomic assignments, and an accepted phylogeny with
#' some of those species in it, will add the missing species in next to a taxonomic 
#' relative.
#'
#' @param tree An ape-style phylogenetic tree
#' @param groupings A data frame with two columns, "species" and "group". Missing species,
#' to be added, are taken as those that do not match a value in the tip labels of tree.
#' @param noRandomTrees The number of desired final trees with all missing species from 
#' groupings added.
#' @param printToScreen Default is TRUE. Will print a list of the missing taxa to screen.
#' @param saveToFile Default is FALSE. To automatically save multiPhylo object to file,
#' set to TRUE. The object will not be stored directly in memory, but it is calculated in
#' memory. The function could probably be sped up by having this save directly to file,
#' appending each tree instead of first saving them all to memory.
#' 
#' @details Given a data frame of two columns, "species" and "group", will take a species
#' that is absent from the phylogeny and bind it in as a sister to one of its taxonomic
#' relatives. A new node is created at half the distance between the existing species and
#' its most recent common ancestor. The new species is assigned
#' a branch length from the new node of half the original distance edge length of the
#' existing species to its original most recent common ancestor. Thus, if the input trees 
#' are ultrametric, the output trees should remain so. Currently, no effort is made to
#' ensure that the taxonomic groups of the missing species are actually to be found in the
#' species in the input tree.
#'
#' @return A multiPhylo object with number of trees as determined by noRandomTrees
#'
#' @export
#'
#' @references Eliot Miller unpublished
#'
#' @examples
#' #create a dummy tree
#' s <- "dummy(((Silly_aluco:4.2,Silly_otus:4.2):3.1,Lesssilly_noctua:7.3):6.3,Forreal_alba:13.5);"
#' cat(s, file = "ex.tre", sep = "\n")
#' dummy.tree <- read.tree("ex.tre")
#' #create a data frame of all taxa, including missing species
#' species <- c("Silly_aluco","Silly_otus","Silly_joke","Lesssilly_noctua","Forreal_alba","Notsosilly_blanca")
#' group <- c(rep("Silly_group", 3), "Lesssilly_group", "Forreal_group", "Forreal_group")
#' dummy.frame <- data.frame(species, group)
#' #run the function
#' test <- randomlyAddTaxa(tree=dummy.tree, groupings=dummy.frame, noRandomTrees=10)

randomlyAddTaxa <- function(tree, groupings, noRandomTrees, printToScreen=TRUE, saveToFile=FALSE)
{
	#set the species and group to character vectors
	groupings$species <- as.character(groupings$species)
	groupings$group <- as.character(groupings$group)

	#subset the groupings to those we have phylogenetic information for	
	realGroupings <- groupings[groupings$species %in% tree$tip.label, ]

	#and those we only have taxonomic information for
	possGroupings <- groupings[!(groupings$species %in% tree$tip.label), ]

	if(printToScreen==TRUE)
	{
		print("You are missing phylogenetic information for these species:", quote=FALSE)
		print(possGroupings$species)
	}

	#set up an empty list to save the results into
	random.trees <- list()
	for(i in 1:noRandomTrees)
	{
		new.tree <- tree
		for(j in 1:dim(possGroupings)[1])
		{
			#subset the real groupings to those instances where the group matches the sp
			#you are adding in. then randomly choose one of those species in that group in
			#the genetic phylogeny. then identify which "node" that tip is
			bindingToList <- realGroupings$species[realGroupings$group == possGroupings$group[j]]
			bindingToSpecies <- sample(bindingToList, 1)
			bindingTo <- which(new.tree$tip.label==bindingToSpecies)

			#this identifies the node that subtends the species you are binding to
			cameFrom <- new.tree$edge[,1][new.tree$edge[,2]==bindingTo]

			#set up this temporary matrix so you can give it row names and subset according
			#to the row name. this row name becomes index, which is what you use to subset
			#to find the right edge length
			tempMatrix <- new.tree$edge
			rownames(tempMatrix) <- 1:dim(tempMatrix)[1]
			index <- rownames(tempMatrix)[tempMatrix[,1]==cameFrom & tempMatrix[,2]==bindingTo]
			index <- as.numeric(index)

			#determine what the distance from that node to the tips is
			wholeDistance <- new.tree$edge.length[index]

			#use bind.tip to add in the tip
			new.tree <- bind.tip(tree=new.tree, tip.label=possGroupings$species[j], edge.length=wholeDistance/2, where=bindingTo, position=wholeDistance/2)
		}		
		#add the last version of new tree in as a new element in the list of random, final trees
		random.trees[[i]] <- new.tree
	}

	class(random.trees) <- "multiPhylo"

	if(saveToFile == FALSE)
	{
		return(random.trees)
	}
	else
	{
		write.tree(random.trees, file="randomized_trees.tre")
		print("Trees saved to working directory")
	}
}
