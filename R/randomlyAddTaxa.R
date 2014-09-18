#' Add species into phylogeny near hypothesized relatives
#'
#' Given a data frame of species and taxonomic assignments, and an accepted phylogeny with
#' some of those species in it, will add the missing species in next to a taxonomic 
#' relative.
#'
#' @param tree An ape-style phylogenetic tree
#' @param groupings A data frame with two columns, "species" and "group". Missing species,
#' to be added, are taken as those that do not match a value in the tip labels of tree.
#' @param fromNode Whether species should be added in a "polytomy" with, "crown" (more recently
#' diverged), "stem" (previously diverged), or "randomly" crown-wards or stem-wards from
#' the tip they are bound to.
#' @param noRandomTrees The number of desired final trees with all missing species from 
#' groupings added.
#' @param printToScreen Default is TRUE. Will print a list of the missing taxa to screen.
#' @param saveToFile Default is FALSE. To automatically save multiPhylo object to file,
#' set to TRUE. The object will not be stored directly in memory, but it is calculated in
#' memory. The function could probably be sped up by having this save directly to file,
#' appending each tree instead of first saving them all to memory.
#' 
#' @details Given a data frame of two columns, "species" and "group", will take a species
#' that is absent from the phylogeny and bind it to a randomly selected taxonomic 
#' relative. Note that the current implementation of the function iteratively builds up
#' the tree, meaning that a species being bound into the tree can be bound to a
#' species that was in the input phylogeny, or to a species that was added during a
#' previous iteration of the function. Previously, added species could only be bound to
#' species that were in the input phylogeny. This new feature slows the function down,
#' but it results in more realistic-looking, less ladderized trees. Four distinct methods
#' of adding new taxa are possible. With "polytomy", the new species is simply assigned to
#' a tip. Both species end up with branch lengths to their most recent common ancestor of
#' 0. With "crown", the new species is bound in at half the distance between the species
#' it is being bound to and that species' original parent node. With "stem", the new
#' species is bound to the parent node of the species it was selected to be bound to. The
#' new species is bound in at half the distance between the parent node and grandparent
#' node. With "randomly", each new species is randomly bound either crown-wards or
#' stem-wards, following the descriptions above. Note that even if "stem" or "randomly" 
#' are chosen, if a species is to be bound to the sister species to all others (the most
#' "basal" species in the phylogeny), it will automatically be bound stem-wards. With all
#' options, if the input tree is ultrametric, the output trees should remain so. 
#' Currently, no effort is made to ensure that the taxonomic groups of the missing species
#' are actually to be found in the species in the input tree.
#'
#' @return A multiPhylo object with number of trees as determined by noRandomTrees
#'
#' @export
#'
#' @references Eliot Miller unpublished
#'
#' @examples
#' #start libraries
#' library(phytools)
#'
#' #load a molecular tree up
#' data(bird.families)
#'
#' #create a data frame of all taxa from the phylogeny, and make up species groups
#' for each.
#' dummy.frame <- data.frame(species=bird.families$tip.label, 
#' group=c(rep("nonPasserine", 95), rep("suboscine", 9), rep("basalOscine", 13), 
#' rep("oscine", 20)))
#' 
#' #now make up a few passerine taxa that are missing and add these into the dummy frame
#' toAdd <- data.frame(species=c("Icteria", "Yuhina", "Pityriasis", "Macgregoria"), 
#' group=c(rep("oscine", 2), rep("basalOscine", 2)))
#' groupsDF <- rbind(dummy.frame, toAdd)
#'
#' #examples of changing the fromNode argument. you can plot or write these trees out to
#' #better see what the differences are
#' crown <- randomlyAddTaxa(tree=bird.families, groupings=groupsDF, fromNode="crown", noRandomTrees=10)
#' stem <- randomlyAddTaxa(tree=bird.families, groupings=groupsDF, fromNode="stem", noRandomTrees=10)
#' polytomy <- randomlyAddTaxa(tree=bird.families, groupings=groupsDF, fromNode="polytomy", noRandomTrees=10)

randomlyAddTaxa <- function(tree, groupings, fromNode, noRandomTrees, printToScreen=TRUE, saveToFile=FALSE)
{
	#this function requires the bind.tip and read.newick functions in phytools
	require(phytools)

	#set the species and group to character vectors
	groupings$species <- as.character(groupings$species)
	groupings$group <- as.character(groupings$group)

	#subset the groupings data frame to those spp we only have taxonomic information for
	possGroupings <- groupings[!(groupings$species %in% tree$tip.label), ]

	if(printToScreen==TRUE)
	{
		print("You are missing phylogenetic information for these species:", quote=FALSE)
		print(possGroupings$species)
	}

	#it is not easily possible to add a missing species below the root. find the root node
	#identity, the species that is connected to it, and the species group that species
	#belongs to. figure out how many species are in the input tree from that group. if
	#only 1, then make a list of species that have to be added above it, even if fromNode
	#is set to "stem" or "randomly"
	sisterToAll <- tree$tip.label[length(tree$tip.label)]
	sisterGroup <- groupings$group[groupings$species==sisterToAll]
	sisterGroupDF <- groupings[groupings$group==sisterGroup,]

	if(sum(!(sisterGroupDF$species %in% tree$tip.label)) == 1)
	{
		dangerList <- possGroupings$species[possGroupings$group==sisterGroup]
		print("The following taxon/taxa could potentially be bound to the taxon that is sister to all others. Cannot add to stem of this taxon, so will always add this/these to crown:", quote=FALSE)
		print(dangerList, quote=FALSE)
	}
	
	else
	{
		dangerList <- "NA"
	}

	#set up an empty list to save the results into
	random.trees <- list()
	for(i in 1:noRandomTrees)
	{
		new.tree <- tree
		for(j in 1:dim(possGroupings)[1])
		{
			#if fromNode was set to "randomly", define whether the current species to be
			#added will be added above or below the node
			if(fromNode=="randomly")
			{
				chooseFrom <- c("crown","stem")
				fromNodeUsed <- sample(chooseFrom,1)
			}
			
			else
			{
				fromNodeUsed <- fromNode
			}

			#add a line to automatically switch fromNode back to "crown" if the species
			#being added is on the danger list
			if(possGroupings$species[j] %in% dangerList & (fromNode=="stem" | fromNode=="randomly"))
			{
				fromNodeUsed <- "crown"
			}

			#subset the groupings to those we have phylogenetic information for	
			#note that this gets updated after each species is added in, so an unsequenced
			#species can get bound to what was also an unsequenced sp on a previous loop
			grouped <- groupings[groupings$species %in% new.tree$tip.label, ]

			#subset the species that are grouped to those whose group matches
			#the species you are adding in. 
			bindingToList <- grouped$species[grouped$group == possGroupings$group[j]]

			#check whether there is more than one species in the list. if so, set this
			#utility object to be 1. otherwise, set it to 0
			
			if(length(bindingToList) > 1)
			{
				bigEnough <- 1
			}
			else
			{
				bigEnough <- 0
			}

			#randomly choose one of the species that is currently in the phylogeny and is
			#in the same taxonomic group as the species being added.
			#identify which "node" that randomly chosen species is

			bindingToSpecies <- sample(bindingToList, 1)
			bindingTo <- which(new.tree$tip.label==bindingToSpecies)

			#identify the node that subtends the species you selected
			parent <- new.tree$edge[,1][new.tree$edge[,2]==bindingTo]

			#identify the node that subtends the parent node
			grandparent <- new.tree$edge[,1][new.tree$edge[,2]==parent]

			#set up a temporary matrix and give it row names. this allows you to pull out
			#the index of the edge in question, and subset the edge.lengths based on that
			#index, to get needed branch lengths later
			tempMatrix <- new.tree$edge
			rownames(tempMatrix) <- 1:dim(tempMatrix)[1]

			#use the matrix to get the index needed below
			parentIndex <- rownames(tempMatrix)[tempMatrix[,1]==parent & tempMatrix[,2]==bindingTo]
			grandparentIndex <- rownames(tempMatrix)[tempMatrix[,1]==grandparent & tempMatrix[,2]==parent]
			parentIndex <- as.numeric(parentIndex)
			grandparentIndex <- as.numeric(grandparentIndex)
			
			#find the edge length between the tip species you randomly selected and
			#its parent node, and the edge length between the parent and grandparent node
			parentDistance <- new.tree$edge.length[parentIndex]
			grandparentDistance <- new.tree$edge.length[grandparentIndex]

			#use bind.tip to add in the tip. if only 1 taxon is present in the phylogeny
			#for the species group in question, the new species has to be bound directly
			#to this to maintain monophyly			
			if(bigEnough == 0)
			{
				new.tree <- bind.tip(tree=new.tree, tip.label=possGroupings$species[j], where=bindingTo, position=parentDistance/2)
			}

			#bind new species directly to randomly selected species if "crown" is selected
			#give the new species a branch length of half the original distance of the
			#selected species and its parent node
			else if(bigEnough == 1 & fromNodeUsed == "crown")
			{
				new.tree <- bind.tip(tree=new.tree, tip.label=possGroupings$species[j], where=bindingTo, position=parentDistance/2)
			}
			
			#if more than one species is present in group, and if "polytomy" is selected
			#bind the new species to the parent node and do not assign it any position
			
			else if(bigEnough == 1 & fromNodeUsed == "polytomy")
			{
				new.tree <- bind.tip(tree=new.tree, tip.label=possGroupings$species[j], where=bindingTo)
			}

			#bind new species to parent node if "stem" is selected. give it an additional
			#branch length of half the distance to the grand parent node (phytools will
			#go ahead and make it ultrametric)
			else if(bigEnough == 1 & fromNodeUsed == "stem")
			{
				new.tree <- bind.tip(tree=new.tree, tip.label=possGroupings$species[j], where=parent, position=grandparentDistance/2)
				workAround <- write.tree(new.tree)
				new.tree <- read.newick(text=workAround)
			}
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
