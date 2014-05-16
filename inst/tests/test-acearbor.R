context("aceArbor works")
test_that("aceArbor marginal ancestral state reconstruction works", {
	require(testthat)
	data(anolis)
	td<-make.treedata(anolis$phy, anolis$dat, name_column=1)
	aaRes<-aceArbor(td, charType="continuous")
	
	# since aceArbor uses fastAnc from phytools, we can compare to ape
	dd<-getVector(td$dat, 1)
	aceRes<-ace(dd, td$phy, method="ML")
	expect_equal(aaRes[[1]]$ancestralStates[,2], aceRes$ace,  tolerance=0.000001)
	
})