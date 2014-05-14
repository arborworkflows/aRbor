context("aceArbor works")
test_that("aceArbor marginal ancestral state reconstruction works", {
	require(testthat)
	data(anolis)
	td<-make.treedata(anolis$phy, anolis$dat, name_column=1)
	td2<-select(td, SVL)
	aaRes<-aceArbor(td2, charType="continuous")
	
	# since aceArbor uses fastAnc from phytools, we can compare to ape
	dd<-td$dat[,1]
	names(dd)<-rownames(td$dat)
	aceRes<-ace(dd, td$phy, method="ML")
	expect_equal(aaRes$ancestralStates[,1], aceRes$ace,  tolerance=0.000001)
	
})