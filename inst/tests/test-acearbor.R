context("aceArbor works")
test_that("aceArbor marginal ancestral state reconstruction works", {
	require(testthat)
	data(anolis)
	anolis$dat$binary <- c(1,0)
	anolis$dat$char <- c("poop", "boobs")
	td<-make.treedata(anolis$phy, anolis$dat, name_column=1)
  tdDiscrete <- checkFactor(td)
  tdContinuous <- checkNumeric(td)
  ##Make sure checkNumeric and checkFactor properly convert the data:
  expect_true(all(sapply(tdDiscrete$dat, is.factor)))
	expect_true(all(sapply(tdContinuous$dat, is.numeric)))
  aaRes<-aceArbor(td, charType="continuous")
	#plot(aaRes, adj=c(0.5, 0.5))
  
	# since aceArbor uses fastAnc from phytools, we can compare to ape
	dd<-getVector(td$dat, 1)
	aceRes<-ace(dd, td$phy, method="ML")
	expect_equal(aaRes[[1]][,2], unname(aceRes$ace),  tolerance=0.000001)
  ddRes <- aceArbor(tdDiscrete, charType="discrete", aceType="marginal")
  #plot(ddRes, type="fan")
})