phyloSimulate<-function(phy, branchFUN, rootFUN) {
	res<-rootFUN() + branchFUN(1, 1)
	return(res)
}

branchPlus<-function(initValue, branchLength) {
	res<-initValue + branchLength
	return(res)
}

root10<-function() {
	res<-10
	return(res)
}
phy<-pbtree(n=100, scale=1)

phyloSimulate(phy, branchPlus, root10)