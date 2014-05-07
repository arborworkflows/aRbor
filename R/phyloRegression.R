phyloRegression<-function(phy, x, y, picPlot=F, stdPlot=T) {
	dd<-data.frame(x, y)
	res<-gls(y~x, correlation=corBrownian(1, phy), data=dd)
	
	if(picPlot==T) {
		picx<-pic(x, phy)
		picy<-pic(y, phy)
		picfm<-lm(picy~picx+0)
		plot(picx, picy, pch=19)
		lx<-c(min(picx), max(picx))
		ly<-c(min(picx)*picfm$coefficients, max(picx)*picfm$coefficients)
		lines(lx, ly, lwd=2, col="red")
	}
	
	if(stdPlot==T) {
		plot(x, y, pch=19)
		pp<-predict(res)
		wmin<-which(x==min(x))
		wmax<-which(x==max(x))
		lx<-c(x[wmin], x[wmax])
		ly<-c(pp[wmin], pp[wmax])
		lines(lx, ly, lwd=2, col="red")


	}
	return(res)
}