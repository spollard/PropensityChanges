
SetupSmallPDF = function(file) {
	pdf(file=file, width=3.25, height=2.5, useDingbats=F)
	par(cex = 0.75)
	par(bg='white')
	par(family='serif')
	
	plot.new() 
	par(mai=c(0.45, 0.55, 0.2, 0.2)) # absolute margins
	
	par(xlog=T)
		plot.window(xlim=c(0.001, 10), ylim=c(0, 2))
	
	title(xlab="Distance Between Branches", cex.lab=1.25, line=1.5)
	title(ylab="Convergence / Divergence", cex.lab=1.25, line=2.25)

# With a tick at 'zero' which is actually 0.0005
	#axis(1, at=c(0.0005, 0.001, 0.01, 0.1, 1, 10), 
	#		labels=c(NA, "0.001", "0.01", "0.1", "1.0", "10"), padj=-1.0)
	
# Without a zero tick
	axis(1, at=c(0.001, 0.01, 0.1, 1, 10), 
			labels=c("0.001", "0.01", "0.1", "1.0", "10"), padj=-1.0)
	
	axis(2, las=1, hadj=0.75)
}

SetupSmallLinearPDF = function(file, ylim=NA) {
	pdf(file=file, width=3.25, height=2.5, useDingbats=F)
	par(cex = 0.75)
	par(bg='white')
	par(family='serif')
	
	plot.new() 
	par(mai=c(0.45, 0.55, 0.2, 0.2)) # absolute margins
	
	if (is.na(ylim)) {
		plot.window(xlim=c(0, 2.1), ylim=c(0, 2))
	} else {
		plot.window(xlim=c(0, 2.1), ylim=ylim)
	}
	
	title(xlab="Distance Between Branches", cex.lab=1.25, line=1.5)
	title(ylab="Convergence / Divergence", cex.lab=1.25, line=2.25)
	
	axis(1, padj=-1.0)
	
	axis(2, las=1, hadj=0.75)
}


SetupMitochondrialPDF = function(file) {
	pdf(file=file, width=3.25, height=2.5, useDingbats=F)
	
	layout(matrix(1:2, 1, byrow=T), widths=c(12,1))

	par(cex = 0.75)
	par(bg='white')
	par(family='serif')

	par(mai=c(0.45, 0.55, 0.2, 0.0))
	plot.new()
	
	plot.window(xlim=c(0, 2.2), ylim=c(0, 2))


	title(xlab="Distance Between Branches", cex.lab=1.25, line=1.5)
	title(ylab="Convergence / Divergence", cex.lab=1.25, line=2.25)
	
	axis(1, padj=-1.0)
	
	axis(2, las=1, hadj=0.75)
}


