
source("SetupSmallPDF.R")

main = function() {
	args <- commandArgs(trailingOnly = TRUE)
	
	# This is how I will tell if the script is being called from the command
	# line or if it is sourced into another script or from the R command line.
	#names(args) = c("file", "output", "with_wag", "with_stokes", "with_trend_line", "with_random", "line_colors", "colors")[1:length(args)]
	if (length(args) > 0) {
		#do.call(ColorDensityDataToPlotPDF, as.list(args))
	}
}



ColorDensityDataToPlotPDF = function(file, 
		output = sub(".density_data", ".pdf", file, fixed=T),
		with_wag = F,
		with_stokes = F,
		with_trend_line = T,
		with_random = F, 
		line_colors = list(line="#000000"),
		colors = colorRampPalette(c("#f0f0f0", "#555555"))(100)
) {
	if (!grepl(".density_data", file, fixed=T)) {
		stop("Input file format does not end with .density_data")
	}
	
	density_data = read.table(file)
	
	color_density_data = AssignColors(density_data, colors)
	SetupMitochondrialPDF(output)
	colorDensityPlot(color_density_data)
	
	if(with_random) 
		PlotTrendLine("../Data/NodePairInfo_randomized.trend_line", line_colors$random)
	if(with_stokes) 
		PlotTrendLine("NewSimulations/stokes_NodePairInfo_highprob_20.trend_line", line_colors$stokes)
	if(with_wag)
		PlotTrendLine("NewSimulations/wag_NodePairInfo_highprob_20.trend_line", line_colors$wag)
	if (with_trend_line){
		trend_line_file = sub(".density_data", ".trend_line", file, fixed=T)
		PlotTrendLine(trend_line_file, line_colors$line)
	}

	
	plotColorBar(min(color_density_data$counts), max(color_density_data$counts), 
			colors)
	
	dev.off()
}

ColorDensityDataToPlotSplitPDF = function(file, 
		with_wag = F,
		with_stokes = F,
		with_trend_line = F,
		with_random = F, 
		line_colors = list(line="#FF0000", stokes="#0000FF", wag="#000000", random="#555555"),
		colors = colorRampPalette(c("#f0f0f0", "#555555"))(100)
) {
	if (!grepl(".density_data", file, fixed=T)) {
		stop("Input file format does not end with .density_data")
	}
	
	density_data = read.table(file)
	
	color_density_data = AssignColors(density_data, colors)
	
	par(cex = 0.75/2)
	par(bg='white')
	par(family='serif')
	
	plot.new()
	par(mai=c(0.45, 0.55, 0.2, 0.1)/2) # absolute margins
	
	plot.window(xlim=c(0, 2.15), ylim=c(0, 2))
	
	colorDensityPlot(color_density_data)
	
	if(with_random) 
		PlotTrendLine("../Data/NodePairInfo_randomized.trend_line", line_colors$random)
	if(with_stokes) 
		PlotTrendLine("../Data/Stokes_cd.trend_line", line_colors$stokes)
	if(with_wag)
		PlotTrendLine("../Data/wag_cd.trend_line", line_colors$wag)
	if (with_trend_line){
		trend_line_file = sub(".density_data", ".trend_line", file, fixed=T)
		PlotTrendLine(trend_line_file, line_colors$line)
	}
	
	
	title(xlab="Distance Between Branches", cex.lab=1.25, line=1.5)
	title(ylab="Convergence / Divergence", cex.lab=1.25, line=2.25)
	
	axis(1, padj=-1.0)
	
	axis(2, las=1, hadj=0.75)
	
	plotColorBarSplit(min(color_density_data$counts), max(color_density_data$counts), 
			colors)
	
}


PlotTrendLine = function(file, color) {
	trend_line = read.table(file)
	lines(trend_line, col=color, lwd=2)
}

AssignColors = function(density_data, colors) {
	count_bins = head(seq(min(density_data$counts), max(density_data$counts), length.out=length(colors)+1), -1)
	
	# This determines the indices of the bin into which each count goes
	zi = findInterval(density_data$counts, count_bins)
	
	density_data$colors = colors[zi]
	return(density_data)
}


colorDensityPlot = function(color_density_data, r=1) { 
	symbols(color_density_data$x_rep, color_density_data$y_rep, 
			circles=color_density_data$radius, inches=F,
			bg=color_density_data$colors, fg=color_density_data$colors, 
			add=T)
}

plotColorBar <- function(min, max, colors) {
	# colors is a vector of colors
	
	paro <- par(no.readonly=T)
	
	
	N <- length(colors)
	
	x <- rep(.5, N) 
	dx <- rep(1, N)
	
	#y[1] <- min; y[N] <- max
	y <- seq(min, max, length.out=N)
	dy <- rep(y[2] - y[1], N)

	par(mai=c(0.45, 0.0, 0.2, 0.2))
	plot.new()
	
	par(cex=0.75)

	plot.window(xlim=c(0,1), ylim=c(min - dy[1] / 2, max + dy[1] / 2))
	symbols(x, y, rectangles=cbind(dx, dy), bg=colors, fg=colors, inches=F, add=T)
	axis(2, las=1)
	
	mtext("Density", side=1, line=1.5, cex=0.75)
	
	par(paro)
}

plotColorBarSplit <- function(min, max, colors) {
	# colors is a vector of colors
	
	paro <- par(no.readonly=T)
	
	
	N <- length(colors)
	
	x <- rep(.5, N) 
	dx <- rep(1, N)
	
	#y[1] <- min; y[N] <- max
	y <- seq(min, max, length.out=N)
	dy <- rep(y[2] - y[1], N)

	par(mai=c(0.4, 0.0, 0.2, 0.2)/2)
	plot.new()
	par(cex=0.5)
	plot.window(xlim=c(0,1), ylim=c(min - dy[1] / 2, max + dy[1] / 2))
	symbols(x, y, rectangles=cbind(dx, dy), bg=colors, fg=colors, inches=F, add=T)
	axis(2, las=1)
	
	mtext("Density", side=1, line=1, cex=0.75)
	
	par(paro)
}

main()