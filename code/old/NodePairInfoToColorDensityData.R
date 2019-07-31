# This might be easier via densCols()
#x  = matrix(rnorm(1e2), ncol = 2)

#dcols = densCols(x)
#plot(x, col = dcols, pch = 20)


NodePairInfoToDensityData = function(file, 	
	ancestor = 'all',  # Can be 'all', 'same', 'different'
	threshold = 0,
	Nx = 60,
	Ny = Nx,
	xlim = c(0, 2), # Can be c(min, max)
	ylim = c(0, 1), 
	sim = F, # Ideally, this should be taken into account earlier. 
	output=paste0(file, ".density")) {

	NodePairInfo = read.table(file=file, header=T)
	print(paste("number of rows of node pair info", nrow(NodePairInfo)))

	if (ancestor == 'all') {
		convergence = NodePairInfo$Convergence
		divergence = NodePairInfo$Divergence
	}
	else if (ancestor == 'same') { 
		convergence = NodePairInfo$CAConvergence
		divergence = NodePairInfo$CADivergence
	}
	else if (ancestor == 'different') { 
		convergence = NodePairInfo$DAConvergence
		divergence = NodePairInfo$DADivergence
	}
	else stop("Ancestor type not recognized") 
	
	# This has always been greater than. It has never been greater than or 
	# equal to.
	indices_to_keep = divergence > threshold
	if (sim) {
		# Must convert time from nucleotide substitutions to amino acid 
		# substitutions
		# The reported 'distance between' is actually one half the distance 
		# between. The multiple of * 2 corrects for this.
# the distance between issue has been resolved in the new simulations
		x = NodePairInfo$Distance[indices_to_keep] * 0.695
	} else {
		x = NodePairInfo$Distance[indices_to_keep]
	}
	y = convergence[indices_to_keep] / (divergence[indices_to_keep])

	# Apply x limits
	indices_to_keep = (x <= xlim[2]) 
	x = x[indices_to_keep]
	y = y[indices_to_keep]



	l = GetBinsAndAssignments(x, Nx, xlim)
	xbins = l$bins
	xi = l$assignments

	l = GetBinsAndAssignments(y, Ny, ylim)
	ybins = l$bins
	yi = l$assignments
	
		
	# This counts the number of points in each 2d bin
	counts = table(factor(xi, levels=1:Nx), factor(yi, levels=1:Ny))
	

	dx = diff(xlim) / Nx
	dy = diff(ylim) / Ny
	
	y_means = tapply(y, xi, mean)
	
	occupied_xbins = xbins[as.numeric(names(y_means))]
	trend_line = data.frame(x = occupied_xbins + dx / 2, y=y_means)
	
	write.table(trend_line, file=paste0(file, ".trend_line"))

	# Construct corresponding x and y matrices 
	x_matrix = matrix(rep(xbins + dx / 2, Nx), nrow=Nx)
	y_matrix = matrix(rep(ybins + dy / 2, 1, each=Ny), nrow=Nx)

	# matrices to keep track of an representative point in a box
	x_rep = matrix(rep(NA, Nx * Ny), nrow=Nx)
	y_rep = matrix(rep(NA, Nx * Ny), nrow=Nx)
	
	# Collect the representative points
	for ( i in seq_along(xi)) {
		
		if (x[i] <= xlim[2] && y[i] <= ylim[2]) {
			x_rep[xi[i], yi[i]] = x[i]
			y_rep[xi[i], yi[i]] = y[i]
		}	
	}
	
	density_data = data.frame(
		counts = as.vector(counts),
		x_rep = as.vector(x_rep), 
		y_rep = as.vector(y_rep), 
		x_matrix = as.vector(x_matrix),
		y_matrix = as.vector(y_matrix), 
		radius = rep(dx/2, length(counts))
	)
	density_data = density_data[ density_data$counts > 0, ]

	write.table(density_data, file=output)
}

GetBinsAndAssignments = function(x, N, lims) {
	if (!is.vector(lims)) {
		print("No lims!")
		lims = range(x)
	}
	bins = head(seq(lims[1], lims[2], length.out=N+1), -1)
	assignments = findInterval(x, bins)
	return(list(bins=bins, assignments=assignments))
}
