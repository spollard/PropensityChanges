
source("NodePairInfoToColorDensityData.R")
source("ColorDensityDataToPlotPDF.R")


# There will be a single command line interface for this script
# It will turn a node pair info file into a graph of C/D vs distance


args <- commandArgs(trailingOnly = TRUE)

npi_filename = args[1]
graph_filename = args[2]
anc = args[3] 
thresh = as.numeric(args[4])
 
 
 
density_filename = paste0(npi_filename,".density_data")

NodePairInfoToDensityData(npi_filename, output=density_filename, anc=anc, thresh=thresh)

colors = colorRampPalette(c("#8f9fff", "#00148f"))(100)
	
ColorDensityDataToPlotPDF(density_filename, output=graph_filename)

print(paste("Saved graph in", graph_filename))


