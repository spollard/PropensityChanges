

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0){
    #stop("Not enough arguments", call.=FALSE)
}
directory = args[1]
burnin = args[2]

a = read.table("../data/simulated_100_sites_100_taxa_g1000_sampled/probs")
png("../data/simulated_100_sites_100_taxa_g1000_sampled/probs.png", 600,600)
plot(density(a$V1, bw=0.005), xlab="Posterior", ylab="Frequency", main="")
