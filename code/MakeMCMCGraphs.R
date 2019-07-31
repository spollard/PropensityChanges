

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0){
    stop("Not enough arguments", call.=FALSE)
}
directory = args[1]
burnin = args[2]

mcmc = read.table(paste0(directory,"mcmc"), T)
mcmcb  = mcmc[-(1:burnin),]
model = read.table(paste0(directory,"model"), T)[-(1:burnin),]
modelb = model[-(1:burnin),]
optimals = read.table(paste0(directory,"optimals"), T)[-(1:burnin),]
optimals = optimals[-(1:burnin),]


png(paste0(directory, "likelihood_burn.png"), 1000,1000)
plot(mcmc$Current_log_likelihood)
png(paste0(directory, "likelihood_no_burn.png"), 1000,1000)
plot(mcmcb$Current_log_likelihood)

png(paste0(directory, "likelihood_posterior.png"), 1000,1000)
plot(density(mcmc$Current_log_likelihood))
png(paste0(directory, "likelihood_posterior_no_burn.png"), 1000,1000)
plot(density(mcmcb$Current_log_likelihood))


png(paste0(directory, "acceptables_size_prior.png"), 1000,1000)
plot(model$Prior_A)
png(paste0(directory, "acceptables_size_prior_no_burn.png"), 1000,1000)
plot(modelb$Prior_A)

png(paste0(directory, "acceptables_size_prior_density.png"), 1000,1000)
plot(density(model$Prior_A))
png(paste0(directory, "acceptables_size_prior_density_no_burn.png"), 1000,1000)
plot(density(modelb$Prior_A))


png(paste0(directory, "switch_prob.png"), 1000,1000)
plot(model$SwitchRate)
png(paste0(directory, "switch_prob_no_burn.png"), 1000,1000)
plot(modelb$SwitchRate)

png(paste0(directory, "switch_prob_density.png"), 1000,1000)
plot(density(model$SwitchRate))
png(paste0(directory, "switch_prob_density_no_burn.png"), 1000,1000)
plot(density(modelb$SwitchRate))

png(paste0(directory, "sub_prob.png"), 1000,1000)
plot(model$SubRate)
png(paste0(directory, "sub_prob.png_no_burn"), 1000,1000)
plot(modelb$SubRate)

png(paste0(directory, "sub_prob_density.png"), 1000,1000)
plot(density(model$SubRate))
png(paste0(directory, "sub_prob_density_no_burn.png"), 1000,1000)
plot(density(modelb$SubRate))


png(paste0(directory, "optimal_acceptable_prior.png"), 1000,1000)
plot(optimals$optimal_acceptable_prior)

png(paste0(directory, "optimal_sub_prob.png"), 1000,1000)
plot(optimals$optimal_sub_prob)

png(paste0(directory, "optimal_switch_prob.png"), 1000,1000)
plot(optimals$optimal_switch_prob)


