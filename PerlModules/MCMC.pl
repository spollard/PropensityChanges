# Refactor out the common MCMC functions.
#
# Stephen Pollard
# 10/21/2014

use strict;
use warnings FATAL => 'all';

sub InitializeMcmcOutput {
	my ($mcmc_out_filename) = @_;
	open my $mcmc_out, ">", $mcmc_out_filename;

	print $mcmc_out join("\t",
		"Generation", "Current_log_likelihood",
		"Proposed_log_likelihood", "Delta_log_likelihood",
		"Acceptance_Probability", "Proposal_accepted", "Parameters_sampled")
	  . "\n";

	return $mcmc_out;
}

sub RecordMcmcState {
	my ($mcmc_out, $generation, $current_log_likelihood,
		$proposed_log_likelihood,$is_accepted, $number_of_parameters_sampled)
	  = @_;
      
    my $acceptance_prob = exp($proposed_log_likelihood - $current_log_likelihood);
	print $mcmc_out join("\t",
		$generation,
		$current_log_likelihood,
		$proposed_log_likelihood,
		$proposed_log_likelihood - $current_log_likelihood,
        $acceptance_prob,
		int($is_accepted),
		$number_of_parameters_sampled)
	  . "\n";

}

# Use a simple Metropolis-Hastings
# prob of accepting proposal = min (1, proposed_likelihood / current_likelhood) 
sub IsProposedLogLikelihoodAccepted {
	my ($current_log_likelihood, $proposed_log_likelihood) = @_;

#	print "Acceptance ratio: " . $proposed_log_likelihood / $current_log_likelihood;
# Cannot take log of 0 so
# Accept if rand == 0
	my $r = rand;
	if ($r == 0) {return 1}

	return log($r) < $proposed_log_likelihood - $current_log_likelihood;
}

1;
