# Simulate a covarion model over a tree
# With a fixed amino acid profile.
#
#
#
# Stephen Pollard 2015-6-4

# Every node in a tree must be in a state. Along a branch, the state can
# change, just like amino acids. This makes sense because a node is an
# instant in time and so there are no state changes there.

# Given the variable state, what is the probability of substitution?

use lib 'C:\Users\Stephen\Documents\Lab work\PerlModules';
require "MultithreadCommands.pl";

unless (caller){

	my $sites_dir = $ARGV[0] // "../data/test/";
	my $number_of_sites = $ARGV[1] // 2;
	mkdir $sites_dir if not -d $sites_dir;

	my $commands = [];

	for my $site (1 .. $number_of_sites) {
		my $command = "perl FixedCovarionMCMC.pl $site";
		push @$commands, $command;
	}
	MultithreadCommands($commands);
}

