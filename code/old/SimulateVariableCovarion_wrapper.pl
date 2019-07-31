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

require "MultithreadCommands.pl";

unless (caller){

	my $sites_dir = $ARGV[0] // "../data/var_sites1/";
	my $number_of_sites = $ARGV[1] // 100;
	mkdir $sites_dir if not -d $sites_dir;

	my $commands = [];

	for my $site (1 .. $number_of_sites) {

		my $tree_filename =  "../data/treeoutfile_averaged.newick";
		my $sequences_filename =
		  $sites_dir . "site$site"."_subTime1_switchTime1_sequences.fasta";
		my $hiddenStates_filename =
		  $sites_dir . "site$site"."_subTime1_switchTime1_hiddenStates.fasta";

		my $command =
"perl SimulateSingleVariableCovarion.pl $tree_filename $sequences_filename $hiddenStates_filename";
		push @$commands, $command;
	}
	MultithreadCommands($commands);
}

