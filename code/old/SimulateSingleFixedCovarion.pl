# Simulate a covarion model over a tree
# With a fixed amino acid profile.
#
#
#
# Stephen Pollard 2015-6-4

use Math::Random;

use Tree;
use Samplers;

# Every node in a tree must be in a state. Along a branch, the state can
# change, just like amino acids. This makes sense because a node is an
# instant in time and so there are no state changes there.

# Given the variable state, what is the probability of substitution?

unless (caller) {
	my $tree_filename = $ARGV[0] // "../data/treeoutfile_averaged.newick";
	my $sequences_filename = $ARGV[1] //"../sites/site$site"
	  ."_subTime$model->{MeanSubTime}_switchTime$model->{MeanSwitchTime}_sequences.fasta";
	my $hiddenStates_filename = $ARGV[2] // "../sites/site$site"
	  ."_subTime$model->{MeanSubTime}_switchTime$model->{MeanSwitchTime}_hiddenStates.fasta";

	my $multiple_subs_filename = $sequences_filename . ".multiple_subs";
	SimulateSingleFixedCovarion(
		$tree_filename, $sequences_filename,
		$hiddenStates_filename, $multiple_subs_filename
	);
}

sub SimulateSingleFixedCovarion {
	my (
		$tree_filename, $sequences_filename,
		$hiddenStates_filename, $multiple_subs_filename
	  )
	  = @_;

	my $model = {
		MeanSwitchTime => 1,
		MeanSubTime => 1,
		Residues => ["A", "B", "C", "D"],
		ResidueProps => {A=>0.25, B=>.25, C=>.25, D=>0.25},
	};
	my $tree = Tree::FromFile($tree_filename);

	# fixed or variable
	$tree->{HiddenState} = Samplers::random_choice(["F", "V"]);

	# Force variable hidden state to test calculating the probabilities
	#$tree->{HiddenState} = "F";

	# residue
	$tree->{Sequence} = SampleNewResidue($model);
	
	open my $multiple_subs_file, ">", $multiple_subs_filename;

	# Force the initial residue to be A
	#$tree->{Sequence} = "A";
	SimulateCovarion($tree, $model, $multiple_subs_file);

	open my $sequences_file, ">", $sequences_filename;
	print $sequences_file Tree::AllSequencesToString($tree);

	open my $hiddenStates_file, ">", $hiddenStates_filename;
	print $hiddenStates_file Tree::AllHiddenStatesToString($tree);
}

sub SimulateCovarion {
	my ($tree, $model, $multiple_subs_file) = @_;

	# Simulate the evolution from ancestor to descendent
	Simulate($tree,  $model, $multiple_subs_file);

	SimulateCovarion($tree->{Left}, $model, $multiple_subs_file) if $tree->{Left};
	SimulateCovarion($tree->{Right}, $model, $multiple_subs_file)if $tree->{Right};
}

sub SampleNewResidue {
	my ($model, $current_residue) = @_;

	# Copy the residue propensities
	my $residue_props = {%{$model->{ResidueProps}}};

	if (defined $current_residue) {
		$residue_props->{$current_residue} = 0;
	}
	$residue_props = [ values %$residue_props];

	my $sampled_residue_index =Samplers::weighted_random_choice($residue_props);
	my $sampled_residue =$model->{Residues}->[$sampled_residue_index];

	return $sampled_residue;
}

sub Simulate {
	my ($descendant, $model, $multiple_subs_file) = @_;

	#print "Working with tree " . $descendant->{Name}, "\n";

	my $ancestor = $descendant->{Up};

	# Root does not have an ancestor
	return if not $ancestor;

	my $current_time = 0;
	my $current_state = $ancestor->{HiddenState};
	my $current_residue = $ancestor->{Sequence};

	my $done = 0;
	my $single_sub = 0;
	my $multiple_subs = 0;

	while (not $done) {
#		print "Starting new hidden state segment\n";
#		print "Current state $current_state\n";
#		print "Current residue $current_residue\n";
#		print "Current time $current_time\n";

		# Sample the next switch time
		my $sampled_time =
		  Math::Random::random_exponential(1, $model->{MeanSwitchTime});
#		print "Sampled time $sampled_time\n";

		if ($current_time + $sampled_time > $descendant->{Distance}) {
			$sampled_time = $descendant->{Distance} - $current_time;
			$done = 1;
		}
		$current_time += $sampled_time;
#		print "Time used $sampled_time\n";

		if ($current_state eq "V") {
#			print "Checking to see if a sub happened\n";
#			print "Time for sub to happen: $sampled_time\n";

			my $done2 = 0;
			my $current_time2 = 0;
			while (not $done2) {
#				print "Current time 2: $current_time2\n";

				# Is this the correct variance?
				my $sampled_time2 =
				  Math::Random::random_exponential(1,$model->{MeanSubTime});
				$current_time2 += $sampled_time2;
#				print "Sampled time 2: $sampled_time2\n";
#				print "Next sub could happen at time $current_time2\n";

				if ($current_time2 > $sampled_time) {
#					print(
#						"No more substitutions in this sampled variable time\n"
#					);
					$done2 = 1;

				}
				else {
					my $chosen_residue =
					  SampleNewResidue($model, $current_residue);
#					print
#					  "SUBSTITUTION from $current_residue to $chosen_residue\n";
					$current_residue = $chosen_residue;
					
					$multiple_subs = 1 if ($single_sub);
					$single_sub = 1 if (not $single_sub);
				}
			}
		}

		# Only flip the current hidden state if we haven't passed the end of the
		# branch.
		if (not $done) {
			$current_state =$current_state eq "V" ? "F" : "V";
		}
	}

	$descendant->{HiddenState} = $current_state;
	$descendant->{Sequence} = $current_residue;

#	print "Final state $current_state\nFinal residue $current_residue\n";
	
	if ($multiple_subs) {
		print $multiple_subs_file $descendant->{Name}, "\n";
	}
}

sub SampleState {
	my ($current_state, $propensities) = @_;

}

