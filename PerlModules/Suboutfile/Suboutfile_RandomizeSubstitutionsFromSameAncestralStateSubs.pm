# Takes all substitutions from a suboutfile and randomly reassigns them to
# sites and branches that had a substitution. It first reads in all 
# substitutions, then randomly reassigns an end state from the list 
# of previously observed end states for a given start state.  


use strict; 
use warnings;

package Suboutfile;

unless(caller) {
	my $suboutfile = $ARGV[0] // "C:/Users/Stephen/Documents/Lab work/Projects/Protein Evolution/Convergence/Data/suboutfile_highprob_ASRV_updatebls7_new_T";
	my $suboutfile_out = $ARGV[0] // "C:/Users/Stephen/Documents/Lab work/Projects/Protein Evolution/Convergence/Data/suboutfile_highprob_ASRV_updatebls7_new_T_randomized";
	
	RandomizeSubs($suboutfile, $suboutfile_out);
}

sub RandomizeSubs {
	my ($suboutfile,$suboutfile_out) = @_;
	
	my $substitutions = ReadSubs($suboutfile);
	
	open my $fh_out, '>', $suboutfile_out;
	print $fh_out "Branch\tSubstitutions\n";
	
	# Reset file handle
	open my $fh, $suboutfile or die $!;
	<$fh>; # Burn header line
	
	# Second pass, replace substitution end state with random end state 
	# drawn from the possible descendents for a given start state
	while (my $line = <$fh>) {
		chomp $line;
		
		my $fields = [split " ", $line];
		print $fh_out shift @$fields;
		
		# This loop is different
		foreach my $sub (@$fields) {
			my $ancestral_state = substr($sub, 0, 1);
			my $number_of_possible_subs = scalar @{$substitutions->{$ancestral_state}};
			my $rand_descendent_state = $substitutions->{$ancestral_state}->[int(rand($number_of_possible_subs))];
			substr($sub, -1) = $rand_descendent_state;
			print $fh_out "\t$sub";
		}
		print $fh_out "\n";
	}
}

# This name will conflict with the other readsubs in randomizesubstitutions.
sub ReadSubs {
	my ($suboutfile) = @_;
	
	open my $fh, $suboutfile or die $!;
	<$fh>; # Burn header line
	
	my $substitutions = {};
	my $number_of_substitutions = 0;
	# First pass through suboutfile, collect probabilities of all substitutions
	while (my $line = <$fh>) { 
		chomp $line;
		
		my $fields = [split " ", $line];
		shift @$fields; # get rid of tag/branch id/node id
		
		foreach my $sub (@$fields) {
			# This line is different
			push @{$substitutions->{substr($sub, 0, 1)}}, substr($sub, -1);
			$number_of_substitutions++;
		}
	}
	
	print "number of subs $number_of_substitutions";
	return $substitutions;
}