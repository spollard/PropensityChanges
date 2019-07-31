# Takes all substitutions from a suboutfile and randomly reassigns them to
# sites and branches that had a substitution. It first reads in all 
# substitutions, then randomly reassigns a start and end state from the list 
# of previously observed substitutions. 

# Could also randomize while keeping the initial amino acid the same. This 
# would begin to take into account the distance separating and would introduce 
# the bias of same starting state vs different. 


use strict; 
use warnings;

use Data::Dumper;

my $number_of_substitutions;
my $keep_ancestral_state;
unless(caller) {
	$keep_ancestral_state = 0;
	
	my $suboutfile = "suboutfile_highprob";
	my $suboutfile_out = "suboutfile_randomized";
	
	
	my $substitutions = ReadSubs($suboutfile);
	
	$number_of_substitutions = $#$substitutions + 1;
	print "number of subs $number_of_substitutions";
	
	RandomizeSubs($suboutfile, $substitutions, $suboutfile_out);
}

sub ReadSubs {
	my ($suboutfile) = @_;
	open my $fh, $suboutfile or die $!;
	<$fh>; # Burn header file
	
	my $substitutions = [];
	# First pass through suboutfile, collect probabilities of all substitutions
	while (my $line = <$fh>) { 
		chomp $line;
		
		my $fields = [split " ", $line];
		shift @$fields; # get rid of tag/branch id/node id
		
		foreach my $sub (@$fields) {
			push $substitutions, substr($sub, 0, 1) . substr($sub, -1);
		}
	}
	return $substitutions
}

sub RandomizeSubs {
	my ($suboutfile, $substitutions, $suboutfile_out) = @_;

	open my $fh_out, '>', $suboutfile_out;
	print $fh_out "Branch\tSubstitutions\n";
	
	# Reset file handle
	open my $fh, $suboutfile or die $!;
	<$fh>; # Burn header file
	
	# second pass, replace substitution start and end state with states randomly 
	# drawn from the substitutions array
	while (my $line = <$fh>) {
		chomp $line;
		
		my $fields = [split " ", $line];
		print $fh_out shift @$fields;
		
		foreach my $sub (@$fields) {
			my $rand_sub = $substitutions->[int(rand($number_of_substitutions))];
			if (not $keep_ancestral_state){
				substr($sub, 0, 1) = substr $rand_sub, 0, 1;
			}
			substr($sub, -1, 1) = substr $rand_sub, -1, 1;
			print $fh_out "\t$sub";
		}
		print $fh_out "\n";
	}

}