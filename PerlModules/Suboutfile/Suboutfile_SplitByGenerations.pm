# Takes all substitutions from a suboutfile and gives data frame of substitions
# including probability. Can include a probability threshold. Choose 0 for no
# threshold.
#
# The posterior probability of the substitutions must be AT LEAST the
# threshold. Substitutions are kept if their posterior probability is greater
# than OR EQUAL TO the threshold.

# Stephen Pollard
# 5/20/2014

use strict;
use warnings;

package Suboutfile;

unless (caller) {
	my $suboutfile = $ARGV[0]
	  // "../Data/pipeline/suboutfile_highprob_ASRV_updatebls7_new_duplicated";
	my $burnin = $ARGV[1] // 2;
	my $gens = $ARGV[2] // 2;
	my $split_dir = $ARGV[4] // $suboutfile . "_splits/";

	mkdir $split_dir if not -d $split_dir;

	SplitByGenerations($suboutfile, $burnin, $gens, $split_dir);
}

sub SplitByGenerations {
	my ($suboutfile, $burnin, $gens, $split_dir) = @_;

	my $fh = initSubs($suboutfile, $burnin);

	my $gen = $burnin;

	open my $split_out, ">", $split_dir . $gen;
	print $split_out "Node\tSubstitutions\n";
	while ($gen < $gens + $burnin) {
		if (eof($fh)) {
			warn "Fewer post-burnin generations found than expected. \n"
			  . "Expected: "
			  . $gens . "\n"
			  . "Found: "
			  . ($gen - $burnin);
			last;
		}
		my $line = <$fh>;
		if ($line =~ /\/\//) {

			# print "Generation $gen" if $gen % 100 == 0;
			$gen++;

			# Don't make a new file for the end-of-last generation
			if ($gen < $gens + $burnin){
				close $split_out;
				open $split_out, ">", $split_dir . $gen;
				print $split_out "Node\tSubstitutions\n";
			}
			next;
		}

		# Print out line
		print $split_out $line;
	}

	print "Number of generations used to when splitting the substitutions: "
	  . ($gen - $burnin) . "\n";

}

sub initSubs {
	my ($file, $burnin) = @_;
	open my $fh, "<", $file or die $! . $file;
	my $g = 0;
	my $line;
	<$fh>; # burn header line
	while ($g != $burnin and defined($line = <$fh>)) {
		if ($line =~ /\/\//) {
			$g++;
		}

	}
	die "End of burnin not reached; check burnin" if $g != $burnin;
	return $fh;
}

1;

