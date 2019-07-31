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
	my $suboutfile = $ARGV[0] // "../Data/pipeline/suboutfile_highprob_ASRV_updatebls7_new_duplicated";
	my $burnin = $ARGV[1] // 1;
	my $gens = $ARGV[2] // 2;
	my $pThreshold = $ARGV[3] // 0.6;
	my $suboutfile_out = $ARGV[4] // $suboutfile . "_highprob";
	
	FilterByProbability($suboutfile, $burnin, $gens, $pThreshold, $suboutfile_out);
}

sub FilterByProbability {
	my ($suboutfile, $burnin, $gens, $pThreshold, $suboutfile_out) = @_;
	
	my $fh = initSubs($suboutfile, $burnin);
	
	my $branches = {};
	my $gen = $burnin;
	
	while ($gen < $gens + $burnin) {
		if (eof($fh)) {
			warn "Fewer post-burnin generations found than expected. \n"
			. "Expected: " . $gens . "\n"
			. "Found: " . ($gen - $burnin);
			last;
		}
		my $line = <$fh>;
		if ($line =~ /\/\//)
		  {
	
			  # print "Generation $gen" if $gen % 100 == 0;
			  $gen++;
			  next;
		}
	
		chomp $line;
	
		my $fields = [split " ", $line];
		my $tag = shift @$fields; # tag contains parent..child or the node name
		# e.g. 132..142 or Node_1
	
		# This is required for taxa that don't have any substitutions.
		# If this is not done, taxa without substitutions won't be printed
		if (@$fields) {
			  foreach my $sub (@$fields) {
				  $branches->{$tag}->{$sub}++;
			  }
		}
		else {
	
			  # This 'substitution' will never pass the threshold
			  $branches->{$tag}->{empty} = -1;
		}
	}
	
	print "Number of generations used to filter the substitutions: " . (
		  $gen - $burnin) . "\n";
	
	my $threshold = ($gen - $burnin) * $pThreshold;
	
	open my $fh_out, '>', $suboutfile_out;
	print $fh_out "Branch\tSubstitutions\n";
	foreach my $tag (keys %$branches) {
		  my $line = $tag;
		  while (my ($sub, $count) = each(%{$branches->{$tag}})) {
			  if ($count >= $threshold) { # This has always been greater than or
				  # equal to the threshold
				  $line .= "\t$sub";
			  }
		  }
		  print $fh_out $line . "\n";
	}
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

