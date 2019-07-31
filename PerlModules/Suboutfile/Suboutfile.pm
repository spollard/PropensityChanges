
use strict;
use warnings;

package Suboutfile;

sub IsOld {
	my ($suboutfile) = @_;
	open my $substitutions_in, "<", $suboutfile or die $!;
	
	scalar <$substitutions_in>; # Burn header line
	my $second_line = scalar <$substitutions_in>;
	my $is_old = ($second_line =~ /\.\./);
	return $is_old;
}

sub ContainsMultipleGenerations {
	my ($suboutfile) = @_;
	open my $substitutions_in, "<", $suboutfile or die $!;
	
	my $contains_multiple_generations = 0;
	my $generation_separators = 0;
	while (<$substitutions_in>) {
		if ($_ eq "//\n") {
			$generation_separators++;
			if ($generation_separators >= 2) {
				$contains_multiple_generations = 1;
			};
		}
	}
	return $contains_multiple_generations;
}


sub Initialize {
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