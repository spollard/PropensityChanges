# Takes all substitutions from a suboutfile and calculates their associations.
# That is the joint probability of all subs with all other subs.
# I an get the total number of subs at a site from a different script.

use strict;
use warnings;

use Data::Dumper;

unless(caller) {
	my $subs = [];
	for my $char1 (split "", "ACDEFGHIKLMNPQRSTVWY"){
		for my $char2 (split "", "ACDEFGHIKLMNPQRSTVWY"){
			push @$subs, $char1.$char2;
		}
	}
	print join("\t", "generation",@$subs) . "\n";

	for my $suboutfile (<../suboutfile_splits/*_new>){
		$suboutfile =~ /(\d*)/;
		next if not $1;
		my $gen = $1;
		next if $gen < 1000;

		#my $suboutfile ="../suboutfile_splits/1000_new";

		my $substitution_count_matrix =
		  CalculateSubstitutionCountMatrix($suboutfile);

		print($gen . "\t" .
			SubstitutionCountMatrixToString($substitution_count_matrix). "\n");
	}

}

sub CalculateSubstitutionCountMatrix {
	my ($suboutfile) = @_;
	open my $fh, $suboutfile or die $!;
	<$fh>; # Burn header line

	my $substitution_count_matrix = {};

	while (my $line = <$fh>) {
		chomp $line;

		my $fields = [split " ", $line];
		shift @$fields; # get rid of tag/branch id/node id

		my $subs = $fields;

		for my $sub (@$subs) {
			$substitution_count_matrix->{ substr($sub, 0, 1)
				  . substr($sub, -1, 1)}++;
		}
	}
	return $substitution_count_matrix;
}

sub SubstitutionCountMatrixToString {
	my ($substitution_count_matrix) = @_;
	my $subs = [];
	for my $char1 (split "", "ACDEFGHIKLMNPQRSTVWY"){
		for my $char2 (split "", "ACDEFGHIKLMNPQRSTVWY"){
			push @$subs, $char1.$char2;
		}
	}
	return join("\t", map {$_//0} @$substitution_count_matrix{@$subs});
}
