# Takes all substitutions from a suboutfile and calculates their associations.
# That is the joint probability of all subs with all other subs.
# I an get the total number of subs at a site from a different script.

use strict;
use warnings;

use Data::Dumper;

my $suboutfile = "../Data/pipeline/suboutfile_highprob_ASRV_updatebls7_new_duplicated_highprob";
my $substitution_frequencies_filename = "../Data/substitution_associations_500";

my $substitution_associations = CalculateSubstitutionAssociations($suboutfile);
PrintSubstitutionAssociationFrequencies($substitution_associations,
	$substitution_frequencies_filename);

sub CalculateSubstitutionAssociations {
	my ($suboutfile) = @_;
	open my $fh, $suboutfile or die $!;
	<$fh>; # Burn header file

	my $substitution_associations = {};

	while (my $line = <$fh>) {
		chomp $line;

		my $fields = [split " ", $line];
		shift @$fields; # get rid of tag/branch id/node id

		my $subs = $fields;

		next if (@$subs < 2);
		$subs = [map {substr $_, 1, -1} @$subs]; 

		# This is N logN time instead of N^2 time
		foreach my $sub_1 (0 .. ($#$subs-1)) {
			my $site_1 = $subs->[$sub_1];
			foreach my $sub_2 (($sub_1+1)..$#$subs) {
				$substitution_associations->{$site_1 . "_"
					  . $subs->[$sub_2]}++;
			}
		}
	}
	return $substitution_associations;
}

sub PrintSubstitutionAssociationFrequencies {
	my ($substitution_associations, $file) = @_;

	open my $fh, ">", $file;
	print $fh "Site_1\tSite_2\tFrequency\n";
	while (my ($sites, $frequency) = each %$substitution_associations){
		my ($site_1, $site_2) = split("_", $sites);
		next if ($site_1 > 500 or $site_2 > 500);
		print $fh join("\t", split("_", $sites), $frequency), "\n";
	}
}
