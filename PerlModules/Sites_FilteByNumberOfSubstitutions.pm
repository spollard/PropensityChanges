# Convert an old suboutfile with branch labels (e.g. ##..##) to new suboutfiles
# with node labels (e.g. Node_136 or Canis_lupus_8234).

use strict;
use warnings;

package Sites;


unless (caller) {
	my $sites_filename = "";
	my $lower_bound = 2;
	my $upper_bound = 6;
	my $sites_output_filename = "";
	
	FilterByNumberOfSubstitutions($sites_filename, $lower_bound, $upper_bound, $sites_output_filename);
	
}


sub FilterByNumberOfSubstitutions {
	my ($sites_filename, $lower_bound, $upper_bound, $sites_output_filename) = @_;
	
	open my $fh, "<", $sites_filename or die $!;
	
	open my $fh_out, ">", $sites_output_filename; 
	
	scalar <$fh>; # Ignore header line
	
	while (<$fh>) {
		my $fields = [split];
		print $fh_out $fields->[0],"\n" if ($fields->[1] >= $lower_bound and $fields->[1] <= $upper_bound); 	
	}	
	
}