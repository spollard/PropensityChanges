# Color schemes for alignments of nucleotides or amino acids
# http://www.bioinformatics.nl/~berndb/aacolour.html

package ColorScheme;

use strict;
use warnings;

use Data::Dumper;



unless (caller) {
	my $color_scheme_file = "AColorScheme.data";
	my $colors = FromFile($color_scheme_file);
	print Dumper($colors);
}


sub FromFile {
	my ($color_scheme_file) = @_;
	open my $fh, "<", $color_scheme_file or die $!;
	
	my $header;
	while (($header = <$fh>) =~ /#/) {};
	chomp $header;
    $header =~ s/\r//g;
	my $expected_header = "Residue\tRed\tGreen\tBlue";
	if ($header ne $expected_header) {
		warn "Color scheme file has the wrong header\n";
		warn "Expected: $expected_header\n";
		warn "Found: $header\n";
	}

	my $colors = {};
	while (<$fh>) {
		my ($character, $R, $G, $B) = split;
		$colors->{$character} = [$R,$G,$B];
	}
	return $colors;
}

sub SetDefaultColorForResidues {
	my ($color_scheme, $default_color, $residues) = @_;
	
	$residues = [split //, "ABCDEFGHIJKLMNOPQRSTUVWXYZ.-?0123456789"] if not $residues; 
	
	for my $residue (@$residues) {
		$color_scheme->{$residue} = $default_color if not $color_scheme->{$residue}; 
	}	
}

1;