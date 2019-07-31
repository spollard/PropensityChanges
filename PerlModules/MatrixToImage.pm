# Turn a matrix into an image

 package MatrixToImage;

use strict;
use warnings;

use GD::Image;
use ColorScheme;

unless (caller) {
	my $png_file = "default_matrix.png";
	my $color_scheme_file = "default.nuc_color_scheme";

	my $matrix = [
		"ACGT", 
		["F", "A", "G", "G", "T"],
	];

	MatrixToPNG($matrix, $png_file, $color_scheme_file);
}

# This works for both arrays of arrays and arrays of strings
sub MatrixToPNG {
	my ($matrix, $png_file, $colors) = @_;
	
	my ($width, $height) = DetermineSize($matrix);

	my $img = GD::Image->new($width, $height);

	# When using GD::Image, you must allocate space for each color
	# colorResolve() allocates if the color has not been allocated yet. 
	$img->colorResolve(255, 255, 255); # This sets white as the default color
	for (values %$colors) {
		$_ = $img->colorResolve(@$_);
	}

	my $residues = [keys $colors];

	# This should not destroy the matrix
	while (my ($i, $row) = each @$matrix) {
		if (ref $row eq "ARRAY") {
			while (my ($j, $char) = each @$row) {
				if ($char ~~ @$residues){ # This saves time on alignments with lots of gaps
					$img->setPixel($j ,$i, $colors->{$char});
				}
			}
		}
		else {
			for my $j (0 .. length($row) - 1) {
				my $char = substr $row, $j, 1;
				if ($char ~~ @$residues){ # This saves time on alignments with lots of gaps
					$img->setPixel($j ,$i, $colors->{$char});
				}
			}
		}
	}

	open my $png, ">", $png_file;
	binmode $png; 
	print $png $img->png;
}

# Works for both arrays of arrays and arrays of sequences
sub DetermineSize {
	my ($matrix) = @_;
	
	my $height = @$matrix;
	my $width = 0;
	for my $row (@$matrix) {
		if (ref $row eq "ARRAY") {
			$width = @$row if @$row > $width;
		}
		else {
			$width = length $row if length $row > $width;
		}
	}
	
	return $width, $height
}