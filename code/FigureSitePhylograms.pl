# Read a tree file, draw a phylogram using GD::Simple
# color subtree

use strict;
use warnings;

use GD::Simple;
use Data::Dumper;

# Can make an SVG instead
#GD::Simple->class('GD::SVG');

use Tree;
use Sequences;
use ColorScheme;

require "Site_toPhylogram_AminoAcidClassFigures.pl";

unless (caller) {
	my $data_dir ="../data/";

	my $tree_filename = $ARGV[0] // $data_dir . "treeoutfile_averaged.newick";
	my $color_scheme_filename = $ARGV[1] // $data_dir . "Stokes.aa_color_scheme";
	
	my $sequences_filename = $ARGV[2] // $data_dir . "sequencesV1-1.fasta";
	#my $sequences_filename = $ARGV[2] // $data_dir . "hiddenStatesV1-1.fasta";
	my $site_phylogram_filename = $ARGV[3] // "../results/s100-1.png";
	
	my $image_width = $ARGV[4] // 3000;
	my $image_height = $ARGV[5] // 1000;
	
	my $site = 1;
	my $highlight = {};
	
	SiteToPNG_Figure($tree_filename, $sequences_filename, $site,
		$color_scheme_filename,$site_phylogram_filename, $image_width,
		$image_height, $highlight);
	
}
