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

require "Site_ToPhylogram.pl";

my $sub_circle_size = 20;

sub SiteToPNG_Figure {
	my ($tree_filename, $sequences_filename, $site, $color_scheme_filename,
		$site_phylogram_filename, $image_width,$image_height, $highlight)
	  = @_;

	my $tree = Tree::FromFile($tree_filename);

	my $color_scheme =ColorScheme::FromFile($color_scheme_filename);
	ColorScheme::SetDefaultColorForResidues($color_scheme, [0, 0, 0]);

	my $sequences = Sequences::FromFasta($sequences_filename);

	Tree::AttachSequences($tree, $sequences);

	CalculateWidths($tree);

	my ($x_expansion, $y_expansion) =
	  CalculateXYExpansion($tree, $image_width, $image_height);

	ResizeTree($tree, $x_expansion, $y_expansion);

	my $img = GD::Simple->new($image_width, $image_height);
	$img->penSize(3,3);
	DrawSiteSubtree_Figure($img, $tree, $tree->{l_width}, 0, $site, $color_scheme, $highlight);

	SaveImage($img, $site_phylogram_filename);
}

sub DrawSiteSubtree_Figure {
	my ($img, $tree, $x_middle, $y_top, $site, $color_scheme, $highlight) = @_;
	my $residue = substr $tree->{Sequence}, $site - 1, 1;

	#print $tree->{Name}, ": ",  $residue;

	$img->fgcolor(@{$color_scheme->{$residue}});
	
	# for site 1070
	#if ($tree->{Name} eq "Node_172" or $tree->{Name} eq "Node_262") {
	# for site 1070 p60
#	if ($tree->{Name} eq "Node_185" or $tree->{Name} eq "Node_263") {
	# for site 1311
	#if ($tree->{Name} eq "Node_1195" or $tree->{Name} eq "Node_972") {
	if ($highlight->{Nodes}->{$tree->{Name}}){
		print "Highlighting $tree->{Name}\n"; 
		$img->moveTo($x_middle, $y_top);
		$img->string("  " . $tree->{Name});
	}
	
	my $y_bottom = $y_top + $tree->{Distance} / 2;
	$img->line($x_middle, $y_top, $x_middle, $y_bottom);



	if ($tree->{Left} and $tree->{Right}) {

		my $x_left = $x_middle - $tree->{l_bar_width};
		my $x_right = $x_middle + $tree->{r_bar_width};

		$img->line($x_left, $y_bottom, $x_right, $y_bottom);

		# if you want to draw angled branches instead...
		# These aren't the kind of angled branches I want
		#		$img->line($x_middle, $y_top, $x_left, $y_bottom);
		#		$img->line($x_middle, $y_top, $x_right, $y_bottom);
		my $y_left = $y_bottom + $tree->{Left}->{Distance} / 2;
		my $y_right = $y_bottom + $tree->{Right}->{Distance} / 2;

		$img->line($x_left, $y_bottom, $x_left, $y_left);
		$img->line($x_right, $y_bottom, $x_right, $y_right);

		# if you see a substitution
		my $l_residue = substr $tree->{Left}->{Sequence}, $site-1, 1;
		if ( $l_residue ne $residue and $l_residue ne "?") {

			$img->moveTo($x_left, $y_left);
			$img->fgcolor(@{$color_scheme->{$l_residue}});
			my $old_bg = $img->bgcolor(@{$color_scheme->{$l_residue}});
			$img->ellipse($sub_circle_size,$sub_circle_size);
			$img->bgcolor($old_bg);
#			
#			$img->fgcolor(@{$color_scheme->{$l_residue}});
			
			#$img->moveTo($x_left-2, $y_left);
			#$img->string($l_residue);
		}

		my $r_residue = substr $tree->{Right}->{Sequence}, $site-1, 1;
		if ( $r_residue ne $residue and $r_residue ne "?") {
			
			# Draw a large filled in circle above here
			$img->moveTo($x_right, $y_right);
			$img->fgcolor(@{$color_scheme->{$r_residue}});
			my $old_bg = $img->bgcolor(@{$color_scheme->{$r_residue}});
			$img->ellipse($sub_circle_size,$sub_circle_size);
			$img->bgcolor($old_bg);
#			
#			$img->fgcolor(@{$color_scheme->{$r_residue}});
			#$img->moveTo($x_right-2, $y_right);

			#$img->string($r_residue);
		}

		DrawSiteSubtree_Figure($img, $tree->{Left}, $x_left, $y_left, $site,
			$color_scheme, $highlight);
		DrawSiteSubtree_Figure($img, $tree->{Right}, $x_right, $y_right, $site,
			$color_scheme, $highlight);
	}
}

