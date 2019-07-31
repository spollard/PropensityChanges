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

# This is for subtrees. 
my $subtree_colors = {
	"Node_838" => [150, 50, 0],
	"Node_410" => [0, 75, 0]
};

# This is for branches. 
my $branch_colors = {
	"Node_838" => [204, 111, 0],
	"Node_410" => [0, 138, 0], 
	"Node_1095" => [0, 0, 200],
	"Node_1026" =>[0, 0, 200]
};

unless (caller) {
	#my $tree_file = $ARGV[0] // "../Data/5taxon/treeoutfile";
	my $tree_file = $ARGV[0] // "../Data/PLEX_data_used_for_paper/treeoutfile_averaged_900";
	my $phylogram_file = $ARGV[1] // "../Data/PLEX_data_used_for_paper/reconstructed_site";
	my $image_width = $ARGV[2] // 2400;
	my $image_height = $ARGV[3] // 1200;
	
	my $sequences_file = "../Data/reconstructed_sequences";
#	$sequences_file ="../Data/PLEX_data_used_for_paper/seq1001";
	
	my $color_scheme_file = "default.aa_color_scheme";
	$color_scheme_file = "Stokes_SIMA.aa_color_scheme";

	my $tree = Tree::FromFile($tree_file);

#	TreeToPNG($tree, $phylogram_file, $image_width, $image_height);
#	exit;
	
	my $color_scheme =
	  ColorScheme::FromFile($color_scheme_file);
	ColorScheme::SetDefaultColorForResidues($color_scheme, [0, 0, 0]);

	my $sequences = Sequences::FromFasta($sequences_file);

	Tree::AttachSequences($tree, $sequences);

}



sub TreeToPNG {
	my ($tree, $phylogram_file, $image_width, $image_height) = @_;

	CalculateWidths($tree);

	my ($x_expansion, $y_expansion) =
	  CalculateXYExpansion($tree, $image_width, $image_height);

	ResizeTree($tree, $x_expansion, $y_expansion);
	
	my $img = GD::Simple->new($image_width, $image_height);

	DrawSubtree($img, $tree, $tree->{l_width}, 0);

	SaveImage($img, $phylogram_file);
}

sub ResizeTree {
	my ($tree, $x_expansion, $y_expansion) = @_;

	$tree->{Distance} = $tree->{Distance} * $y_expansion;

	$tree->{l_bar_width} = $tree->{l_bar_width} * $x_expansion;
	$tree->{r_bar_width} = $tree->{r_bar_width} * $x_expansion;

	# These aren't used in drawing the tree and so are extra work right now.
	# But this function resizes the entire tree and not bits and pieces of it.
	$tree->{l_width} = $tree->{l_width} * $x_expansion;
	$tree->{r_width} = $tree->{r_width} * $x_expansion;

	if ($tree->{Left} and $tree->{Right}) {
		ResizeTree($tree->{Left}, $x_expansion, $y_expansion);
		ResizeTree($tree->{Right}, $x_expansion, $y_expansion);
	}
}

# Cannot calculate exact X and Y locations for every node recursively because
# the x location depends on the rest of the tree and not just the subtree.
sub CalculateWidths {
	my ($tree) = @_;

	if ($tree->{Left} and $tree->{Right}) {
		CalculateWidths($tree->{Left});
		CalculateWidths($tree->{Right});

		$tree->{l_width} = $tree->{Left}->{l_width} + $tree->{Left}->{r_width};
		$tree->{r_width} =
		  $tree->{Right}->{l_width} + $tree->{Right}->{r_width};

		$tree->{l_bar_width} = $tree->{Left}->{r_width};
		$tree->{r_bar_width} = $tree->{Right}->{l_width};
		
	}
	else {
		$tree->{l_width} = 1; # this number must be positive
		$tree->{r_width} = 1; # this number must be positive
		$tree->{l_bar_width} = 0;
		$tree->{r_bar_width} = 0;
	}

	return $tree->{l_width}, $tree->{r_width};
}

sub DrawSubtree {
	my ($img, $tree, $x_middle, $y_top) = @_;

	my $y_bottom = $y_top + $tree->{Distance};
	
	CheckPen($img, $tree->{Name});
	
	$img->line($x_middle, $y_top, $x_middle, $y_bottom);
	ResetBranchPen($img, $tree->{Name});
	$img->line($x_middle, $y_top, $x_middle, $y_bottom);
	
	if ($tree->{Left} and $tree->{Right}) {
		my $x_left = $x_middle - $tree->{l_bar_width};
		my $x_right = $x_middle + $tree->{r_bar_width};

		$img->line($x_left, $y_bottom, $x_right, $y_bottom);

		# if you want to draw angled branches instead...
		# These aren't the kind of angled branches I want
		#		$img->line($x_middle, $y_top, $x_left, $y_bottom);
		#		$img->line($x_middle, $y_top, $x_right, $y_bottom);

		DrawSubtree($img, $tree->{Left}, $x_left, $y_bottom);
		DrawSubtree($img, $tree->{Right}, $x_right, $y_bottom);
	}
	ResetSubtreePen($img, $tree->{Name});
}

sub CheckPen {
	my ($img, $tree_name) = @_;
	if (defined $subtree_colors->{$tree_name}) {
		print "Found a color for $tree_name\n";
		$img->penSize(10,10);
		
		if (ref $subtree_colors->{$tree_name} eq "ARRAY") {
			$img->fgcolor(@{$subtree_colors->{$tree_name}});
		}
		else {
			$img->fgcolor($subtree_colors->{$tree_name});
		}	
	}
	if (defined $branch_colors->{$tree_name}) {
		print "Found a color for $tree_name\n";
		$img->penSize(10,10);
		
		if (ref $branch_colors->{$tree_name} eq "ARRAY") {
			$img->fgcolor(@{$branch_colors->{$tree_name}});
		}
		else {
			$img->fgcolor($branch_colors->{$tree_name});
		}	
	}
}

sub ResetSubtreePen {
	my ($img, $tree_name) = @_;
	if (defined $subtree_colors->{$tree_name}) {
		$img->fgcolor('black');
		$img->penSize(1,1);
	}
}

sub ResetBranchPen {
	my ($img, $tree_name) = @_;
	if (defined $branch_colors->{$tree_name}) {
		$img->fgcolor('black');
		$img->penSize(1,1);
	}
}


sub CalculateXYExpansion {
	my ($tree, $image_width, $image_height) = @_;

	my $tree_width = $tree->{l_width} + $tree->{r_width};
	my $tree_height = Tree::FindDeepestDescendant($tree);

	print "Tree width $tree_width and height $tree_height\n";
	
	my $x_expansion = $image_width / $tree_width;
	my $y_expansion = $image_height / $tree_height;

	return $x_expansion, $y_expansion;
}

sub SaveImage {
	my ($img, $phylogram_file) = @_;

	open my $phylogram_out, '>', $phylogram_file;
	binmode $phylogram_out;

	print $phylogram_out $img->png();
	#print $phylogram_out $img->svg();
}
