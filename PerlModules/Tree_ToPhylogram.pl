# Read a tree file, draw a phylogram using GD::Simple

use strict;
use warnings FATAL => 'all';


use GD::Simple;
use Data::Dumper;

# Can make an SVG instead
#GD::Simple->class('GD::SVG');


use Tree;
use Sequences;
use ColorScheme;


unless (caller) {
    #my $tree_file = $ARGV[0] // "../Data/5taxon/treeoutfile";
    my $tree_file = $ARGV[0] // "../data/1257.tree";
    my $phylogram_file = $ARGV[1] // "../data/simple_tree";
    my $image_width = $ARGV[2] // 8000;
    my $image_height = $ARGV[3] // 6000;
    
    my $sequences_file = "../Data/simple.seqoutfile";
    my $color_scheme_file = "TS.aa_color_scheme";

    my $tree = Tree::FromFile($tree_file);

    TreeToPNG($tree, $phylogram_file, $image_width, $image_height);
    exit;
    
    my $sequences = Sequences::FromFasta($sequences_file);
    my $color_scheme = ColorScheme::FromFile($color_scheme_file);
    ColorScheme::SetDefaultColorForResidues($color_scheme, [0, 0, 0]);

    Tree::AttachSequences($tree, $sequences);
    
    for my $site (1 .. length $tree->{Sequence}){
        print $site . "\n";
        SiteToPNG($tree, $site, $color_scheme, $phylogram_file . $site . ".png", $image_width,
            $image_height);
    }
}

sub SiteToPNG {
    my ($tree, $site, $color_scheme, $phylogram_file, $image_width,
        $image_height) = @_;

    Tree::recurse_post($tree, \&CalculateWidths);
    
    my ($x_expansion, $y_expansion) =
      CalculateXYExpansion($tree, $image_width, $image_height);
      
    Tree::recurse_pre($tree, \&ResizeTree, $x_expansion, $y_expansion);

    my $img = GD::Simple->new($image_width, $image_height);
    
    DrawSiteSubtree($img, $tree, $tree->{l_width}, 0, $site, $color_scheme);

    SaveImage($img, $phylogram_file);
}

sub TreeToPNG {
    my ($tree, $phylogram_file, $image_width, $image_height) = @_;

    Tree::recurse_post($tree, \&CalculateWidths);

    my ($x_expansion, $y_expansion) =
      CalculateXYExpansion($tree, $image_width, $image_height);

    Tree::recurse_pre($tree, \&ResizeTree, $x_expansion, $y_expansion);
    
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
}

# Cannot calculate exact X and Y locations for every node recursively because
# the x location depends on the rest of the tree and not just the subtree.
sub CalculateWidths {
    my ($tree) = @_;
    #print "$tree->{Name}\n";

    if ($tree->{Left} and $tree->{Right}) {
        $tree->{l_width} = $tree->{Left}->{l_width} + $tree->{Left}->{r_width};
        $tree->{r_width} = $tree->{Right}->{l_width} + $tree->{Right}->{r_width};

        $tree->{l_bar_width} = $tree->{Left}->{r_width};
        $tree->{r_bar_width} = $tree->{Right}->{l_width};
    }
    elsif ($tree->{Left} and not $tree->{Right}) {
        $tree->{l_width} = $tree->{Left}->{l_width};
        $tree->{r_width} = $tree->{Left}->{r_width};
        $tree->{l_bar_width} = 0;
        $tree->{r_bar_width} = 0;
    }
    else {
        $tree->{l_width} = 1; # this number must be positive
        $tree->{r_width} = 1; # this number must be positive
        $tree->{l_bar_width} = 0;
        $tree->{r_bar_width} = 0;
    }
}

sub DrawSubtree {
    my ($img, $tree, $x_middle, $y_top) = @_;

    my $y_bottom = $y_top + $tree->{Distance};
    
    $img->moveTo($x_middle, $y_top);
    $img->font('Times');
    $img->fontsize(50);
    $img->string($tree->{Name});
    
    $img->line($x_middle, $y_top, $x_middle, $y_bottom);

    if ($tree->{Left} and $tree->{Right}) {
        my $x_left = $x_middle - $tree->{l_bar_width};
        my $x_right = $x_middle + $tree->{r_bar_width};

        $img->line($x_left, $y_bottom, $x_right, $y_bottom);

        # if you want to draw angled branches instead...
        # These aren't the kind of angled branches I want
        #       $img->line($x_middle, $y_top, $x_left, $y_bottom);
        #       $img->line($x_middle, $y_top, $x_right, $y_bottom);

        DrawSubtree($img, $tree->{Left}, $x_left, $y_bottom);
        DrawSubtree($img, $tree->{Right}, $x_right, $y_bottom);
    }
    elsif ($tree->{Left} and not $tree->{Right}) {
        my $y_bottom = $y_top + $tree->{Distance};
        $img->line($x_middle, $y_top, $x_middle, $y_bottom);
        
        DrawSubtree($img, $tree->{Left}, $x_middle, $y_bottom);
    }
}

sub DrawSiteSubtree {
    my ($img, $tree, $x_middle, $y_top, $site, $color_scheme) = @_;
    my $residue = substr $tree->{Sequence}, $site - 1, 1;
    
    $img->fgcolor(@{$color_scheme->{$residue}});
     
    
    
    my $y_bottom = $y_top + $tree->{Distance} / 2;
    $img->line($x_middle, $y_top, $x_middle, $y_bottom);

    if ($tree->{Left} and $tree->{Right}) {
        
        my $x_left = $x_middle - $tree->{l_bar_width};
        my $x_right = $x_middle + $tree->{r_bar_width};

        $img->line($x_left, $y_bottom, $x_right, $y_bottom);

        # if you want to draw angled branches instead...
        # These aren't the kind of angled branches I want
        #       $img->line($x_middle, $y_top, $x_left, $y_bottom);
        #       $img->line($x_middle, $y_top, $x_right, $y_bottom);
        my $y_left = $y_bottom + $tree->{Left}->{Distance} / 2; 
        my $y_right = $y_bottom + $tree->{Right}->{Distance} / 2; 
        
        $img->line($x_left, $y_bottom, $x_left, $y_left);
        $img->line($x_right, $y_bottom, $x_right, $y_right); 

        # if you see a substitution
        my $l_residue = substr $tree->{Left}->{Sequence}, $site-1, 1;
        if ( $l_residue ne $residue) {
            $img->fgcolor(@{$color_scheme->{$l_residue}});
            $img->moveTo($x_left, $y_left);
            #$img->string("$tree->{Left}->{Name}:$l_residue");
        }
        
        my $r_residue = substr $tree->{Right}->{Sequence}, $site-1, 1;
        if ( $r_residue ne $residue) {
            $img->fgcolor(@{$color_scheme->{$r_residue}});
            $img->moveTo($x_right, $y_right);
            #$img->string("$tree->{Right}->{Name}:$r_residue");
        }
    


        DrawSiteSubtree($img, $tree->{Left}, $x_left, $y_left, $site, $color_scheme);
        DrawSiteSubtree($img, $tree->{Right}, $x_right, $y_right, $site, $color_scheme);
    }
    elsif ($tree->{Left} and not $tree->{Right}) {
        my $y_bottom = $y_top + $tree->{Distance} / 2 + $tree->{Left}->{Distance} / 2;
        $img->line($x_middle, $y_top, $x_middle, $y_bottom);
        
        my $l_residue = substr $tree->{Left}->{Sequence}, $site-1, 1;
        if ( $l_residue ne $residue) {
            $img->fgcolor(@{$color_scheme->{$l_residue}});
            $img->moveTo($x_middle, $y_bottom);
            #$img->string("$tree->{Left}->{Name}:$l_residue");
        }
        
        DrawSiteSubtree($img, $tree->{Left}, $x_middle, $y_bottom, $site, $color_scheme);
    }
}

sub CalculateXYExpansion {
    my ($tree, $image_width, $image_height) = @_;

    my $tree_width = $tree->{l_width} + $tree->{r_width};
    my $tree_height = Tree::FindDeepestDescendant($tree);

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
