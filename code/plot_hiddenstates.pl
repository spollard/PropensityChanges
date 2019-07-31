use strict;
use warnings;
use diagnostics;

use Data::Dumper;

use lib 'C:\Users\Stephen\Documents\Lab work\PerlModules';

use Tree;
use Sequences;
use ColorScheme;

require "Tree_ToPhylogram.pl";

unless (caller) {
    my $tree_filename = $ARGV[0] // "../data/test_tree3.newick";
    my $sequences_filename = $ARGV[1] // "rounded_hidden_states.txt";
    my $phylogram_filename = $ARGV[2] // "../data/test/most_likely_hiddenstates.png";
    my $image_width = $ARGV[3] // 800;
    my $image_height = $ARGV[4] // 600;

    my $color_scheme_filename = $ARGV[5] // "../data/hiddenstate.color_scheme";

    plot_hiddenstates($tree_filename, $sequences_filename,1,$phylogram_filename,$image_width,$image_height,$color_scheme_filename);
}

sub plot_hiddenstates {
    my ($tree_filename, $sequences_filename,$site,$phylogram_filename,$image_width,$image_height,$color_scheme_filename) = @_;
    my $color_scheme = ColorScheme::FromFile($color_scheme_filename);
    ColorScheme::SetDefaultColorForResidues($color_scheme, [0,0,0]);

    my $tree = Tree::FromFile($tree_filename);
    my $sequences = Sequences::FromFasta($sequences_filename);
    Tree::AttachSequences($tree, $sequences);

    SiteToPNG($tree, $site, $color_scheme, $phylogram_filename, $image_width, $image_height);
}