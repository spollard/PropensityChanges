use strict;
use warnings;
use diagnostics;

use Data::Dumper;

our $PerlModulesDir;
BEGIN { require "PerlModulesDir.pl"; }
use lib $PerlModulesDir;

require "plot_hiddenstates.pl";

require "amino_acids.pl";
our $amino_acids;

unless (caller) {
    my $tree_filename = $ARGV[0] // "../data/test_tree3.newick";
    my $sequences_dir = $ARGV[1] // "../data/acceptable/";
    my $phylogram_dir = $ARGV[2] // "../data/acceptable/";
    my $image_width = $ARGV[3] // 400;
    my $image_height = $ARGV[4] // 300;

    my $color_scheme_filename = $ARGV[5] // "../data/Acceptables.aa_color_scheme";
    
    for my $aa (@$amino_acids){
        my $sequences_filename = $sequences_dir . $aa . ".fasta";
        
        my $sequences = Sequences::FromFasta($sequences_filename);
        my $sites = length($sequences->{(keys %$sequences)[0]});
        
        for my $site (1..$sites){
            my $phylogram_filename = $phylogram_dir . $aa . $site . ".png";
            plot_hiddenstates($tree_filename, $sequences_filename,$site,$phylogram_filename,$image_width,$image_height,$color_scheme_filename);
        }
    }
}
