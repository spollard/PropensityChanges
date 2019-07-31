use strict;
use warnings FATAL => 'all';

use Data::Dumper;

our $PerlModulesDir;
BEGIN { require "PerlModulesDir.pl"; }
use lib $PerlModulesDir;

require "plot_hiddenstates.pl";

unless (caller) {
    my $tree_filename = $ARGV[0] // "../data/test_tree3.newick";
    my $sequences_filename = $ARGV[1] // "../data/test_seqs6.fasta";
    my $phylogram_dir = $ARGV[2] // "../data/acceptable14nogibbs/sites/";
    my $image_width = $ARGV[3] // 400;
    my $image_height = $ARGV[4] // 300;
    my $color_scheme_filename = $ARGV[5] // "../data/Acceptables.aa_color_scheme";
    
    mkdir $phylogram_dir if not -d $phylogram_dir;
        
    my $sequences = Sequences::FromFasta($sequences_filename);
    my $sites = length($sequences->{(keys %$sequences)[0]});
    
    for my $site (1..$sites){
        my $phylogram_filename = $phylogram_dir . $site . ".png";
        plot_hiddenstates($tree_filename, $sequences_filename,$site,$phylogram_filename,$image_width,$image_height,$color_scheme_filename);
       
        my $gravity = "northwest";
        my $width = "2";
        my $fontsize = "50";
        my $where = "+0+0";

        system("magick convert $phylogram_filename -gravity $gravity -strokewidth $width -pointsize $fontsize -annotate $where $site $phylogram_filename");

    }
}
