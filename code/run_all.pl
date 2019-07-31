# Read a tree file, draw a phylogram using GD::Simple
# color subtree

use strict;
use warnings;

our $PerlModulesDir;
BEGIN { require "PerlModulesDir.pl"; }
use lib $PerlModulesDir;


require "UniqueDirectory.pl";
require "SequentialCommands.pl";
require "MultithreadCommands.pl";



# Default options
my $options = Options::Read("AcceptableMCMC.options");

# Optionally provide options files on the command line. Later files overwrite newer files
for my $file (@$ARGV){
    %$options = (%$options, Options::Read($file));
}



# my $tree_filename = "../data/treeoutfile_averaged.newick";
# my $sequences_filename = "../data/sequences.fasta";

#my $out_dir = "../data/mlsites_site30_g1000_selected_fixed_swp0.1/";
# my $out_dir = "../data/acceptable_31/";
my $out_dir = "../data/simulated_correct_tightness/";
mkdir $out_dir if not -d $out_dir;

my $generations = 100;

# my $tree_filename = "../data/treeoutfile_averaged.newick";
#my $tree_filename = "../data/node_99.newick";
my $tree_filename = "../data/treeoutfile_new";
# my $sequences_filename = "../data/seqml_400k.fasta";
my $sequences_filename = "../data/simulated/sequences.fasta";

my $acceptabilities_filename = $out_dir . "acceptabilities.fasta";

# my $tree_filename = "../data/test_tree4.newick";
# my $sequences_filename = "../data/test_seqs8.fasta";

my $tree_out_filename = $out_dir . "tree.newick";


my $Nsites = 100;

require "amino_acids.pl";
our $amino_acids;


my $width = 200;
my $height = 200;

my $big_width = 2 * $width;
my $big_height = 2 * $height;

#AcceptableMCMC.pl $tree_filename $sequences_filename $tree_out_filename",

my $model_filename = $out_dir."/model";
my $sites_dir = $out_dir."/sites/";
my $summary_filename = $out_dir."/summary.txt";

my $images_dir = $out_dir . "/images/";
mkdir $images_dir if not -d $images_dir;

my $commands1 = [
    "perl AcceptableMCMC.pl $tree_filename $sequences_filename $generations $out_dir $tree_out_filename",
    "perl stats.pl $model_filename $sites_dir $summary_filename",
    "perl collect_probable_acceptabilities.pl $sites_dir $acceptabilities_filename",
    # "perl prob_to_aa_fastas.pl $sites_dir $out_dir",
    # "perl plot_acceptables.pl $tree_out_filename $out_dir $images_dir $width $height",
    #"perl plot_tree_residues.pl $tree_out_filename $tree_out_filename.fasta $sites_dir $width $height",
    ];
my $commands2 = [
    map {"magick montage \"$images_dir/$_%01d.png[1-$Nsites]\" -tile 1x$Nsites -geometry $width"."x$height \"$images_dir/$_.png\""} @$amino_acids
];
my $commands3= [
    "magick convert \"$images_dir/A.png\" -gravity northwest -strokewidth 2 -pointsize 50 -annotate +0+0 A \"$images_dir/A.png\"",
    "magick convert \"$images_dir/B.png\" -gravity northwest -strokewidth 2 -pointsize 50 -annotate +0+0 B \"$images_dir/B.png\"",
    "magick convert \"$images_dir/C.png\" -gravity northwest -strokewidth 2 -pointsize 50 -annotate +0+0 C \"$images_dir/C.png\"",
    "magick convert \"$images_dir/D.png\" -gravity northwest -strokewidth 2 -pointsize 50 -annotate +0+0 D \"$images_dir/D.png\"",
    "magick convert \"$images_dir/E.png\" -gravity northwest -strokewidth 2 -pointsize 50 -annotate +0+0 E \"$images_dir/E.png\"",
];
my $commands4 = [
    "magick montage $sites_dir/[1-5].png -tile 1x5 -geometry $width"."x$height $out_dir/All_sites.png",
    "magick montage $out_dir/All_sites.png $images_dir/[A-E].png -tile 6x1 -geometry $width"."x$big_height \"$out_dir/All_sites_and_acceptables.png\"",
    
];

SequentialCommands($commands1);
#MultithreadCommands($commands2);
#MultithreadCommands($commands3);
#SequentialCommands($commands4);

