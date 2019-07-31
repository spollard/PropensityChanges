# Read a tree file, draw a phylogram using GD::Simple
# color subtree

use strict;
use warnings;

use lib 'C:\Users\Stephen\Documents\Lab work\PerlModules';


require "UniqueDirectory.pl";
require "SequentialCommands.pl";

#my $tree_filename = "../data/treeoutfile_averaged.newick";
#my $sequences_filename = "../data/sequences.fasta";

my $out_dir = "../data/test/3/";
mkdir $out_dir if not -d $out_dir;

my $tree_filename = "../data/test_tree3.newick";
my $sequences_filename = "../data/test_seqs3.fasta";
my $tree_out_filename = $out_dir . "tree.newick";
my $probs = $out_dir . "probs.txt";
my $rounded = $out_dir . "rounded_hidden_states.txt";
my $png = $out_dir . "most_likely_hiddenstates.png";
my $res_png = $out_dir . "tree_residues.png";

my $width = 1600;
my $height = 1000;

my $commands = [
    # The "1" in the next line is for which site to select from the 
    # sequences file
    "perl FixedCovarionMCMC.pl 2 $tree_filename $sequences_filename $tree_out_filename",
    "Rscript main.R > $probs",
    "perl prob_to_fasta.pl $probs $rounded",
    "perl plot_hiddenstates.pl $tree_out_filename $rounded $png $width $height",
    "perl plot_tree_residues.pl $tree_filename $sequences_filename $res_png $width $height",
    
    #"open $png",
];

SequentialCommands($commands);


