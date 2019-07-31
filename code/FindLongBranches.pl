# Read a tree file, draw a phylogram using GD::Simple
# color subtree

use strict;
use warnings;

use Tree;




unless (caller) {
	my $tree_file = $ARGV[0] // "../data/treeoutfile_averaged.newick" ;
    my $threshold = $ARGV[1] // 0.5;
    
    my $tree = Tree::FromFile($tree_file);
    #my $long_branches = Tree::FindLongBranches($tree, $threshold);
    FindLong($tree);
    #print join "\n", @$long_branches;
}

sub addIfLong{
	my ($tree, $long, $threshold) = @_;
	push @$long, $tree->{Name} if $tree->{Distance} > $threshold; 
}


sub FindLong {
	my ($tree) = @_;
	my $long = [];
	my $threshold = 0.5;
	Tree::recurse_pre($tree, \&addIfLong, $long, $threshold);
	print join "\n", @$long;
}