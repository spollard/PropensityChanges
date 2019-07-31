# Find the average tree (same topology) in a treeoutfile after burnin

use strict;
use warnings;
use FindBin;

package Treeoutfile;

use lib "$FindBin::Bin/Treeoutfile";
use Treeoutfile;
use TreeoutfileToNodes;
use TreeoutfileToSingleNewickTree;

unless (caller) {
	print "This package brings together the inferface for how to work with 
	 treeoutfiles\n";
	 print(Treeoutfile::IsOld("../Data/treeoutfile_old"));
}

1;